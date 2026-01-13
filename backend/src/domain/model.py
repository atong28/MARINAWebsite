import logging
import pytorch_lightning as pl
import torch
import torch.nn as nn
import torch.distributed as dist

from src.domain.const import SELF_ATTN_INPUTS
from src.domain.settings import MARINAArgs
from src.domain.ranker import RankingSet
from src.domain.fingerprint.fp_loader import FPLoader
from src.domain.encoder import build_encoder

logger = logging.getLogger("lightning")
if dist.is_initialized():
    rank = dist.get_rank()
    if rank != 0:
        logger.setLevel(logging.WARNING)
logger_should_sync_dist = torch.cuda.device_count() > 1

class CrossAttentionBlock(nn.Module):
    """
    Query (global CLS) attends to Key/Value (the spectral peaks + other tokens).
    """

    def __init__(self, dim_model, num_heads, ff_dim, dropout=0.1):
        super().__init__()
        self.attn = nn.MultiheadAttention(
            embed_dim=dim_model,
            num_heads=num_heads,
            dropout=dropout,
            bias=True,
            batch_first=True,
        )
        self.norm1 = nn.LayerNorm(dim_model)
        self.ff = nn.Sequential(
            nn.Linear(dim_model, ff_dim),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
            nn.Linear(ff_dim, dim_model),
        )
        self.norm2 = nn.LayerNorm(dim_model)

    def forward(self, query, key, value, key_padding_mask=None):
        attn_out, _ = self.attn(
            query,
            key,
            value,
            key_padding_mask=key_padding_mask,
        )
        q1 = self.norm1(query + attn_out)
        ff_out = self.ff(q1)
        out = self.norm2(q1 + ff_out)
        return out

class MARINA(pl.LightningModule):
    def __init__(self, args: MARINAArgs, fp_loader: FPLoader):
        super().__init__()

        self.args = args
        self.fp_loader = fp_loader
        if self.global_rank == 0:
            logger.info("[MARINA] Started Initializing")
        self.fp_length = args.out_dim
        self.out_dim = args.out_dim
        self.batch_size = args.batch_size
        self.lr = args.lr
        self.weight_decay = args.weight_decay
        self.heads = args.heads
        self.layers = args.layers
        self.ff_dim = args.ff_dim
        self.dropout = args.dropout
        self.scheduler = args.scheduler
        self.dim_model = args.dim_model
        self.use_jaccard = args.use_jaccard
        self.freeze_weights = args.freeze_weights

        self.enc_nmr = build_encoder(
            args.dim_model,
            args.nmr_dim_coords,
            [args.c_wavelength_bounds, args.h_wavelength_bounds],
            args.nmr_is_sign_encoding
        )
        self.enc_c_nmr = build_encoder(
            args.dim_model,
            args.c_nmr_dim_coords,
            [args.c_wavelength_bounds],
            args.c_nmr_is_sign_encoding
        )
        self.enc_h_nmr = build_encoder(
            args.dim_model,
            args.h_nmr_dim_coords,
            [args.h_wavelength_bounds],
            args.h_nmr_is_sign_encoding
        )
        self.enc_ms = build_encoder(
            args.dim_model,
            args.ms_dim_coords,
            [args.mz_wavelength_bounds, args.intensity_wavelength_bounds],
            args.ms_is_sign_encoding
        )
        self.enc_mw = build_encoder(
            args.dim_model,
            args.mw_dim_coords,
            [args.mw_wavelength_bounds],
            args.mw_is_sign_encoding
        )
        self.encoders = {
            "hsqc": self.enc_nmr,
            "h_nmr": self.enc_h_nmr,
            "c_nmr": self.enc_c_nmr,
            "mass_spec": self.enc_ms,
            "mw": self.enc_mw
        }
        self.encoders = nn.ModuleDict(
            {k: v for k, v in self.encoders.items() if k in self.args.input_types})
        self.self_attn = nn.ModuleDict({
            modality: nn.TransformerEncoder(
                nn.TransformerEncoderLayer(
                    d_model=self.dim_model, nhead=self.heads,
                    dim_feedforward=self.ff_dim,
                    batch_first=True, dropout=self.dropout
                ),
                num_layers=args.self_attn_layers[modality]
            )
            for modality in self.encoders
        })
        self.mod_tokens = nn.ParameterDict({
            modality: nn.Parameter(torch.randn(1, 1, self.dim_model))
            for modality in self.encoders
        })
        self.cross_blocks = nn.ModuleList([
            CrossAttentionBlock(
                dim_model=self.dim_model,
                num_heads=self.heads,
                ff_dim=self.ff_dim,
                dropout=self.dropout
            )
            for _ in range(self.layers)
        ])
        self.global_cls = nn.Parameter(torch.randn(1, 1, self.dim_model))
        self.fc = nn.Linear(self.dim_model, self.out_dim)
        self._val_mm = torch.nn.ModuleDict()
        self._test_mm = torch.nn.ModuleDict()
        if self.freeze_weights:
            for parameter in self.parameters():
                parameter.requires_grad = False
        self.ranker = None
        self.spectral_types = ['all_inputs'] + ['_'.join(types) for types in self.args.additional_test_types]
        if self.global_rank == 0:
            logger.info("[MARINA] Initialized")

    def forward(self, batch, batch_idx=None, return_representations=False):
        B = next(iter(batch.values())).size(0)
        all_points = []
        all_masks = []
        for m, x in batch.items():
            if m not in SELF_ATTN_INPUTS:
                continue
            B, L, D_in = x.shape
            mask = (x.abs().sum(-1) == 0)
            enc_seq = self.encoders[m](x.view(B * L, D_in)).view(B, L, self.dim_model)
            mod_token = self.mod_tokens[m].to(enc_seq.device).expand(B, 1, -1)
            enc_seq = torch.cat([mod_token, enc_seq], dim=1)
            mask = torch.cat([torch.zeros(B, 1, dtype=torch.bool, device=enc_seq.device), mask], dim=1)
            attended = self.self_attn[m](enc_seq, src_key_padding_mask=mask)
            all_points.append(attended)
            all_masks.append(mask)
        joint_seq = torch.cat(all_points, dim=1)
        joint_mask = torch.cat(all_masks, dim=1)
        global_token = self.global_cls.expand(B, 1, -1)
        for block in self.cross_blocks:
            global_token = block(
                query=global_token,
                key=joint_seq,
                value=joint_seq,
                key_padding_mask=joint_mask
            )
        out = self.fc(global_token.squeeze(1))
        if return_representations:
            return global_token.squeeze(1).detach().cpu().numpy()
        return out

    def log(self, name, value, *args, **kwargs):
        if kwargs.get('sync_dist') is None:
            kwargs['sync_dist'] = logger_should_sync_dist
        super().log(name, value, *args, **kwargs)

    def setup_ranker(self):
        store = self.fp_loader.load_rankingset(self.args.fp_type)
        metric = "cosine"

        self.ranker = RankingSet(store=store, metric=metric)
