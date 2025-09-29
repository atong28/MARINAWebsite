import json
import logging
import threading
from typing import Dict, Any, Optional
import pickle
import torch
import torch.nn.functional as F
import os

from src.settings import MARINAArgs
from src.model import MARINA
from src.fp_loader import make_fp_loader
from src.data import collate
from src.ranker import RankingSet

logger = logging.getLogger(__name__)
_model_lock = threading.Lock()
_model: Optional[torch.nn.Module] = None
_args: Optional[MARINAArgs] = None
_fp_loader = None
_rankingset_cache = None  # Cache the rankingset to avoid reloading


def _strip_removed_keys(state_dict: Dict[str, Any]) -> Dict[str, Any]:
    # remove any keys that reference removed embeddings
    bad_keywords = ("elem_embed", "cnt_embed")
    return {k: v for k, v in state_dict.items() if not any(b in k for b in bad_keywords)}


def _try_normalize_prefixes(sd: Dict[str, Any], model: torch.nn.Module) -> Dict[str, Any]:
    # Try common lightning prefixes
    prefixes = ("model.", "marina.", "net.")
    model_keys = set(model.state_dict().keys())
    for p in prefixes:
        if any(k.startswith(p) for k in sd.keys()):
            stripped = { (k[len(p):] if k.startswith(p) else k): v for k, v in sd.items() }
            # if many keys now match model keys, return stripped
            match_count = sum(1 for k in stripped.keys() if k in model_keys)
            if match_count > 0:
                logger.info(f"Prefix '{p}' stripped; {match_count} keys match model")
                return stripped
    return sd


def load_model(ckpt_path: str = "data/best.ckpt", params_path: str = "data/params.json") -> torch.nn.Module:
    """Load model and keep it in module global. Safe to call multiple times."""
    global _model, _args, _fp_loader
    if _model is not None:
        return _model

    logger.info("Loading model params from %s", params_path)
    with open(params_path, "r") as f:
        params = json.load(f)
    _args = MARINAArgs(**params)

    # create fp loader according to args
    _fp_loader = make_fp_loader(_args.fp_type, entropy_out_dim=_args.out_dim)

    logger.info("Instantiating MARINA model")
    model = MARINA(_args, _fp_loader)

    logger.info("Loading checkpoint from %s (cpu)", ckpt_path)
    ckpt = torch.load(ckpt_path, map_location="cpu")
    sd = ckpt.get("state_dict", ckpt)
    sd = _strip_removed_keys(sd)
    sd = _try_normalize_prefixes(sd, model)

    # load (allow missing/unexpected)
    res = model.load_state_dict(sd, strict=False)
    logger.info("Model load_state_dict result: %s", res)

    model.to(torch.device("cpu"))
    model.eval()
    torch.set_grad_enabled(False)

    _model = model
    return _model


## Index-based prediction helpers removed for MVP (dataset not hosted on server)


def predict_from_raw(raw_inputs: Dict[str, Any], k: int = 5):
    """Accepts a raw input dict (same spec as SpectralInputLoader.load output) and returns top-k."""
    logger.info(f"predict_from_raw called with k={k} (type: {type(k)})")
    logger.info(f"raw_inputs keys: {list(raw_inputs.keys())}")
    
    global _model
    if _model is None:
        load_model()

    # expect raw_inputs to be a dict of modality->arrays/tensors matching training shapes
    # we assume caller provides numpy lists; convert to tensors
    processed = {}
    for k_mod, v in raw_inputs.items():
        if isinstance(v, list):
            tensor = torch.tensor(v, dtype=torch.float)
            # Reshape sequence data to 2D (num_peaks, features_per_peak)
            if k_mod in ['hsqc', 'h_nmr', 'c_nmr', 'mass_spec']:
                if k_mod == 'hsqc':
                    # HSQC: [H1, C13, intensity, H1, C13, intensity, ...] -> [C13, H1, intensity, C13, H1, intensity, ...]
                    if len(v) % 3 == 0:
                        tensor = tensor.view(-1, 3)
                        # Reorder from [H, C, intensity] to [C, H, intensity] (canonical order)
                        tensor = tensor[:, [1, 0, 2]]  # Swap columns 0 and 1
                    else:
                        logger.warning(f"HSQC data length {len(v)} is not divisible by 3, skipping")
                        continue
                elif k_mod == 'c_nmr':
                    # C NMR: [shift1, shift2, shift3, ...] -> (N,1) -> (N,3) with zero padding
                    tensor = tensor.view(-1, 1)
                    tensor = F.pad(tensor, (0, 2), "constant", 0)  # -> (N,3)
                elif k_mod == 'h_nmr':
                    # H NMR: [shift1, shift2, shift3, ...] -> (N,1) -> (N,3) with zero padding
                    tensor = tensor.view(-1, 1)
                    tensor = F.pad(tensor, (1, 1), "constant", 0)  # -> (N,3)
                elif k_mod == 'mass_spec':
                    # Mass spec: [mz1, intensity1, mz2, intensity2, ...]
                    if len(v) % 2 == 0:
                        tensor = tensor.view(-1, 2)
                        tensor = F.pad(tensor, (0, 1), "constant", 0)  # -> (N,3)
                    else:
                        logger.warning(f"Mass spec data length {len(v)} is not divisible by 2, skipping")
                        continue
            processed[k_mod] = tensor
        elif isinstance(v, torch.Tensor):
            processed[k_mod] = v
        else:
            try:
                processed[k_mod] = torch.tensor(v, dtype=torch.float)
            except Exception:
                processed[k_mod] = v

    logger.info(f"Processed inputs: {list(processed.keys())}")
    for k_proc, v in processed.items():
        if isinstance(v, torch.Tensor):
            logger.info(f"  {k_proc}: tensor shape {v.shape}, dtype {v.dtype}")
        else:
            logger.info(f"  {k_proc}: {type(v)} = {v}")

    batch_inputs, batch_fps = collate([(processed, torch.zeros((_args.out_dim,), dtype=torch.float))])
    
    logger.info(f"Batch inputs after collate: {list(batch_inputs.keys())}")
    for k_batch, v in batch_inputs.items():
        if isinstance(v, torch.Tensor):
            logger.info(f"  {k_batch}: tensor shape {v.shape}, dtype {v.dtype}")
        else:
            logger.info(f"  {k_batch}: {type(v)} = {v}")

    with _model_lock:
        with torch.no_grad():
            logger.info("Calling model forward...")
            out = _model(batch_inputs)
            logger.info(f"Model output shape: {out.shape}")
            preds = out.detach().cpu()
            pred_fp = preds[0]
            
            # Ensure k is an integer and pred_fp.numel() is an integer
            try:
                k_int = int(k)
            except (ValueError, TypeError) as e:
                logger.error(f"Failed to convert k={k} (type: {type(k)}) to int: {e}")
                k_int = 5  # fallback to default
            
            try:
                numel_int = int(pred_fp.numel())
            except (ValueError, TypeError) as e:
                logger.error(f"Failed to convert pred_fp.numel()={pred_fp.numel()} to int: {e}")
                numel_int = pred_fp.numel()  # fallback to original value
            
            # Use RankingSet.retrieve_idx() for consistency with existing codebase
            # Apply sigmoid to get probabilities
            pred = torch.sigmoid(pred_fp)
            logger.info(f"Pred FP sum: {torch.where(pred.squeeze()>0.5, 1, 0).sum()}")
            
            # Use the existing RankingSet to retrieve results (with caching)
            global _rankingset_cache
            if _rankingset_cache is None:
                _rankingset_cache = _fp_loader.load_rankingset(_args.fp_type)
            ranker = RankingSet(store=_rankingset_cache, metric="cosine")
            
            # Retrieve top-k indices using the existing method
            idxs = ranker.retrieve_idx(pred.unsqueeze(0), n=k_int)
            
            # Get similarity scores for the retrieved indices
            sims = ranker._sims(pred.unsqueeze(0))  # (N, 1)
            sims_sorted, _ = torch.topk(sims.squeeze(), k=k_int, dim=0)
            
            # Convert to lists for return
            vals = sims_sorted.tolist()
            idxs = idxs.squeeze().tolist()

    return vals, idxs, pred.tolist()
