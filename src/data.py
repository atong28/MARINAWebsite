import os
import pickle
import torch
import traceback
import sys
import logging
import random
from itertools import islice
from typing import Any, Optional

from torch.utils.data import DataLoader, Dataset
from torch.nn.utils.rnn import pad_sequence
import pytorch_lightning as pl

from src.const import DEBUG_LEN, DROP_PERCENTAGE, INPUTS_CANONICAL_ORDER, DATASET_ROOT
from src.settings import MARINAArgs
from src.fp_loader import FPLoader
from src.inputs import SpectralInputLoader, MFInputLoader
# from src.modality_dropout_scheduler import ModalityDropoutScheduler

logger = logging.getLogger('lightning')

class MARINADataset(Dataset):
    def __init__(self, args: MARINAArgs, fp_loader: FPLoader, split: str = 'train', override_input_types: Optional[list[str]] = None):
        try:
            self.args = args  # <-- store args

            if split != 'train':
                args.requires = args.input_types
            logger.info(f'[MARINADataset] Initializing {split} dataset with input types {args.input_types} and required inputs {args.requires}')
            self.input_types = args.input_types if override_input_types is None else override_input_types
            self.requires = args.requires if override_input_types is None else override_input_types

            with open(os.path.join(DATASET_ROOT, 'index.pkl'), 'rb') as f:
                data: dict[int, Any] = pickle.load(f)
            data = {
                idx: entry for idx, entry in data.items()
                if entry['split'] == split and
                any(
                    entry[f'has_{input_type}']
                    for input_type in self.input_types
                    if input_type not in ('mw', 'formula')
                )
            }
            data_len = len(data)
            logger.info(f'[MARINADataset] Requiring the following items to be present: {self.requires}')
            data = {
                idx: entry for idx, entry in data.items()
                if all(entry[f'has_{dtype}'] for dtype in self.requires)
            }
            logger.info(f'[MARINADataset] Purged {data_len - len(data)}/{data_len} items. {len(data)} items remain')
            print(f'[MARINADataset] Dataset size: {len(data)}')
            if args.debug and len(data) > DEBUG_LEN:
                logger.info(f'[MARINADataset] Debug mode activated. Data length set to {DEBUG_LEN}')
                data = dict(islice(data.items(), DEBUG_LEN))

            if len(data) == 0:
                raise RuntimeError(f'[MARINADataset] Dataset split {split} is empty!')
            
            self.jittering = args.jittering
            self.spectral_loader = SpectralInputLoader(DATASET_ROOT, data)
            self.mfp_loader = MFInputLoader(fp_loader)
            
            self.data = list(data.items())

            # Initialize modality dropout scheduler if requested
            if self.args.modality_dropout_scheduler:
                # ModalityDropoutScheduler is not available, skip for now
                self.drop_scheduler = None
            else:
                self.drop_scheduler = None

            logger.info('[MARINADataset] Setup complete!')
        
        except Exception:
            logger.error(traceback.format_exc())
            logger.error('[MARINADataset] While instantiating the dataset, ran into the above error.')
            sys.exit(1)
    
    def __len__(self):
        return len(self.data)

    def set_phase(self, phase: float):
        """Forward curriculum phase to the modality dropout scheduler, if enabled."""
        if hasattr(self, "drop_scheduler") and self.drop_scheduler is not None:
            if hasattr(self.drop_scheduler, "set_phase"):
                self.drop_scheduler.set_phase(phase)
    
    def __getitem__(self, idx):
        data_idx, data_obj = self.data[idx]
        if self.args.modality_dropout_scheduler and self.drop_scheduler is not None:
            input_types = set(self.drop_scheduler.sample_kept_modalities(
                data_obj, requires=set(self.requires)      # e.g. {'hsqc'} or empty set
            ))
        else:
            available_types = {
                'hsqc': data_obj['has_hsqc'],
                'c_nmr': data_obj['has_c_nmr'],
                'h_nmr': data_obj['has_h_nmr'],
                'mass_spec': data_obj['has_mass_spec']
            }
            drop_candidates = [k for k, v in available_types.items() if k in self.input_types and v]
            assert len(drop_candidates) > 0, 'Found an empty entry!'
            
            always_keep = random.choice(drop_candidates)
            input_types = set(self.input_types)
            for input_type in self.input_types:
                if not data_obj[f'has_{input_type}']:
                    input_types.remove(input_type)
                elif (input_type != always_keep and 
                    input_type not in self.requires and 
                    random.random() < DROP_PERCENTAGE[input_type]):
                    input_types.remove(input_type)
        return self.spectral_loader.load(data_idx, input_types, jittering = self.jittering), self.mfp_loader.load(data_idx)

def collate(batch):
    """
    batch: list of (data_inputs: dict, mfp: Tensor)
    returns: (batch_inputs: dict[str→Tensor], batch_fps: Tensor)
    """
    logger = logging.getLogger(__name__)
    
    logger.info(f"Collate called with batch size: {len(batch)}")
    dicts, fps = zip(*batch)
    
    # Debug the input dictionaries
    for i, d in enumerate(dicts):
        logger.info(f"  Sample {i}: {list(d.keys())}")
        for k, v in d.items():
            if isinstance(v, torch.Tensor):
                logger.info(f"    {k}: tensor shape {v.shape}, dtype {v.dtype}")
                logger.info(str(v)) 
            else:
                logger.info(f"    {k}: {type(v)} = {v}")
    
    batch_inputs = {}

    # 1) Handle all the *sequence* modalities
    for mod in INPUTS_CANONICAL_ORDER:
        if mod == "mw":
            # skip MW here—handle below
            continue

        seqs = [d.get(mod) for d in dicts]
        # if none of the samples have this modality, skip it entirely
        if all(x is None for x in seqs):
            continue

        # replace missing with empty (0×D) tensors
        # find the first real tensor to infer D
        real_tensors = [x for x in seqs if isinstance(x, torch.Tensor) and x.ndim == 2]
        if real_tensors:
            D = real_tensors[0].shape[1]
            seqs = [
                x if (isinstance(x, torch.Tensor) and x.ndim == 2) else torch.zeros((0, D), dtype=torch.float)
                for x in seqs
            ]
        else:
            # If no real tensors found, skip this modality entirely
            continue
        # now pad them into a (B, L_mod, D) tensor
        batch_inputs[mod] = pad_sequence(seqs, batch_first=True)

    # 2) Handle MW *scalar* specially
    mw_vals = [d.get("mw") for d in dicts]
    if any(v is not None for v in mw_vals):
        # replace None with 0.0 (or another sentinel if you like)
        mw_floats = [float(v) if v is not None else 0.0 for v in mw_vals]
        # create a (B,) tensor of scalars
        batch_inputs["mw"] = torch.tensor(mw_floats, dtype=torch.float)

    # 3) Handle element‐group tokens (idx + count)
    elem_idx_seqs = [d.get('elem_idx') for d in dicts]
    if any(x is not None for x in elem_idx_seqs):
        # pad element‐ID sequences (pad_value=0)
        batch_inputs['elem_idx'] = pad_sequence(
            [x if x is not None else torch.zeros(0,dtype=torch.long)
             for x in elem_idx_seqs],
            batch_first=True,
            padding_value=0
        )
        # pad count sequences (pad_value=0)
        cnt_seqs = [d.get('elem_cnt') for d in dicts]
        batch_inputs['elem_cnt'] = pad_sequence(
            [x if x is not None else torch.zeros(0,dtype=torch.long)
             for x in cnt_seqs],
            batch_first=True,
            padding_value=0
        )

    # 4) Stack your fingerprints
    batch_fps = torch.stack(fps, dim=0)
    
    # Debug the final batch
    logger.info(f"Final batch_inputs keys: {list(batch_inputs.keys())}")
    for k, v in batch_inputs.items():
        if isinstance(v, torch.Tensor):
            logger.info(f"  {k}: tensor shape {v.shape}, dtype {v.dtype}")
        else:
            logger.info(f"  {k}: {type(v)} = {v}")
    
    if not batch_inputs:
        logger.warning("WARNING: batch_inputs is empty!")
    
    return batch_inputs, batch_fps

class MARINADataModule(pl.LightningDataModule):
    def __init__(self, args: MARINAArgs, fp_loader: FPLoader):
        super().__init__()
        self.args = args
        self.batch_size = args.batch_size
        self.num_workers = args.num_workers
        self.collate_fn = collate
        self.persistent_workers = bool(args.persistent_workers and self.num_workers > 0)
        self.fp_loader = fp_loader
        if args.hybrid_early_stopping:
            self.test_types = [args.input_types] + [[input_type] for input_type in (set(args.input_types) - {'mw', 'formula'})]
        else:
            self.test_types = [args.input_types]
        
        self._fit_is_setup = False
        self._test_is_setup = False
    
    def setup(self, stage):
        if (stage == "fit" or stage == "validate" or stage is None) and not self._fit_is_setup:
            self.train = MARINADataset(self.args, self.fp_loader, split='train')
            self.val = [MARINADataset(self.args, self.fp_loader, split='val', override_input_types=input_type) for input_type in self.test_types]
            self._fit_is_setup = True
        if (stage == "test") and not self._test_is_setup:
            self.test = [MARINADataset(self.args, self.fp_loader, split='test', override_input_types=input_type) for input_type in self.test_types]
            self._test_is_setup = True
        if stage == "predict":
            raise NotImplementedError("Predict setup not implemented")
    
    def __getitem__(self, idx):
        if not self._fit_is_setup:
            self.setup(stage = 'fit')
        return self.train[idx]
    
    def train_dataloader(self):
        if not self._fit_is_setup:
            self.setup(stage = 'fit')
        return DataLoader(
            self.train,
            shuffle=True,
            batch_size=self.batch_size,
            collate_fn=self.collate_fn,
            num_workers=self.num_workers,
            pin_memory=True, 
            persistent_workers=self.persistent_workers
        )

    def val_dataloader(self):
        if not self._fit_is_setup:
            self.setup(stage = 'fit')
        return [DataLoader(
            val_dl,
            batch_size=self.batch_size,
            collate_fn=self.collate_fn, 
            num_workers=self.num_workers,
            pin_memory=True,
            persistent_workers=self.persistent_workers
        ) for val_dl in self.val]

    def test_dataloader(self):
        if not self._test_is_setup:
            self.setup(stage = 'test')
        return [DataLoader(
            test_dl,
            batch_size=self.batch_size,
            collate_fn=self.collate_fn, 
            num_workers=self.num_workers,
            pin_memory=True,
            persistent_workers=self.persistent_workers
        ) for test_dl in self.test]
