"""
Domain data module - imports MARINA dataset directly from marina.
"""
import torch
from torch.nn.utils.rnn import pad_sequence
from src.marina.src.modules.marina.dataset import MARINADataset, MARINADataModule
from src.marina.src.modules.core.const import INPUTS_CANONICAL_ORDER

def collate(batch):
    """
    Collate function extracted from marina's MARINADataModule._collate_fn.
    """
    dicts, fps = zip(*batch)
    batch_inputs = {}

    for mod in INPUTS_CANONICAL_ORDER:
        seqs = [d.get(mod) for d in dicts]
        if all(x is None for x in seqs):
            continue
        D = next(x.shape[1] for x in seqs if isinstance(
            x, torch.Tensor) and x.ndim == 2)
        seqs = [
            x if (isinstance(x, torch.Tensor) and x.ndim ==
                  2) else torch.zeros((0, D), dtype=torch.float)
            for x in seqs
        ]
        batch_inputs[mod] = pad_sequence(seqs, batch_first=True)

    batch_fps = torch.stack(fps, dim=0)
    return batch_inputs, batch_fps

__all__ = ['MARINADataset', 'MARINADataModule', 'collate']
