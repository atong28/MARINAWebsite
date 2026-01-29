"""
Domain fingerprint utilities - imports from marina and adds domain-specific functions.
"""
import math
import torch

# Import utility functions directly from marina
from src.marina.src.modules.data.fp_utils import (
    count_circular_substructures,
    get_bitinfos,
    _mk_rdkit,
    BitInfo,
)

def tanimoto_similarity(a: torch.Tensor, b: torch.Tensor) -> float:
    """
    Generalized Tanimoto similarity in [0, 1] for real-valued vectors.

    For binary vectors this reduces to the standard Jaccard/Tanimoto:
      T(a, b) = (a · b) / (||a||^2 + ||b||^2 - a · b)

    Inputs are converted to 1D float tensors, truncated to the shorter length
    if shapes differ, and the result is clamped to [0.0, 1.0].
    """
    if a is None or b is None:
        return 0.0

    a = a.detach().to(dtype=torch.float32).view(-1)
    b = b.detach().to(dtype=torch.float32).view(-1)

    if a.numel() == 0 or b.numel() == 0:
        return 0.0

    if a.numel() != b.numel():
        n = min(a.numel(), b.numel())
        a = a[:n]
        b = b[:n]

    dot = torch.dot(a, b)
    denom = a.pow(2).sum() + b.pow(2).sum() - dot
    if denom <= 0:
        return 0.0

    val = (dot / denom).item()
    if not math.isfinite(val):
        return 0.0

    return max(0.0, min(1.0, float(val)))

__all__ = [
    'count_circular_substructures',
    'get_bitinfos',
    '_mk_rdkit',
    'BitInfo',
    'tanimoto_similarity',
]
