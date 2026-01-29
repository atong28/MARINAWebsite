"""
Domain fingerprint loader module - imports directly from marina.
Adds backend-only factory for model-root–scoped loaders.
"""
import os

from src.marina.src.modules.data.fp_loader import (
    FPLoader,
    EntropyFPLoader,
    make_fp_loader,
)


def make_fp_loader_for_model_root(
    model_root: str,
    fp_type: str = "RankingEntropy",
    entropy_out_dim: int = 16384,
    max_radius: int = 6,
) -> EntropyFPLoader:
    """
    Backend-only factory: create EntropyFPLoader with dataset_root and retrieval
    both under model_root. Use this instead of make_fp_loader when loading from
    a per-model directory (e.g. data/marina_best/).
    """
    retrieval_path = os.path.join(model_root, "retrieval.pkl")
    loader = EntropyFPLoader(dataset_root=model_root, retrieval_path=retrieval_path)
    loader.setup(entropy_out_dim, max_radius, retrieval_path=retrieval_path)
    return loader


__all__ = ["FPLoader", "EntropyFPLoader", "make_fp_loader", "make_fp_loader_for_model_root"]
