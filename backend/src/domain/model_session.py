import bisect
import json
import logging
import os
import threading
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import torch

from src.config import CKPT_PATH, MODEL_ROOT, PARAMS_PATH, RETRIEVAL_PATH
from src.domain.fingerprint.fp_loader import FPLoader, make_fp_loader_for_model_root
from src.domain.model import MARINA, SPECTRE
from src.domain.ranker import RankingSet
from src.domain.settings import MARINAArgs, SPECTREArgs


logger = logging.getLogger(__name__)


def _load_metadata(metadata_path: str) -> Dict[str, Any]:
    with open(metadata_path, "r") as f:
        return json.load(f)


def _mw_from_entry(entry: Optional[Dict[str, Any]]) -> Optional[float]:
    if not entry:
        return None
    m = entry.get("mw")
    if isinstance(m, (int, float)):
        try:
            return float(m)
        except Exception:
            pass
    return None


@dataclass
class ModelSession:
    """
    Encapsulates the loaded model (MARINA or SPECTRE), fingerprint loader,
    rankingset store, metadata, and MW index for O(1) range filtering.
    """
    # Model metadata
    type: str  # e.g. "marina" or "spectre"

    # Underlying model and configuration
    args: MARINAArgs | SPECTREArgs
    model: torch.nn.Module
    fp_loader: FPLoader
    fp_type: str
    ckpt_path: str
    params_path: str
    retrieval_path: str
    model_root: str
    _metadata: Dict[str, Any]
    _metadata_path: str

    _rankingset_store: Optional[torch.Tensor] = None
    _rankingset_wrapper: Optional[RankingSet] = None
    _filtered_rankingsets: Dict[Tuple[Optional[float], Optional[float]], RankingSet] = field(default_factory=dict, repr=False)
    _mw_sorted: Optional[List[tuple]] = None  # (mw, idx) sorted by mw
    _mw_by_idx: Optional[Dict[int, float]] = None
    _num_rows: Optional[int] = None
    _lock: threading.RLock = field(default_factory=threading.RLock, repr=False)

    @classmethod
    def from_model_root(cls, model_root: str, model_type: Optional[str] = None) -> "ModelSession":
        """
        Load model and metadata from a per-model directory (e.g. data/marina_best/).
        Uses make_fp_loader_for_model_root so dataset_root and retrieval both
        point under model_root.
        """
        ckpt_path = os.path.join(model_root, "best.ckpt")
        params_path = os.path.join(model_root, "params.json")
        retrieval_path = os.path.join(model_root, "retrieval.pkl")
        metadata_path = os.path.join(model_root, "metadata.json")

        # Default to MARINA for backwards compatibility when type is not provided.
        model_type = model_type or "marina"

        logger.debug("Loading model params from %s", params_path)
        with open(params_path, "r") as f:
            params = json.load(f)

        if model_type == "spectre":
            args: MARINAArgs | SPECTREArgs = SPECTREArgs(**params)
        else:
            # Treat any non-spectre type as MARINA-compatible by default.
            args = MARINAArgs(**params)

        fp_loader = make_fp_loader_for_model_root(
            model_root,
            fp_type=getattr(args, "fp_type", "RankingEntropy"),
            entropy_out_dim=args.out_dim,
            max_radius=6,
        )

        if model_type == "spectre":
            logger.debug("Instantiating SPECTRE model")
            model = SPECTRE(args, fp_loader)
        else:
            logger.debug("Instantiating MARINA model")
            model = MARINA(args, fp_loader)

        logger.debug("Loading checkpoint from %s (cpu)", ckpt_path)
        ckpt = torch.load(ckpt_path, map_location="cpu")
        sd = ckpt.get("state_dict", ckpt)
        res = model.load_state_dict(sd, strict=False)
        logger.debug("Model load_state_dict result: %s", res)

        model.to(torch.device("cpu"))
        model.eval()
        torch.set_grad_enabled(False)

        metadata = _load_metadata(metadata_path)

        return cls(
            type=model_type,
            args=args,
            model=model,
            fp_loader=fp_loader,
            fp_type=args.fp_type,
            ckpt_path=ckpt_path,
            params_path=params_path,
            retrieval_path=retrieval_path,
            model_root=model_root,
            _metadata=metadata,
            _metadata_path=metadata_path,
        )

    @classmethod
    def from_paths(
        cls,
        ckpt_path: str = CKPT_PATH,
        params_path: str = PARAMS_PATH,
        retrieval_path: str = RETRIEVAL_PATH,
    ) -> "ModelSession":
        """
        Factory that loads from explicit paths. Derives model_root from
        ckpt_path directory and delegates to from_model_root.
        """
        model_root = os.path.dirname(os.path.abspath(ckpt_path))
        return cls.from_model_root(model_root)

    def _ensure_mw_index(self) -> None:
        with self._lock:
            if self._mw_sorted is not None:
                return
            store = self.fp_loader.load_rankingset(self.fp_type)
            self._rankingset_store = store
            n = store.shape[0]
            self._num_rows = n
            mw_by_idx: Dict[int, float] = {}
            for i in range(n):
                entry = self._metadata.get(str(i))
                mw = _mw_from_entry(entry)
                if mw is not None:
                    mw_by_idx[i] = mw
            self._mw_by_idx = mw_by_idx
            # (mw, idx) sorted by mw for range queries
            self._mw_sorted = sorted((mw, idx) for idx, mw in mw_by_idx.items())

    def get_rankingset(self) -> RankingSet:
        """Lazily load the retrieval rankingset and MW index, then wrap in RankingSet."""
        if self._rankingset_store is None:
            self._ensure_mw_index()
        if self._rankingset_wrapper is None:
            assert self._rankingset_store is not None
            self._rankingset_wrapper = RankingSet(store=self._rankingset_store, metric="cosine")
        return self._rankingset_wrapper

    def get_filtered_rankingset(
        self,
        mw_min: Optional[float],
        mw_max: Optional[float],
    ) -> RankingSet:
        """
        Get a cached RankingSet filtered by MW range.
        Creates and caches filtered RankingSet instances per MW range to avoid
        repeated instantiation on every API call.
        """
        cache_key = (mw_min, mw_max)
        
        with self._lock:
            if cache_key in self._filtered_rankingsets:
                return self._filtered_rankingsets[cache_key]
            
            # Get base rankingset store
            base_ranker = self.get_rankingset()
            base_store = base_ranker.data
            
            # Compute kept indices for MW range
            kept_indices = self.indices_in_mw_range(mw_min, mw_max)
            
            # Filter the store
            from src.domain.ranker import filter_rankingset_rows
            filtered_store = filter_rankingset_rows(base_store, kept_indices)
            
            # Create and cache filtered RankingSet
            filtered_ranker = RankingSet(store=filtered_store, metric="cosine")
            self._filtered_rankingsets[cache_key] = filtered_ranker
            
            return filtered_ranker

    def get_metadata(self) -> Dict[str, Any]:
        return self._metadata

    def get_entry(self, idx: int) -> Optional[Dict[str, Any]]:
        return self._metadata.get(str(idx))

    def get_smiles(self, idx: int) -> Optional[str]:
        entry = self.get_entry(idx)
        if not entry:
            return None
        smi = entry.get("canonical_3d_smiles")
        if smi and smi != "N/A":
            return smi
        return entry.get("canonical_2d_smiles")

    def get_exact_mass(self, idx: int) -> Optional[float]:
        return _mw_from_entry(self.get_entry(idx))

    def within_mw_range(
        self,
        idx: int,
        mw_min: Optional[float],
        mw_max: Optional[float],
    ) -> bool:
        if mw_min is None and mw_max is None:
            return True
        mass = self.get_exact_mass(idx)
        if mass is None:
            return True
        if mw_min is not None and mass < mw_min:
            return False
        if mw_max is not None and mass > mw_max:
            return False
        return True

    def indices_in_mw_range(
        self,
        mw_min: Optional[float],
        mw_max: Optional[float],
    ) -> List[int]:
        """
        Return list of indices whose MW is in [mw_min, mw_max]. Uses
        precomputed MW index (O(n) once per model load). No filter -> all indices.
        """
        self._ensure_mw_index()
        assert self._num_rows is not None
        if mw_min is None and mw_max is None:
            out = list(range(self._num_rows))
            return out
        assert self._mw_sorted is not None and self._mw_by_idx is not None
        if not self._mw_sorted:
            out = list(range(self._num_rows))
            return out
        mw_vals = [t[0] for t in self._mw_sorted]
        lo = bisect.bisect_left(mw_vals, mw_min) if mw_min is not None else 0
        hi = (
            bisect.bisect_right(mw_vals, mw_max)
            if mw_max is not None
            else len(mw_vals)
        )
        kept = [self._mw_sorted[i][1] for i in range(lo, hi)]
        # Indices with no mw are treated as passing the filter
        no_mw = [i for i in range(self._num_rows) if i not in self._mw_by_idx]
        kept = sorted(set(kept) | set(no_mw))
        return kept
