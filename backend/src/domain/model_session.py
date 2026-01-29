import bisect
import json
import logging
import os
import threading
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import torch

from src.config import CKPT_PATH, MODEL_ROOT, PARAMS_PATH, RETRIEVAL_PATH
from src.domain.fingerprint.fp_loader import FPLoader, make_fp_loader_for_model_root
from src.domain.model import MARINA
from src.domain.ranker import RankingSet
from src.domain.settings import MARINAArgs


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
    Encapsulates the loaded MARINA model, fingerprint loader, rankingset store,
    metadata, and MW index for O(1) range filtering.
    """

    args: MARINAArgs
    model: MARINA
    fp_loader: FPLoader
    fp_type: str
    ckpt_path: str
    params_path: str
    retrieval_path: str
    model_root: str
    _metadata: Dict[str, Any]
    _metadata_path: str

    _rankingset_store: Optional[torch.Tensor] = None
    _mw_sorted: Optional[List[tuple]] = None  # (mw, idx) sorted by mw
    _mw_by_idx: Optional[Dict[int, float]] = None
    _num_rows: Optional[int] = None
    _lock: threading.Lock = field(default_factory=threading.Lock, repr=False)

    @classmethod
    def from_model_root(cls, model_root: str) -> "ModelSession":
        """
        Load model and metadata from a per-model directory (e.g. data/marina_best/).
        Uses make_fp_loader_for_model_root so dataset_root and retrieval both
        point under model_root.
        """
        ckpt_path = os.path.join(model_root, "best.ckpt")
        params_path = os.path.join(model_root, "params.json")
        retrieval_path = os.path.join(model_root, "retrieval.pkl")
        metadata_path = os.path.join(model_root, "metadata.json")

        logger.info("Loading model params from %s", params_path)
        with open(params_path, "r") as f:
            params = json.load(f)
        args = MARINAArgs(**params)

        fp_loader = make_fp_loader_for_model_root(
            model_root,
            fp_type=args.fp_type,
            entropy_out_dim=args.out_dim,
            max_radius=6,
        )

        logger.info("Instantiating MARINA model")
        model = MARINA(args, fp_loader)

        logger.info("Loading checkpoint from %s (cpu)", ckpt_path)
        ckpt = torch.load(ckpt_path, map_location="cpu")
        sd = ckpt.get("state_dict", ckpt)
        res = model.load_state_dict(sd, strict=False)
        logger.info("Model load_state_dict result: %s", res)

        model.to(torch.device("cpu"))
        model.eval()
        torch.set_grad_enabled(False)

        metadata = _load_metadata(metadata_path)

        return cls(
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
        self._ensure_mw_index()
        assert self._rankingset_store is not None
        return RankingSet(store=self._rankingset_store, metric="cosine")

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
            return list(range(self._num_rows))
        assert self._mw_sorted is not None and self._mw_by_idx is not None
        if not self._mw_sorted:
            return list(range(self._num_rows))
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
