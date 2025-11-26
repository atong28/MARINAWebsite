import json
import threading
from typing import Any, Dict, Optional

from src.config import METADATA_PATH

class MetadataService:
    """Service for accessing molecule metadata. No caching - read on demand."""
    _instance: Optional['MetadataService'] = None
    _lock = threading.Lock()

    def __init__(self, path: str = METADATA_PATH) -> None:
        self._path = path
        self._data: Optional[Dict[str, Any]] = None
        self._data_lock = threading.Lock()

    @classmethod
    def instance(cls) -> "MetadataService":
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = MetadataService()
        return cls._instance

    def _ensure_loaded(self) -> None:
        """Load metadata file if not already loaded."""
        if self._data is None:
            with self._data_lock:
                if self._data is None:
                    with open(self._path, "r") as f:
                        self._data = json.load(f)

    def get_entry(self, idx: int) -> Optional[Dict[str, Any]]:
        """Get metadata entry for a given index."""
        self._ensure_loaded()
        assert self._data is not None
        return self._data.get(str(idx))

    def get_smiles(self, idx: int) -> Optional[str]:
        """Get SMILES string for a given index."""
        entry = self.get_entry(idx)
        if not entry:
            return None
        smi = entry.get("canonical_3d_smiles")
        if smi and smi != "N/A":
            return smi
        return entry.get("canonical_2d_smiles")

    def get_exact_mass(self, idx: int) -> Optional[float]:
        """
        Get exact mass / MW for a given index.

        Assumes each metadata entry has a root-level key \"mw\" storing
        the molecular weight. Returns None if it cannot be determined.
        """
        entry = self.get_entry(idx)
        if not entry:
            return None
        
        mass = entry.get("mw")
        if isinstance(mass, (int, float)):
            try:
                return float(mass)
            except Exception:
                return None
        return None

    def within_mw_range(
        self,
        idx: int,
        mw_min: Optional[float],
        mw_max: Optional[float],
    ) -> bool:
        """
        Check whether the molecule at idx lies within the given molecular weight range.

        If no bounds are provided, always returns True.
        If mass cannot be determined, the entry is treated as passing the filter
        to avoid silently discarding valid structures.
        """
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

