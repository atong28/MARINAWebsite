import json
import threading
from typing import Any, Dict, Optional

from rdkit import Chem
from rdkit.Chem import Descriptors

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
        Get or compute exact mass for a given index.

        Prefers any precomputed exact_mass field in metadata, but will fall back
        to computing from the best-available SMILES using RDKit. Returns None
        if mass cannot be determined.
        """
        entry = self.get_entry(idx)
        if not entry:
            return None

        mass = entry.get("exact_mass")
        if isinstance(mass, (int, float)):
            try:
                return float(mass)
            except Exception:
                pass

        smiles = self.get_smiles(idx)
        if not smiles:
            return None

        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None
            return float(Descriptors.ExactMolWt(mol))
        except Exception:
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

