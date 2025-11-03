import json
import threading
from typing import Any, Dict, Optional
from src.config import METADATA_PATH


class MetadataService:
    _instance = None
    _lock = threading.Lock()

    def __init__(self, path: str = METADATA_PATH) -> None:
        self._path = path
        self._cache: Optional[Dict[str, Any]] = None
        self._cache_lock = threading.Lock()

    @classmethod
    def instance(cls) -> "MetadataService":
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = MetadataService()
        return cls._instance

    def _ensure_loaded(self) -> None:
        if self._cache is None:
            with self._cache_lock:
                if self._cache is None:
                    with open(self._path, "r") as f:
                        self._cache = json.load(f)

    def get_entry(self, idx: int) -> Optional[Dict[str, Any]]:
        self._ensure_loaded()
        assert self._cache is not None
        return self._cache.get(str(idx))

    def get_smiles(self, idx: int) -> Optional[str]:
        entry = self.get_entry(idx)
        if not entry:
            return None
        smi = entry.get("canonical_3d_smiles")
        if smi and smi != "N/A":
            return smi
        return entry.get("canonical_2d_smiles")

    def preload(self) -> None:
        self._ensure_loaded()

    def reload(self) -> None:
        with self._cache_lock:
            with open(self._path, "r") as f:
                self._cache = json.load(f)


