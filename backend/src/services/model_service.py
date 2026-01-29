import threading
from typing import Optional

import torch

from src.services.model_manifest import get_default_model_id
from src.services.model_registry import ensure_loaded, get

# Avoid circular import: ModelSession used for type hints only
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from src.domain.model_session import ModelSession


class ModelService:
    """Service for managing ML model and related resources via model registry."""
    _instance: Optional["ModelService"] = None
    _lock = threading.Lock()

    def __init__(self) -> None:
        self._init_lock = threading.Lock()

    @classmethod
    def instance(cls) -> "ModelService":
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = ModelService()
        return cls._instance

    def _session(self, model_id: Optional[str] = None) -> "ModelSession":
        mid = model_id or get_default_model_id()
        s = get(mid)
        if s is None:
            raise RuntimeError(f"Model {mid} not loaded. Call ensure_loaded first.")
        return s

    def is_ready(self, model_id: Optional[str] = None) -> bool:
        mid = model_id or get_default_model_id()
        try:
            s = get(mid)
            return s is not None and s.model is not None
        except Exception:
            return False

    def ensure_loaded(self, model_id: Optional[str] = None) -> None:
        mid = model_id or get_default_model_id()
        if not self.is_ready(mid):
            with self._init_lock:
                if not self.is_ready(mid):
                    ensure_loaded(mid)

    def get_session(self, model_id: Optional[str] = None) -> "ModelSession":
        self.ensure_loaded(model_id)
        return self._session(model_id)

    def get_model(self, model_id: Optional[str] = None) -> torch.nn.Module:
        return self.get_session(model_id).model

    def get_fp_loader(self, model_id: Optional[str] = None):
        return self.get_session(model_id).fp_loader

    def get_rankingset(self, model_id: Optional[str] = None):
        self.ensure_loaded(model_id)
        session = self._session(model_id)
        ranker = session.get_rankingset()
        return ranker.data
