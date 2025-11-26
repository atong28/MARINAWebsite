import threading
from typing import Optional
import torch

from src.domain.predictor import get_session

class ModelService:
    """Service for managing ML model and related resources. No caching - compute on demand."""
    _instance: Optional['ModelService'] = None
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

    def is_ready(self) -> bool:
        """Check if model is loaded and ready."""
        try:
            session = get_session()
            return session.model is not None
        except Exception:
            return False

    def ensure_loaded(self) -> None:
        """Ensure model is loaded. Loads if not already loaded."""
        if not self.is_ready():
            with self._init_lock:
                if not self.is_ready():
                    get_session()

    def get_model(self) -> torch.nn.Module:
        """Get the loaded model. Loads if necessary."""
        self.ensure_loaded()
        session = get_session()
        return session.model

    def get_fp_loader(self):
        """Get the fingerprint loader. Loads model if necessary."""
        self.ensure_loaded()
        session = get_session()
        return session.fp_loader

    def get_rankingset(self):
        """Get the rankingset. Loads model and rankingset if necessary."""
        self.ensure_loaded()
        session = get_session()
        ranker = session.get_rankingset()
        return ranker.data

