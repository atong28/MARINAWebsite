import threading
from typing import Optional
import torch

import src.domain.predictor as predictor_module

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
            return predictor_module._model is not None
        except Exception:
            return False

    def ensure_loaded(self) -> None:
        """Ensure model is loaded. Loads if not already loaded."""
        if not self.is_ready():
            with self._init_lock:
                if not self.is_ready():
                    predictor_module.load_model()

    def get_model(self) -> torch.nn.Module:
        """Get the loaded model. Loads if necessary."""
        self.ensure_loaded()
        assert predictor_module._model is not None
        return predictor_module._model

    def get_fp_loader(self):
        """Get the fingerprint loader. Loads model if necessary."""
        self.ensure_loaded()
        return predictor_module._fp_loader

    def get_rankingset(self):
        """Get the rankingset. Loads model and rankingset if necessary."""
        self.ensure_loaded()
        # Rankingset is loaded on-demand, not cached
        if predictor_module._rankingset_cache is None:
            if predictor_module._fp_loader is not None and predictor_module._args is not None:
                rankingset = predictor_module._fp_loader.load_rankingset(predictor_module._args.fp_type)
                return rankingset
        return predictor_module._rankingset_cache

