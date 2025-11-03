import threading
from typing import Optional

import torch

from src import predictor


class ModelService:
    _instance = None
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
        try:
            return predictor._model is not None
        except Exception:
            return False

    def ensure_loaded(self) -> None:
        if not self.is_ready():
            with self._init_lock:
                if not self.is_ready():
                    predictor.load_model()

    def get_model(self) -> torch.nn.Module:
        self.ensure_loaded()
        assert predictor._model is not None
        return predictor._model

    def get_fp_loader(self):
        self.ensure_loaded()
        return predictor._fp_loader

    def get_rankingset(self):
        self.ensure_loaded()
        if predictor._rankingset_cache is None:
            predictor._rankingset_cache = predictor._fp_loader.load_rankingset(predictor._args.fp_type)
        return predictor._rankingset_cache


