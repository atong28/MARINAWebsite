"""
Model registry: map model_id -> ModelSession for multi-model support.
"""
import logging
import threading
from typing import Dict, Optional

from src.config import DEFAULT_MODEL_ID, MODEL_ROOT
from src.services.model_manifest import get_model_info, load_models_json

logger = logging.getLogger(__name__)

_registry: Dict[str, "ModelSession"] = {}
_registry_lock = threading.Lock()


def get(model_id: str) -> Optional["ModelSession"]:
    """Return the ModelSession for model_id, or None if not loaded."""
    with _registry_lock:
        return _registry.get(model_id)


def ensure_loaded(model_id: str, model_root: Optional[str] = None) -> "ModelSession":
    """
    Load model_id if not already in registry. Uses model_root if provided.
    Otherwise uses model manifest (models.json): only registers type=marina.
    If manifest missing or model_id not in manifest, falls back to MODEL_ROOT
    for default (DEFAULT_MODEL_ID).
    """
    with _registry_lock:
        if model_id in _registry:
            return _registry[model_id]

    root: Optional[str] = model_root
    if root is None:
        entries = load_models_json()
        if entries:
            info = get_model_info(model_id)
            if info is not None:
                if info.type != "marina":
                    raise RuntimeError(
                        f"Model type {info.type!r} not yet supported (model_id={model_id})"
                    )
                root = info.root
        if root is None and model_id == DEFAULT_MODEL_ID:
            root = MODEL_ROOT
        if root is None:
            raise RuntimeError(
                f"Unknown model_id {model_id!r} and no model_root provided"
            )

    return register(model_id, root)


def register(model_id: str, model_root: str) -> "ModelSession":
    """
    Load model from model_root, register under model_id, and return the session.
    Caller must not hold _registry_lock.
    """
    from src.domain.model_session import ModelSession

    session = ModelSession.from_model_root(model_root)
    with _registry_lock:
        _registry[model_id] = session
    logger.info("Registered model %s from %s", model_id, model_root)
    return session
