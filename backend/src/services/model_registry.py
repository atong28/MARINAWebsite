"""
Model registry: map model_id -> ModelSession for multi-model support.
Supports optional MAX_LOADED_MODELS with LRU eviction; default model is pinned.
"""
from __future__ import annotations

import logging
import threading
from typing import Dict, List, Optional, Set, TYPE_CHECKING

from src.config import DEFAULT_MODEL_ID, MAX_LOADED_MODELS, MODEL_ROOT
from src.services.model_manifest import get_model_info, load_models_json

if TYPE_CHECKING:
    from src.domain.model_session import ModelSession

logger = logging.getLogger(__name__)

_registry: Dict[str, "ModelSession"] = {}
_registry_lock = threading.Lock()
# LRU order: oldest first (evict from front). Updated on get() and register().
_lru: List[str] = []
# Pinned model ids (e.g. default) are never evicted.
_pinned: Set[str] = set()


def _touch_lru(model_id: str) -> None:
    """Move model_id to end of LRU (most recently used). Call with _registry_lock held."""
    if model_id in _lru:
        _lru.remove(model_id)
    _lru.append(model_id)


def _evict_lru_if_needed(requested_id: str) -> None:
    """
    If MAX_LOADED_MODELS is set and we're at capacity, evict the least recently
    used model that is not pinned. Call with _registry_lock held (only when
    we are about to load a new model and requested_id is not in _registry).
    """
    if MAX_LOADED_MODELS <= 0:
        return
    while len(_registry) >= MAX_LOADED_MODELS and _lru:
        candidate = _lru.pop(0)
        if candidate in _pinned or candidate not in _registry:
            if candidate in _registry:
                _lru.append(candidate)
            continue
        session = _registry.pop(candidate)
        logger.info("Evicted model %s from registry (LRU, max_loaded=%s)", candidate, MAX_LOADED_MODELS)
        del session
        break


def get(model_id: str) -> Optional["ModelSession"]:
    """Return the ModelSession for model_id, or None if not loaded."""
    with _registry_lock:
        session = _registry.get(model_id)
        if session is not None:
            _touch_lru(model_id)
        return session


def ensure_loaded(model_id: str, model_root: Optional[str] = None) -> "ModelSession":
    """
    Load model_id if not already in registry. Uses model_root if provided.
    Otherwise uses model manifest (models.json). If at MAX_LOADED_MODELS,
    evicts the least recently used non-pinned model first.
    """
    with _registry_lock:
        if model_id in _registry:
            _touch_lru(model_id)
            return _registry[model_id]
        _evict_lru_if_needed(model_id)

    root: Optional[str] = model_root
    model_type: Optional[str] = None
    if root is None:
        entries = load_models_json()
        if entries:
            info = get_model_info(model_id)
            if info is not None:
                # Allow all manifest-supported types (e.g., 'marina', 'spectre').
                model_type = info.type
                root = info.root
        if root is None and model_id == DEFAULT_MODEL_ID:
            root = MODEL_ROOT
        if root is None:
            raise RuntimeError(
                f"Unknown model_id {model_id!r} and no model_root provided"
            )

    return register(model_id, root, model_type=model_type)


def pin_model(model_id: str) -> None:
    """Mark model_id as pinned so it is never evicted by LRU. Call after preload."""
    with _registry_lock:
        _pinned.add(model_id)


def register(model_id: str, model_root: str, model_type: Optional[str] = None) -> "ModelSession":
    """
    Load model from model_root, register under model_id, and return the session.
    Caller must not hold _registry_lock.
    """
    from src.domain.model_session import ModelSession

    session = ModelSession.from_model_root(model_root, model_type=model_type)
    with _registry_lock:
        _registry[model_id] = session
        _touch_lru(model_id)
        if model_id == DEFAULT_MODEL_ID:
            _pinned.add(model_id)
    logger.info("Registered model %s from %s", model_id, model_root)
    return session
