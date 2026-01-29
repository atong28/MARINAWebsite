"""
Model manifest: load and validate data/models.json, resolve model roots, expose default.
"""
from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass
from typing import List, Optional

from src.config import DATA_DIR, DEFAULT_MODEL_ID, MODELS_JSON_PATH

logger = logging.getLogger(__name__)

SUPPORTED_TYPES = frozenset({"marina", "spectre"})


@dataclass
class ModelEntry:
    id: str
    root: str  # absolute path
    root_rel: str  # relative (as in JSON), for API
    type: str
    default: bool


_manifest: Optional[List[ModelEntry]] = None


def load_models_json(path: Optional[str] = None) -> List[ModelEntry]:
    """
    Load and validate models.json. Resolve root to absolute path via DATA_DIR.
    Returns empty list on missing file, invalid schema, or validation error.
    """
    global _manifest
    p = path if path is not None else MODELS_JSON_PATH
    if not os.path.isfile(p):
        logger.info("models.json not found at %s", p)
        return []

    try:
        with open(p, "r") as f:
            data = json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        logger.warning("Failed to load models.json: %s", e)
        return []

    models_raw = data.get("models")
    if not isinstance(models_raw, list) or len(models_raw) == 0:
        logger.warning("models.json: missing or empty 'models' array")
        return []

    entries: List[ModelEntry] = []
    seen_ids: set[str] = set()
    default_count = 0

    for i, m in enumerate(models_raw):
        if not isinstance(m, dict):
            logger.warning("models.json: models[%d] is not an object", i)
            continue
        mid = m.get("id")
        root = m.get("root")
        typ = m.get("type")
        default = m.get("default", False)

        if not isinstance(mid, str) or not mid:
            logger.warning("models.json: models[%d] missing or invalid 'id'", i)
            continue
        if mid in seen_ids:
            logger.warning("models.json: duplicate id %r", mid)
            continue
        seen_ids.add(mid)

        if not isinstance(root, str) or not root:
            logger.warning("models.json: model %r missing or invalid 'root'", mid)
            continue
        root_rel = root
        if os.path.isabs(root):
            root_abs = root
        else:
            root_abs = os.path.normpath(os.path.join(DATA_DIR, root))

        if not isinstance(typ, str) or typ not in SUPPORTED_TYPES:
            logger.warning(
                "models.json: model %r has invalid 'type' (expected marina|spectre)",
                mid,
            )
            continue
        if not isinstance(default, bool):
            default = False
        if default:
            default_count += 1

        entries.append(
            ModelEntry(
                id=mid,
                root=root_abs,
                root_rel=root_rel,
                type=typ,
                default=bool(default),
            )
        )

    if default_count != 1:
        logger.warning(
            "models.json: exactly one model must have default=true (got %d)",
            default_count,
        )
        return []

    _manifest = entries
    return entries


def _ensure_loaded() -> List[ModelEntry]:
    if _manifest is not None:
        return _manifest
    return load_models_json()


def get_default_model_id() -> str:
    """Default model id from manifest; else DEFAULT_MODEL_ID from config."""
    entries = _ensure_loaded()
    for e in entries:
        if e.default:
            return e.id
    return DEFAULT_MODEL_ID


def get_model_info(model_id: str) -> Optional[ModelEntry]:
    """Return manifest entry for model_id, or None."""
    entries = _ensure_loaded()
    for e in entries:
        if e.id == model_id:
            return e
    return None


def list_models() -> List[ModelEntry]:
    """Return all manifest entries (after load)."""
    return _ensure_loaded()


def resolve_and_validate_model_id(
    model_id: Optional[str],
) -> tuple[str, Optional[tuple[int, str]]]:
    """
    Resolve model_id from request (use default if omitted). Validate via manifest.
    Returns (model_id, None) on success, or (model_id, (status_code, detail)) on error.
    When manifest is missing, always use default and never error.
    """
    entries = _ensure_loaded()
    default_id = get_default_model_id()
    mid = model_id or default_id

    if not entries:
        return (default_id, None)

    info = get_model_info(mid)
    if info is None:
        return (mid, (400, f"Unknown model_id '{mid}'"))
    if info.type != "marina":
        return (mid, (501, f"Model type '{info.type}' is not yet supported."))
    return (mid, None)
