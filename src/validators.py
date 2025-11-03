from typing import Tuple

from src.config import DEFAULT_TOP_K, MAX_TOP_K


def validate_k(raw_k) -> Tuple[int, str]:
    """
    Validate and clamp k, returning (k_value, error_message).
    error_message is empty string if no error.
    """
    try:
        k = int(raw_k)
    except Exception:
        k = DEFAULT_TOP_K
    if k <= 0:
        return k, 'k must be positive integer'
    if k > MAX_TOP_K:
        k = MAX_TOP_K
    return k, ''


def validate_smiles(smiles: str) -> Tuple[str, str]:
    """
    Basic validation for SMILES input (non-empty, trimmed). Returns (value, error_message).
    """
    if smiles is None:
        return '', 'SMILES string required'
    s = str(smiles).strip()
    if not s:
        return '', 'SMILES string cannot be empty'
    return s, ''


def require_fields(payload: dict, *fields) -> Tuple[bool, str]:
    missing = [f for f in fields if (f not in payload or payload.get(f) in (None, ''))]
    if missing:
        return False, f"Missing required fields: {', '.join(missing)}"
    return True, ''


def require_json_fields(*fields):
    """Decorator to enforce JSON body and required fields; returns JSON error on failure."""
    def decorator(fn):
        def wrapper(*args, **kwargs):
            try:
                from flask import request
                payload = request.get_json(silent=True) or {}
                ok, msg = require_fields(payload, *fields) if fields else (True, '')
                if not ok:
                    from src.http import json_error_response
                    return json_error_response(msg, status=400)
                return fn(*args, **kwargs)
            except Exception as e:
                from src.http import json_error_response
                return json_error_response(str(e) or 'Invalid request', status=400)
        wrapper.__name__ = getattr(fn, '__name__', 'wrapped')
        return wrapper
    return decorator


