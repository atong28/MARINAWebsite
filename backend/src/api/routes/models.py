"""
Model discovery: GET /api/models returns available models from models.json.
"""
from fastapi import APIRouter, status

from src.services.model_manifest import get_default_model_id, list_models
from src.services.model_registry import get

router = APIRouter()


@router.get("/models", status_code=status.HTTP_200_OK)
async def models():
    """
    List available models from models.json with loaded status.
    """
    entries = list_models()
    default_id = get_default_model_id()
    models_out = []
    for e in entries:
        loaded = get(e.id) is not None
        models_out.append(
            {
                "id": e.id,
                "root": e.root_rel,
                "type": e.type,
                "default": e.default,
                "loaded": loaded,
                "display_name": e.display_name or e.id,
            }
        )
    return {"models": models_out, "default_model_id": default_id}
