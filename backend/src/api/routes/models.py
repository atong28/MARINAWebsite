"""
Model discovery: GET /api/models returns available models from models.json.
Model loading: POST /api/models/{model_id}/load explicitly loads a model.
"""
import logging
from fastapi import APIRouter, HTTPException, status

from src.services.model_manifest import get_default_model_id, get_model_info, list_models
from src.services.model_registry import get
from src.services.model_service import ModelService

logger = logging.getLogger(__name__)
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
            }
        )
    return {"models": models_out, "default_model_id": default_id}


@router.post("/models/{model_id}/load", status_code=status.HTTP_200_OK)
async def load_model(model_id: str):
    """
    Explicitly load a model and preload its resources (rankingset, MW index).
    Returns status indicating if model was loaded, already loaded, or if an error occurred.
    """
    # Validate model exists in manifest
    info = get_model_info(model_id)
    if info is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Model '{model_id}' not found in manifest"
        )

    # Check if already loaded
    model_service = ModelService.instance()
    if model_service.is_ready(model_id):
        return {
            "model_id": model_id,
            "status": "already_loaded",
            "message": f"Model '{model_id}' is already loaded"
        }

    try:
        # Load model and preload resources
        logger.info("Loading model %s and preloading resources...", model_id)
        model_service.preload_resources(model_id)
        logger.info("Successfully loaded model %s", model_id)
        return {
            "model_id": model_id,
            "status": "loaded",
            "message": f"Model '{model_id}' loaded successfully"
        }
    except Exception as e:
        logger.error("Failed to load model %s: %s", model_id, e, exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to load model '{model_id}': {str(e)}"
        )
