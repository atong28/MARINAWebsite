"""
Health check endpoint.
"""
import time

from fastapi import APIRouter, status
from fastapi.responses import JSONResponse
from src.services.model_service import ModelService

router = APIRouter()

# Track server start time (set at module import)
_server_start_time = time.time()


@router.get("/health", status_code=status.HTTP_200_OK)
async def health():
    """Health check endpoint with model status."""
    uptime = time.time() - _server_start_time
    base_response = {
        'uptime_seconds': round(uptime, 1),
        'server_start_time': _server_start_time
    }
    
    model_service = ModelService.instance()
    is_ready = model_service.is_ready()
    
    if is_ready:
        return {
            **base_response,
            'status': 'ok',
            'model_loaded': True,
            'message': 'Model is ready for predictions'
        }
    else:
        return JSONResponse(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            content={
                **base_response,
                'status': 'initializing',
                'model_loaded': False,
                'message': 'Model initialization in progress'
            }
        )

