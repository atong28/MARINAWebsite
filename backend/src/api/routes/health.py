"""
Health check endpoint.
"""
import asyncio
import time

from fastapi import APIRouter, status
from src.services.model_service import ModelService

router = APIRouter()

# Track server start time (set at module import)
_server_start_time = time.time()


@router.get("/health/live", status_code=status.HTTP_200_OK)
async def health_live():
    """
    Minimal liveness probe: process is up and event loop is responding.
    No model registry or is_ready() — use for Docker/orchestrator healthchecks
    so the container is not marked unhealthy when busy or during model load.
    """
    return {"status": "alive"}


@router.get("/health", status_code=status.HTTP_200_OK)
async def health():
    """Health check endpoint with model status."""
    uptime = time.time() - _server_start_time
    base_response = {
        'uptime_seconds': round(uptime, 1),
        'server_start_time': _server_start_time
    }
    
    model_service = ModelService.instance()
    is_ready = await asyncio.to_thread(model_service.is_ready)
    
    if is_ready:
        return {
            **base_response,
            'status': 'ok',
            'model_loaded': True,
            'message': 'Model is ready for predictions'
        }
    else:
        return {
            **base_response,
            'status': 'initializing',
            'model_loaded': False,
            'message': 'Model initialization in progress'
        }

