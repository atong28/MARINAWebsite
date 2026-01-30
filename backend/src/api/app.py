"""
FastAPI application entry point for MARINA backend.
"""
import asyncio
import logging
import multiprocessing as mp
import os
import signal
import time
from fastapi import FastAPI, Request

from src.api.middleware.cors import setup_cors
from src.api.middleware.error_handler import setup_error_handlers
from src.api.middleware.rate_limit import setup_rate_limiting
from src.api.middleware.request_id import setup_request_id_middleware
from src.api.middleware.request_logging import setup_request_logging

from src.services.model_service import ModelService
from src.services.compute_pool import get_compute_pool, shutdown_compute_pool
from src.config import (
    MAX_COMPUTE_QUEUE,
    MAX_COMPUTE_WORKERS,
    MAX_CONCURRENT_HEAVY_OPS,
    PORT,
    HOST,
)

import uvicorn

# Initialize logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - [%(request_id)s] - %(levelname)s - %(message)s'
)
_old_factory = logging.getLogRecordFactory()


def _record_factory(*args, **kwargs):
    record = _old_factory(*args, **kwargs)
    if not hasattr(record, "request_id"):
        record.request_id = "-"
    return record


logging.setLogRecordFactory(_record_factory)
logger = logging.getLogger(__name__)

# Create FastAPI app
# Note: OpenAPI docs URLs are set after routers are included to ensure /api prefix
app = FastAPI(
    title="MARINA API",
    description="Molecular Structure Annotator API",
    version="2.0.0"
)

# Setup middleware
setup_cors(app)
setup_error_handlers(app)
limiter = setup_rate_limiting(app)
setup_request_id_middleware(app)
setup_request_logging(app)

# Optional semaphore to limit concurrent heavy ops (predict/analyze/ablation)
app.state.heavy_semaphore = (
    asyncio.Semaphore(MAX_CONCURRENT_HEAVY_OPS) if MAX_CONCURRENT_HEAVY_OPS > 0 else None
)
app.state.compute_pool = get_compute_pool(MAX_COMPUTE_WORKERS, MAX_COMPUTE_QUEUE)


def _handle_signal(signum, frame) -> None:
    logger.error("Received signal %s; killing compute workers immediately", signum)
    for child in mp.active_children():
        try:
            child.terminate()
            child.join(timeout=0.5)
        except Exception:
            pass
        if child.is_alive() and child.pid:
            os.kill(child.pid, signal.SIGKILL)
    os._exit(0)


signal.signal(signal.SIGINT, _handle_signal)
signal.signal(signal.SIGTERM, _handle_signal)


async def run_heavy(request: Request, coro):
    """Run coroutine under heavy_semaphore if set (for predict/analyze/ablation)."""
    sem = getattr(request.app.state, "heavy_semaphore", None)
    if sem is not None:
        async with sem:
            return await coro
    return await coro


from src.api.routes import (
    health,
    logs,
    models as models_router,
    predict,
    smiles_search,
    analyze,
    secondary_retrieval,
    ablation,
    fingerprints,
)

# Include routers with /api prefix
app.include_router(health.router, prefix="/api", tags=["health"])
app.include_router(logs.router, prefix="/api", tags=["logs"])
app.include_router(models_router.router, prefix="/api", tags=["models"])
app.include_router(predict.router, prefix="/api", tags=["prediction"])
app.include_router(smiles_search.router, prefix="/api", tags=["search"])
app.include_router(analyze.router, prefix="/api", tags=["analysis"])
app.include_router(secondary_retrieval.router, prefix="/api", tags=["retrieval"])
app.include_router(ablation.router, prefix="/api", tags=["analysis"])
app.include_router(fingerprints.router, prefix="/api", tags=["fingerprints"])

# Track server start time
_server_start_time = time.time()

# Initialize model in background on startup
@app.on_event("startup")
async def startup_event():
    """Bootstrap from models.json; preload models according to PRELOAD_MODELS."""
    logger.info("Starting MARINA backend...")
    from src.config import PRELOAD_MODELS
    from src.services.model_manifest import (
        get_default_model_id,
        load_models_json,
        list_models,
    )

    load_models_json()
    entries = list_models()
    if not entries:
        logger.warning("No models configured in models.json; skipping model preload.")
        return

    # Resolve which entries to preload
    preload_val = PRELOAD_MODELS.lower()
    if not PRELOAD_MODELS or preload_val == "default":
        default_id = get_default_model_id()
        to_preload = [e for e in entries if e.id == default_id]
        if not to_preload:
            to_preload = [entries[0]]
            logger.info("Default model %s not in list; preloading first: %s", default_id, entries[0].id)
    elif preload_val == "all":
        to_preload = entries
    else:
        # Comma-separated list of model ids
        want = {s.strip() for s in PRELOAD_MODELS.split(",") if s.strip()}
        to_preload = [e for e in entries if e.id in want]

    from src.services.model_registry import pin_model

    model_service = ModelService.instance()
    for e in to_preload:
        try:
            logger.info("Preloading model %s (type=%s)...", e.id, e.type)
            model_service.preload_resources(e.id)
            pin_model(e.id)  # Pinned so LRU never evicts preloaded models
            logger.info("Model %s (type=%s) loaded successfully", e.id, e.type)
        except Exception as exc:
            logger.error(
                "Failed to preload model %s (type=%s): %s", e.id, e.type, exc, exc_info=True
            )
            raise


@app.on_event("shutdown")
async def shutdown_event():
    """Ensure compute workers are terminated on shutdown."""
    await shutdown_compute_pool()

if __name__ == "__main__":
    uvicorn.run(
        "src.api.app:app",
        host=HOST,
        port=PORT,
        reload=True,  # Enable hot reload for development
        log_level="info"
    )

