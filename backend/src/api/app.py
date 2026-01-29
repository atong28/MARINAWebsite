"""
FastAPI application entry point for MARINA backend.
"""
import logging
import time
from fastapi import FastAPI

from src.api.middleware.cors import setup_cors
from src.api.middleware.error_handler import setup_error_handlers
from src.api.middleware.rate_limit import setup_rate_limiting

from src.services.model_service import ModelService
from src.config import PORT, HOST

import uvicorn

# Initialize logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
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

from src.api.routes import (
    health,
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
    """Bootstrap from models.json when present; load all configured models."""
    logger.info("Starting MARINA backend...")
    from src.services.model_manifest import load_models_json, list_models

    load_models_json()
    entries = list_models()
    if not entries:
        logger.warning("No models configured in models.json; skipping model preload.")
        return

    model_service = ModelService.instance()
    for e in entries:
        try:
            logger.info("Preloading model %s (type=%s)...", e.id, e.type)
            model_service.preload_resources(e.id)
            logger.info("Model %s (type=%s) loaded successfully", e.id, e.type)
        except Exception as exc:
            logger.error(
                "Failed to preload model %s (type=%s): %s", e.id, e.type, exc, exc_info=True
            )
            # Fail-fast: if a configured model can't load, don't start the app.
            raise

if __name__ == "__main__":
    uvicorn.run(
        "src.api.app:app",
        host=HOST,
        port=PORT,
        reload=True,  # Enable hot reload for development
        log_level="info"
    )

