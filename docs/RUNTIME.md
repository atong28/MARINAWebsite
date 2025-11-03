# Runtime and Deployment

## Environment Variables
- PORT: HTTP port (default 5000)
- WORKERS: Gunicorn workers (default 4)
- THREADS: Gunicorn threads per worker (default 4)
- DEFAULT_TOP_K: Default top-k (default 10)
- MAX_TOP_K: Max top-k (default 50)
- MOLECULE_IMG_SIZE: Renderer image size (default 400)
- ENABLE_RENDERING_CACHE: Enable RDKit render cache (true/false)
- RENDERING_CACHE_SIZE: LRU cache entries for renderer (default 1000)

## Performance
- Compression: If `flask-compress` installed, responses are compressed automatically.
- Rate limiting: If `flask-limiter` installed, per-route limits are enabled.
- Metadata cache: Preloaded on startup; `MetadataService.reload()` to refresh.
- Renderer cache: LRU keyed by size and SMILES; call `clear_cache()` to drop.

## Concurrency
- Gunicorn workers × threads recommended starting point: 4×4.
- Adjust based on CPU and RDKit availability.

## API Docs
- Swagger UI: `/docs`
- Spec: `/openapi.json`


