# MARINA TypeScript Rebuild

This is the rebuilt version of MARINA using TypeScript, FastAPI, and React with modern design principles.

## Structure

```
├── backend/          # FastAPI backend (Python)
│   ├── src/         # Source code
│   ├── tests/       # Test files
│   ├── pixi.toml    # Pixi environment configuration
│   └── requirements.txt  # Python dependencies (for Docker)
├── frontend/         # React + TypeScript frontend
│   ├── src/         # Source code
│   └── package.json # Node dependencies
├── shared/           # Shared types and schemas
└── docker-compose.yml
```

## Development Setup

### Prerequisites

- Docker and Docker Compose (for containerized development)
- Node.js 20+ (for local frontend development)
- Python 3.11+ (for local backend development)
- Pixi (recommended for backend dependency management): `curl -fsSL https://pixi.sh/install.sh | bash`

### Initial Setup

1. **Create environment file (Single Source of Truth for Ports):**
   ```bash
   cp .env.example .env
   # Edit .env to customize ports if needed (defaults: BACKEND_PORT=5000, FRONTEND_PORT=3000)
   ```
   
   This `.env` file at the project root is the **single source of truth** for all port configuration:
   - `BACKEND_PORT`: Used by Docker Compose, Pixi tasks, and backend config
   - `FRONTEND_PORT`: Used by Vite dev server and Docker Compose
   - `VITE_API_BASE`: Used by frontend to connect to backend
   
   All components automatically read from this file - no need to configure ports in multiple places!

### Backend Development

1. **With Pixi (recommended for local development):**
   ```bash
   # Set up environment (create .env from .env.example if needed)
   cd backend
   pixi install          # Install dependencies
   pixi run dev          # Start with hot reload (uses BACKEND_PORT from .env)
   ```
   The backend will be available at `http://localhost:${BACKEND_PORT:-5000}`
   - Uses Pixi for reproducible environment
   - All dependencies (PyTorch, RDKit, etc.) are managed by Pixi
   - Hot reload enabled automatically
   - Port configured via `BACKEND_PORT` in root `.env` file

2. **With Docker (for containerized development):**
   ```bash
   # Ensure .env file exists at project root
   docker-compose up backend
   ```
   The backend will be available at `http://localhost:${BACKEND_PORT:-5000}`
   - Source code is mounted for hot reload
   - PyTorch and RDKit are pre-installed in the image
   - Port configured via `BACKEND_PORT` in root `.env` file

3. **With pip (fallback):**
   ```bash
   cd backend
   pip install -r requirements.txt
   # Port is read from BACKEND_PORT env var or defaults to 5000
   uvicorn src.api.app:app --reload --host ${BACKEND_HOST:-0.0.0.0} --port ${BACKEND_PORT:-5000}
   ```

### Frontend Development

1. **With Docker:**
   ```bash
   # Ensure .env file exists at project root
   docker-compose up frontend
   ```
   The frontend will be available at `http://localhost:${FRONTEND_PORT:-3000}`
   - Port configured via `FRONTEND_PORT` in root `.env` file

2. **Local development (recommended for hot reload):**
   ```bash
   cd frontend
   npm install
   # Port is read from FRONTEND_PORT env var or defaults to 3000
   npm run dev
   ```
   The frontend will be available at `http://localhost:${FRONTEND_PORT:-3000}` with Vite hot reload
   - Port and API proxy target configured via `.env` file at project root

### Full Stack Development

```bash
docker-compose up
```

This starts both backend and frontend services.

## Key Features

### Backend
- **FastAPI** with automatic OpenAPI documentation
- **Pydantic** models for request/response validation
- **No caching** - all computation is on-demand
- **Type hints** throughout
- **Hot reload** via volume mounts in Docker

### Frontend
- **React 18** with TypeScript
- **Vite** for fast builds and HMR
- **React Query** for data fetching
- **Zustand** for state management
- **x-data-spreadsheet** for spreadsheet functionality
- **Strict TypeScript** mode

## Migration Notes

- All caching has been removed from the backend
- Custom spreadsheet implementation replaced with x-data-spreadsheet
- Manual state management replaced with Zustand
- API client uses React Query instead of manual fetch wrappers
- Error handling uses FastAPI exception handlers and React Error Boundaries

## Production Build

```bash
# Build backend
cd backend
docker build -t marina-backend:latest .

# Build frontend
cd frontend
npm run build
# Frontend static files can be served by nginx or the FastAPI backend
```

## Testing

```bash
# Backend tests
cd backend
pixi run test          # With Pixi
# OR
pytest                 # With pip

# Frontend tests (when implemented)
cd frontend
npm test
```

## Environment Variables

**Single Source of Truth**: All port configuration is managed through a root `.env` file.

Create a `.env` file in the project root (see `.env.example`):

### Port Configuration (Single Source of Truth)

- `BACKEND_PORT`: Backend server port (default: `5000`)
- `BACKEND_HOST`: Backend server host (default: `0.0.0.0`)
- `FRONTEND_PORT`: Frontend dev server port (default: `3000`)
- `VITE_API_BASE`: Backend API URL for frontend (default: `http://localhost:5000`)

**Note**: These ports are used consistently across:
- Docker Compose (`docker-compose.yml`)
- Vite dev server (`frontend/vite.config.ts`)
- Frontend API client (`frontend/src/services/api.ts`)
- Pixi tasks (`backend/pixi.toml`)
- Backend config (`backend/src/config.py`)

### Other Configuration

- **Required:**
  - `DATA_DIR`: Path to data directory (default: `data`)
  - `METADATA_PATH`: Path to metadata.json (default: `data/metadata.json`)

- **Optional:**
  - `MOLECULE_IMG_SIZE`: Molecule image size in pixels (default: `400`)
  - `DEFAULT_TOP_K`: Default number of results (default: `10`)
  - `MAX_TOP_K`: Maximum number of results (default: `50`)
  - `HIGHLIGHT_DEBUG`: Enable debug logging for highlighting (default: `false`)
  - `RDKIT_ENABLED`: Enable RDKit rendering (default: `true`)

## Troubleshooting

### Import Errors

If you see import errors when running the backend:
- Ensure you're running from the correct directory
- Check that `sys.path` adjustments in route files are working
- Verify all `__init__.py` files exist in package directories

### Docker Issues

**Port conflicts:**
- Change `BACKEND_PORT` or `FRONTEND_PORT` in root `.env` file
- Or stop other services using those ports
- All port references will automatically use the new values

**Volume mount issues:**
- Ensure paths in `docker-compose.yml` are correct relative to the compose file
- Check file permissions on mounted directories

**Model loading fails:**
- Verify `data/` directory is mounted correctly
- Check that model files exist in the data directory
- Review logs: `docker-compose logs backend`

### Model Loading

If the model fails to load:
- Check that `data/best.ckpt` exists
- Verify `data/params.json` is present
- Ensure `data/metadata.json` is accessible
- Check logs for specific error messages

### Hot Reload Not Working

**Backend:**
- Verify volume mount: `./backend/src:/app/src:rw` (should be `rw` not `ro`)
- Check that uvicorn is running with `--reload` flag
- Restart container: `docker-compose restart backend`

**Frontend:**
- For development, run `npm run dev` locally instead of using Docker
- Or use a dev Dockerfile that runs Vite dev server

## Architecture Overview

### Backend Structure

- `src/api/` - FastAPI routes and middleware
- `src/services/` - Business logic (model, rendering, metadata)
- `src/domain/` - Domain models and core logic
  - `models/` - Pydantic request/response models
  - `fingerprint/` - Fingerprint loading and utilities
  - `drawing/` - Molecule drawing and visualization
- `src/config.py` - Configuration management

### Frontend Structure

- `src/pages/` - Top-level page components
- `src/components/` - Reusable UI components
  - `common/` - Shared components (Button, Card, ErrorBoundary)
  - `spectral/` - Spectral input components
  - `spreadsheet/` - Spreadsheet table component
  - `results/` - Results display components
  - `analysis/` - Analysis page components
- `src/services/` - API client and React Query hooks
- `src/store/` - Zustand state management

### Data Flow

1. User inputs spectral data or SMILES
2. Frontend sends request to FastAPI backend
3. Backend processes with ML model
4. Results returned with molecule images and metadata
5. Frontend displays results and allows analysis
6. Analysis page shows fingerprint visualization and overlays

## Contributing

### Code Style

- **Backend**: Follow PEP 8, use type hints
- **Frontend**: Follow ESLint rules, use TypeScript strict mode

### Commit Messages

Use conventional commits:
- `feat:` for new features
- `fix:` for bug fixes
- `docs:` for documentation
- `refactor:` for code refactoring

### Pull Request Process

1. Create feature branch from `main`
2. Make changes with tests
3. Ensure all tests pass
4. Update documentation if needed
5. Submit PR with clear description

## Deployment

### Production Build

**Backend:**
```bash
cd backend
docker build -t marina-backend:latest .
```

**Frontend:**
```bash
cd frontend
npm run build
# Static files in dist/ can be served by nginx or FastAPI
```

### Environment Setup

- Set all required environment variables
- Ensure data directory is accessible
- Configure reverse proxy (nginx) if needed
- Set up SSL certificates for HTTPS

### Scaling Considerations

- Backend can be scaled horizontally (stateless)
- Use load balancer for multiple backend instances
- Consider Redis for rate limiting in production
- Frontend is static and can be served via CDN

