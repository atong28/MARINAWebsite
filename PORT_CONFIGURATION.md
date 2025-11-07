# Port Configuration - Single Source of Truth

This document describes how port configuration works across the MARINA project.

## Overview

All port configuration is managed through a **root `.env` file** at the project root. This is the single source of truth for all ports.

## Environment Variables

Create a `.env` file at the project root with the following variables:

```bash
# Backend
BACKEND_PORT=5000
BACKEND_HOST=0.0.0.0

# Frontend
FRONTEND_PORT=3000

# API Connection
VITE_API_BASE=http://localhost:5000
```

## Where Ports Are Used

### 1. Docker Compose (`docker-compose.yml`)

- **Backend port mapping**: `${BACKEND_PORT:-5000}:5000`
- **Frontend port mapping**: `${FRONTEND_PORT:-3000}:80`
- **Backend environment**: `PORT=${BACKEND_PORT:-5000}`, `HOST=${BACKEND_HOST:-0.0.0.0}`
- **Frontend environment**: `VITE_API_BASE=${VITE_API_BASE:-...}`

### 2. Backend Config (`backend/src/config.py`)

- Reads `BACKEND_PORT` (falls back to `PORT` for compatibility)
- Reads `BACKEND_HOST` (falls back to `HOST` for compatibility)
- Used by `backend/src/api/app.py` when running with uvicorn

### 3. Pixi Tasks (`backend/pixi.toml`)

- `start` and `dev` tasks use `${BACKEND_PORT:-5000}` and `${BACKEND_HOST:-0.0.0.0}`
- Environment variables should be exported before running or set in `.env` file

### 4. Vite Config (`frontend/vite.config.ts`)

- Reads `FRONTEND_PORT` for dev server port (default: 3000)
- Reads `BACKEND_PORT` to construct API proxy target
- Reads `VITE_API_BASE` for API base URL
- Loads `.env` from project root (parent of `frontend/`)

### 5. Frontend API Client (`frontend/src/services/api.ts`)

- Reads `VITE_API_BASE` from environment (set at build time)
- Falls back to `http://localhost:5000` if not set

## How to Change Ports

1. **Edit root `.env` file:**
   ```bash
   BACKEND_PORT=8080
   FRONTEND_PORT=3001
   VITE_API_BASE=http://localhost:8080
   ```

2. **All components will automatically use the new ports:**
   - Docker Compose will map the new ports
   - Vite dev server will start on the new frontend port
   - Frontend API client will connect to the new backend port
   - Pixi tasks will use the new backend port

## Default Values

If `.env` file doesn't exist or variables are not set:
- `BACKEND_PORT`: 5000
- `BACKEND_HOST`: 0.0.0.0
- `FRONTEND_PORT`: 3000
- `VITE_API_BASE`: http://localhost:5000

## Verification

To verify ports are configured correctly:

1. Check `.env` file exists at project root
2. Run `docker-compose config` to see resolved port mappings
3. Check Vite dev server starts on expected port
4. Verify frontend can connect to backend API

