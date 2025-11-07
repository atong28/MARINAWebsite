# Multi-stage Dockerfile for MARINA backend
# Stage 1: Base dependencies (PyTorch, RDKit) - rarely changes
FROM python:3.11-slim as base-deps

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install heavy dependencies that rarely change
RUN pip install --no-cache-dir \
    torch>=2.0.0 \
    rdkit-pypi>=2023.3.1 \
    numpy>=1.24.0

# Stage 2: Application dependencies - changes more often
FROM base-deps as app-deps

WORKDIR /app

# Copy requirements and install
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Stage 3: Application code - changes frequently (mounted in dev)
FROM app-deps

WORKDIR /app

# Copy application code (will be overridden by volume mount in dev)
COPY src/ ./src/

# Expose port
EXPOSE 5000

# Default command (use uvicorn with reload for development)
CMD ["uvicorn", "src.api.app:app", "--host", "0.0.0.0", "--port", "5000", "--reload"]

