# Multi-stage Dockerfile for MARINA backend using Pixi
# Stage 1: Install Pixi and dependencies
FROM ubuntu:24.04 as pixi-base

# Install system dependencies required for pixi and conda/mamba
RUN apt-get update && apt-get install -y \
    curl \
    ca-certificates \
    git \
    bzip2 \
    && rm -rf /var/lib/apt/lists/*

# Install Pixi
RUN curl -fsSL https://pixi.sh/install.sh | bash
ENV PATH="/root/.pixi/bin:$PATH"

WORKDIR /app

# Copy pixi configuration files
COPY pixi.toml ./

# Install dependencies using pixi
# This creates the conda environment with all dependencies
RUN pixi install

# Stage 2: Application code
FROM pixi-base

WORKDIR /app

# Copy application code (will be overridden by volume mount in dev)
COPY src/ ./src/

# Ensure pixi is in PATH
ENV PATH="/root/.pixi/bin:$PATH"

# CRITICAL: Set LD_LIBRARY_PATH at Docker level to prioritize pixi's libraries
# This ensures conda libraries are found BEFORE system libraries
ENV LD_LIBRARY_PATH="/app/.pixi/envs/default/lib:${LD_LIBRARY_PATH:-}"

# Also set KMP_DUPLICATE_LIB_OK to handle any remaining conflicts
ENV KMP_DUPLICATE_LIB_OK="TRUE"

# Expose port
EXPOSE 5000

# Default command using pixi run task
CMD ["pixi", "run", "start"]