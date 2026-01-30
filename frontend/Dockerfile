# Multi-stage Dockerfile for MARINA frontend
# Stage 1: Build
FROM node:20-slim as build

WORKDIR /app

# Copy package files
COPY package*.json ./
COPY tsconfig*.json ./
COPY vite.config.ts ./

# Install dependencies
RUN npm ci

# Copy source code
COPY src/ ./src/
COPY index.html ./

# Set production API base URLs for build
ARG VITE_API_BASE=/api
ARG VITE_API_BASE_MARINA=/api/marina
ARG VITE_API_BASE_SPECTRE=/api/spectre
ENV VITE_API_BASE=${VITE_API_BASE}
ENV VITE_API_BASE_MARINA=${VITE_API_BASE_MARINA}
ENV VITE_API_BASE_SPECTRE=${VITE_API_BASE_SPECTRE}

# Build
RUN npm run build

# Stage 2: Copy dist to volume mount point
# The dist directory will be shared with nginx via docker volume
FROM alpine:latest
RUN apk add --no-cache bash
COPY --from=build /app/dist /dist
# Keep container running so volume stays mounted
CMD ["sh", "-c", "while true; do sleep 3600; done"]

