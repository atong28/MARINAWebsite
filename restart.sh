#!/bin/bash
set -euo pipefail

PROJECT_ROOT="/Users/atong/Documents/UCSD/Projects/GURU/MARINABackend"
COMPOSE="docker compose"

cd "$PROJECT_ROOT"

echo "Stopping containers and removing volumes..."
$COMPOSE down -v --remove-orphans

echo "Rebuilding images and starting services..."
$COMPOSE up -d --build

echo "Done."
