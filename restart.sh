#!/bin/bash
set -euo pipefail
COMPOSE="docker compose"
echo "Stopping containers and removing volumes..."
$COMPOSE down -v --remove-orphans
echo "Rebuilding images and starting services..."
$COMPOSE up -d --build
echo "Done."
