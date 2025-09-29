#!/usr/bin/env sh
set -eu

PORT="${PORT:-5000}"

# Ensure logs directory exists
mkdir -p /app/logs

echo "Starting server on port ${PORT}"

GUNICORN_ARGS="-w 1 -b 0.0.0.0:${PORT} \
  --access-logfile /app/logs/gunicorn.access.log \
  --error-logfile /app/logs/gunicorn.error.log \
  --log-level warning"

if command -v pixi >/dev/null 2>&1; then
  echo "pixi detected, running inside pixi"
  if pixi help >/dev/null 2>&1 && pixi help | grep -q "run"; then
    exec pixi run gunicorn $GUNICORN_ARGS app:app
  else
    exec pixi shell -c "gunicorn $GUNICORN_ARGS app:app"
  fi
else
  echo "pixi not found; running command directly"
  exec gunicorn $GUNICORN_ARGS app:app
fi
