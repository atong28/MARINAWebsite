#!/usr/bin/env sh
set -eu

PORT="${PORT:-5000}"
# Ensure logs directory exists
mkdir -p /app/logs

# Gunicorn logging: access and error logs to files in logs/
CMD="gunicorn -w 1 -b 0.0.0.0:${PORT} \
  --access-logfile /app/logs/gunicorn.access.log \
  --error-logfile /app/logs/gunicorn.error.log \
  --log-level warning \
  app:app"

echo "Starting server on port ${PORT}"

# If pixi is present, prefer to run inside the pixi environment so installed binaries are available
if command -v pixi >/dev/null 2>&1; then
  echo "pixi detected, running inside pixi"
  # prefer pixi run if available, otherwise fallback to pixi shell
  if pixi help >/dev/null 2>&1 && pixi help | grep -q "run"; then
    exec pixi run -- sh -c "$CMD"
  else
    exec pixi shell -c "$CMD"
  fi
else
  echo "pixi not found; running command directly"
  exec sh -c "$CMD"
fi
