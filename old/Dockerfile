FROM ghcr.io/prefix-dev/pixi:latest

# working directory
WORKDIR /app

# copy only what's needed
COPY pixi.toml /app/pixi.toml
COPY pixi.lock /app/pixi.lock
COPY requirements.txt /app/requirements.txt

# install runtime deps via pip (pixi image should already contain Python)
RUN if command -v pixi >/dev/null 2>&1; then \
			echo "pixi found, installing via pixi..." && pixi install || true; \
		else \
			pip install --no-cache-dir -r /app/requirements.txt; \
		fi

# Copy application code (exclude ckpt, large files; use .dockerignore)
COPY src /app/src
COPY app.py /app/app.py

# copy frontend assets so Flask can serve them
COPY static /app/static

# copy startup script that prefers running inside pixi
COPY start.sh /app/start.sh
RUN chmod +x /app/start.sh

# Default data path is mounted at runtime; do not COPY data into image
ENV DATA_DIR=/app/data
EXPOSE 5000

ENV PORT=5000

# Use start.sh which will exec gunicorn inside pixi if available
ENTRYPOINT ["/app/start.sh"]
