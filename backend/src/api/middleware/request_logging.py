import logging
import time
from typing import Optional

from fastapi import Request

from src.api.middleware.request_id import get_request_id
from src.services.request_log_store import get_log_store


class RequestIdFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        request_id = get_request_id()
        record.request_id = request_id or "-"
        return True


class RequestLogHandler(logging.Handler):
    def __init__(self) -> None:
        super().__init__()
        self._store = get_log_store()

    def emit(self, record: logging.LogRecord) -> None:
        request_id = getattr(record, "request_id", None)
        if not request_id or request_id == "-":
            return
        entry = {
            "ts": time.time(),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
        }
        self._store.add(request_id, entry)


def setup_request_logging(app) -> None:
    root_logger = logging.getLogger()
    root_logger.addFilter(RequestIdFilter())
    root_logger.addHandler(RequestLogHandler())

    @app.middleware("http")
    async def log_request(request: Request, call_next):
        start = time.time()
        logger = logging.getLogger("request")
        request_id = getattr(request.state, "request_id", None)
        logger.info(
            "request_start method=%s path=%s request_id=%s",
            request.method,
            request.url.path,
            request_id,
        )
        response = None
        try:
            response = await call_next(request)
        finally:
            duration_ms = int((time.time() - start) * 1000)
            logger.info(
                "request_end method=%s path=%s status=%s duration_ms=%s request_id=%s",
                request.method,
                request.url.path,
                getattr(response, "status_code", "unknown"),
                duration_ms,
                request_id,
            )
        return response
