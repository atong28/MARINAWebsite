import uuid
from contextvars import ContextVar
from typing import Optional

from fastapi import Request

request_id_var: ContextVar[Optional[str]] = ContextVar("request_id", default=None)


def get_request_id() -> Optional[str]:
    return request_id_var.get()


def setup_request_id_middleware(app) -> None:
    @app.middleware("http")
    async def add_request_id(request: Request, call_next):
        incoming = request.headers.get("X-Request-ID")
        request_id = incoming or str(uuid.uuid4())
        request.state.request_id = request_id
        token = request_id_var.set(request_id)
        try:
            response = await call_next(request)
        finally:
            request_id_var.reset(token)
        response.headers["X-Request-ID"] = request_id
        return response
