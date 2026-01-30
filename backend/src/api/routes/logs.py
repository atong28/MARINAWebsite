import asyncio
import json

from fastapi import APIRouter, Request
from fastapi.responses import StreamingResponse

from src.services.request_log_store import get_log_store

router = APIRouter()


@router.get("/logs/stream/{request_id}")
async def stream_logs(request_id: str, request: Request):
    store = get_log_store()
    queue, initial = store.subscribe(request_id)

    async def event_stream():
        for entry in initial:
            yield f"data: {json.dumps(entry)}\n\n"
        try:
            while True:
                if await request.is_disconnected():
                    break
                try:
                    entry = await asyncio.wait_for(queue.get(), timeout=10)
                    yield f"data: {json.dumps(entry)}\n\n"
                except asyncio.TimeoutError:
                    yield "event: ping\ndata: {}\n\n"
        finally:
            store.unsubscribe(request_id, queue)

    return StreamingResponse(event_stream(), media_type="text/event-stream")
