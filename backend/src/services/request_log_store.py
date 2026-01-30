import asyncio
import threading
import time
from collections import deque
from typing import Deque, Dict, List, Optional, Tuple


class RequestLogStore:
    def __init__(self, max_entries: int = 500, ttl_seconds: int = 3600) -> None:
        self._max_entries = max_entries
        self._ttl_seconds = ttl_seconds
        self._lock = threading.Lock()
        self._store: Dict[str, Tuple[float, Deque[dict]]] = {}
        self._subscribers: Dict[str, List[asyncio.Queue]] = {}

    def _prune_locked(self) -> None:
        now = time.time()
        expired = [
            request_id
            for request_id, (created_at, _) in self._store.items()
            if now - created_at > self._ttl_seconds
        ]
        for request_id in expired:
            self._store.pop(request_id, None)
            self._subscribers.pop(request_id, None)

    def add(self, request_id: str, entry: dict) -> None:
        if not request_id:
            return
        with self._lock:
            created_at, logs = self._store.get(request_id, (time.time(), deque()))
            logs.append(entry)
            while len(logs) > self._max_entries:
                logs.popleft()
            self._store[request_id] = (created_at, logs)
            self._prune_locked()
            queues = list(self._subscribers.get(request_id, []))

        for q in queues:
            try:
                q.put_nowait(entry)
            except Exception:
                pass

    def get(self, request_id: str) -> List[dict]:
        with self._lock:
            self._prune_locked()
            _, logs = self._store.get(request_id, (time.time(), deque()))
            return list(logs)

    def subscribe(self, request_id: str) -> Tuple[asyncio.Queue, List[dict]]:
        queue: asyncio.Queue = asyncio.Queue()
        with self._lock:
            created_at, logs = self._store.get(request_id, (time.time(), deque()))
            self._store[request_id] = (created_at, logs)
            self._subscribers.setdefault(request_id, []).append(queue)
            initial = list(logs)
        return queue, initial

    def unsubscribe(self, request_id: str, queue: asyncio.Queue) -> None:
        with self._lock:
            queues = self._subscribers.get(request_id, [])
            if queue in queues:
                queues.remove(queue)
            if not queues:
                self._subscribers.pop(request_id, None)


_store_lock = threading.Lock()
_store: Optional[RequestLogStore] = None


def get_log_store() -> RequestLogStore:
    global _store
    with _store_lock:
        if _store is None:
            _store = RequestLogStore()
        return _store
