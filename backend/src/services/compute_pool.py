import asyncio
import logging
import multiprocessing as mp
import os
import signal
import threading
import traceback
from dataclasses import dataclass
from typing import Any, Callable, Optional
import time

logger = logging.getLogger(__name__)


class ComputeOverloadedError(RuntimeError):
    """Raised when the compute pool is at capacity."""


class ComputeTimeoutError(TimeoutError):
    """Raised when a compute job times out."""


@dataclass
class ComputeResult:
    ok: bool
    value: Any = None
    error: Optional[str] = None
    trace: Optional[str] = None


def _worker_entry(func: Callable, args: tuple, kwargs: dict, result_queue: mp.Queue) -> None:
    try:
        value = func(*args, **kwargs)
        result_queue.put(ComputeResult(ok=True, value=value))
    except Exception as exc:  # pragma: no cover - best-effort capture
        result_queue.put(
            ComputeResult(
                ok=False,
                error=f"{type(exc).__name__}: {exc}",
                trace=traceback.format_exc(),
            )
        )


def _wait_for_result(result_queue: mp.Queue, process: mp.Process) -> ComputeResult:
    """Wait for a result or detect worker exit without a result."""
    while True:
        try:
            return result_queue.get_nowait()
        except Exception:
            if not process.is_alive():
                process.join(timeout=0.1)
                return ComputeResult(ok=False, error="Worker exited without result")
            time.sleep(0.01)


class ComputeWorkerPool:
    def __init__(self, max_workers: int, max_queue: int = 0) -> None:
        if max_workers <= 0:
            raise ValueError("max_workers must be > 0")
        self._ctx = mp.get_context("spawn")
        self._semaphore = asyncio.Semaphore(max_workers)
        self._max_queue = max_queue
        self._pending = 0
        self._pending_lock = asyncio.Lock()
        self._active_lock = asyncio.Lock()
        self._active: set[mp.Process] = set()

    async def _acquire_slot(self) -> None:
        if self._max_queue <= 0:
            if self._semaphore.locked():
                raise ComputeOverloadedError("Compute pool is at capacity")
            await self._semaphore.acquire()
            return
        async with self._pending_lock:
            if self._pending >= self._max_queue:
                raise ComputeOverloadedError("Compute queue is full")
            self._pending += 1
        await self._semaphore.acquire()
        async with self._pending_lock:
            self._pending = max(0, self._pending - 1)

    @staticmethod
    def _terminate_process(process: mp.Process) -> None:
        if process.is_alive():
            process.terminate()
            process.join(timeout=1.0)
        if process.is_alive() and process.pid:
            os.kill(process.pid, signal.SIGKILL)

    async def run(
        self,
        func: Callable,
        *args: Any,
        timeout: Optional[float] = None,
        **kwargs: Any,
    ) -> Any:
        await self._acquire_slot()
        result_queue = self._ctx.Queue()
        process = self._ctx.Process(
            target=_worker_entry,
            args=(func, args, kwargs, result_queue),
        )
        process.start()

        async with self._active_lock:
            self._active.add(process)

        try:
            wait_task = asyncio.to_thread(_wait_for_result, result_queue, process)
            result = await asyncio.wait_for(wait_task, timeout=timeout)
        except asyncio.TimeoutError as exc:
            self._terminate_process(process)
            raise ComputeTimeoutError("Compute job timed out") from exc
        except asyncio.CancelledError:
            self._terminate_process(process)
            raise
        finally:
            self._semaphore.release()
            async with self._active_lock:
                self._active.discard(process)
            try:
                result_queue.close()
            except Exception:
                pass

        if not result.ok:
            logger.error("Compute worker error: %s\n%s", result.error, result.trace or "")
            raise RuntimeError(result.error or "Compute worker failed")

        return result.value

    async def shutdown(self) -> None:
        async with self._active_lock:
            processes = list(self._active)
            self._active.clear()
        for proc in processes:
            self._terminate_process(proc)


_pool_lock = threading.Lock()
_pool: Optional[ComputeWorkerPool] = None


def get_compute_pool(max_workers: int, max_queue: int = 0) -> ComputeWorkerPool:
    global _pool
    with _pool_lock:
        if _pool is None:
            _pool = ComputeWorkerPool(max_workers=max_workers, max_queue=max_queue)
        return _pool


async def shutdown_compute_pool() -> None:
    global _pool
    with _pool_lock:
        pool = _pool
        _pool = None
    if pool is not None:
        await pool.shutdown()
