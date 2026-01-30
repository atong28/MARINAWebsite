import asyncio
import logging
import multiprocessing as mp
import os
import signal
import threading
import traceback
import time
from dataclasses import dataclass
from typing import Any, Optional


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


def _worker_loop(request_queue: mp.Queue, response_queue: mp.Queue) -> None:
    """Worker process: preload model once, then handle jobs."""
    from src.services.model_service import ModelService
    from src.domain.predictor import predict_from_raw
    import torch
    from src.domain.drawing.draw import (
        compute_bit_environments_batch,
        draw_similarity_comparison,
        render_molecule_with_change_overlays,
    )

    model_service = ModelService.instance()
    # Preload default model once in this worker process
    model_service.preload_resources()
    fp_loader = model_service.get_fp_loader()

    while True:
        job = request_queue.get()
        if job is None:
            break
        job_id = job["job_id"]
        op = job["op"]
        payload = job["payload"]
        try:
            if op == "predict":
                value = predict_from_raw(
                    raw_inputs=payload["raw_inputs"],
                    k=payload["k"],
                    model_id=payload.get("model_id"),
                )
            elif op == "bit_envs":
                value = compute_bit_environments_batch(
                    smiles=payload["smiles"],
                    fp_indices=payload["fp_indices"],
                    fp_loader=fp_loader,
                )
            elif op == "similarity_map":
                fp_tensor = torch.tensor(payload["predicted_fp"], dtype=torch.float32)
                value = draw_similarity_comparison(
                    predicted_fp=fp_tensor,
                    retrieval_smiles=payload["smiles"],
                    fp_loader=fp_loader,
                    img_size=payload["img_size"],
                )
            elif op == "change_overlay":
                value = render_molecule_with_change_overlays(
                    smiles=payload["smiles"],
                    original_fp=payload["original_fp"],
                    new_fp=payload["new_fp"],
                    fp_loader=fp_loader,
                    threshold=payload.get("threshold", 0.5),
                    img_size=payload["img_size"],
                )
            else:
                raise ValueError(f"Unknown op '{op}'")

            response_queue.put((job_id, ComputeResult(ok=True, value=value)))
        except Exception as exc:  # pragma: no cover
            response_queue.put(
                (
                    job_id,
                    ComputeResult(
                        ok=False,
                        error=f"{type(exc).__name__}: {exc}",
                        trace=traceback.format_exc(),
                    ),
                )
            )


class ComputeWorkerPool:
    """Persistent worker pool that preloads model resources once per worker."""

    def __init__(self, max_workers: int, max_queue: int = 0) -> None:
        if max_workers <= 0:
            raise ValueError("max_workers must be > 0")

        self._ctx = mp.get_context("spawn")
        self._request_queue = self._ctx.Queue()
        self._response_queue = self._ctx.Queue()
        self._workers: list[mp.Process] = []

        self._max_queue = max_queue
        self._pending = 0
        self._pending_lock = asyncio.Lock()

        self._jobs: dict[str, asyncio.Future] = {}
        self._jobs_lock = asyncio.Lock()
        self._response_task: Optional[asyncio.Task] = None
        self._restart_lock = asyncio.Lock()

        self._spawn_workers(max_workers)

    def _spawn_workers(self, count: int) -> None:
        for _ in range(count):
            proc = self._ctx.Process(
                target=_worker_loop,
                args=(self._request_queue, self._response_queue),
            )
            proc.start()
            self._workers.append(proc)

    async def _ensure_response_task(self) -> None:
        if self._response_task is None:
            self._response_task = asyncio.create_task(self._response_loop())

    async def _response_loop(self) -> None:
        """Background task that delivers worker results to waiting Futures."""
        while True:
            try:
                job_id, result = await asyncio.to_thread(self._response_queue.get)
            except Exception:
                continue

            async with self._jobs_lock:
                future = self._jobs.pop(job_id, None)

            if future is None:
                continue

            if result.ok:
                if not future.done():
                    future.set_result(result.value)
            else:
                logger.error("Compute worker error: %s\n%s", result.error, result.trace or "")
                if not future.done():
                    future.set_exception(RuntimeError(result.error or "Compute worker failed"))

    async def _acquire_slot(self) -> None:
        """Apply simple queue/backpressure limits."""
        if self._max_queue <= 0:
            if self._pending >= len(self._workers):
                raise ComputeOverloadedError("Compute pool is at capacity")
            self._pending += 1
            return

        async with self._pending_lock:
            if self._pending >= self._max_queue:
                raise ComputeOverloadedError("Compute queue is full")
            self._pending += 1

    async def _restart_workers(self) -> None:
        """Kill and respawn all workers (used on timeout)."""
        async with self._restart_lock:
            for _ in self._workers:
                self._request_queue.put(None)

            for proc in self._workers:
                if proc.is_alive():
                    proc.terminate()
                    proc.join(timeout=1.0)
                if proc.is_alive() and proc.pid:
                    os.kill(proc.pid, signal.SIGKILL)

            count = len(self._workers)
            self._workers = []
            self._spawn_workers(count)

    async def run(self, op: str, payload: dict, timeout: Optional[float] = None) -> Any:
        """Submit a job to the pool and await the result."""
        await self._ensure_response_task()
        await self._acquire_slot()

        job_id = f"{int(time.time() * 1000)}-{os.getpid()}-{id(payload)}"
        loop = asyncio.get_running_loop()
        future: asyncio.Future = loop.create_future()

        async with self._jobs_lock:
            self._jobs[job_id] = future

        self._request_queue.put({"job_id": job_id, "op": op, "payload": payload})

        try:
            return await asyncio.wait_for(future, timeout=timeout)
        except asyncio.TimeoutError as exc:
            # Drop the future (if still present) and restart workers
            async with self._jobs_lock:
                self._jobs.pop(job_id, None)
            await self._restart_workers()
            raise ComputeTimeoutError("Compute job timed out") from exc
        finally:
            async with self._pending_lock:
                self._pending = max(0, self._pending - 1)

    async def shutdown(self) -> None:
        """Terminate all workers and stop background response task."""
        for _ in self._workers:
            self._request_queue.put(None)

        for proc in self._workers:
            if proc.is_alive():
                proc.terminate()
                proc.join(timeout=1.0)
            if proc.is_alive() and proc.pid:
                os.kill(proc.pid, signal.SIGKILL)

        if self._response_task is not None:
            self._response_task.cancel()


_pool_lock = threading.Lock()
_pool: Optional[ComputeWorkerPool] = None


def get_compute_pool(max_workers: int, max_queue: int = 0) -> ComputeWorkerPool:
    """Singleton accessor for ComputeWorkerPool."""
    global _pool
    with _pool_lock:
        if _pool is None:
            _pool = ComputeWorkerPool(max_workers=max_workers, max_queue=max_queue)
        return _pool


async def shutdown_compute_pool() -> None:
    """Shutdown the global compute pool if it exists."""
    global _pool
    with _pool_lock:
        pool = _pool
        _pool = None
    if pool is not None:
        await pool.shutdown()