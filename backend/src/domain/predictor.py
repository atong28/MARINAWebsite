import json
import logging
import threading
from typing import Dict, Any

import torch
import torch.nn.functional as F

from src.domain.data import collate
from src.domain.model_session import ModelSession

logger = logging.getLogger(__name__)
_session_lock = threading.Lock()
_session: ModelSession | None = None


def get_session() -> ModelSession:
    """
    Get or create the singleton ModelSession for this process.
    """
    global _session
    if _session is not None:
        return _session
    with _session_lock:
        if _session is None:
            _session = ModelSession.from_paths()
    return _session


def load_model(ckpt_path: str = "data/best.ckpt", params_path: str = "data/params.json") -> torch.nn.Module:
    """
    Backward-compatible shim that returns the underlying model.
    The path arguments are currently ignored in favor of central config.
    """
    session = get_session()
    return session.model

def _preprocess_raw_inputs(raw_inputs: Dict[str, Any]) -> Dict[str, torch.Tensor]:
    """
    Convert raw input dict into the tensor dict expected by the model.
    """
    processed: Dict[str, torch.Tensor] = {}
    for k_mod, v in raw_inputs.items():
        print(k_mod, v)
        if isinstance(v, list):
            tensor = torch.tensor(v, dtype=torch.float)
            # Reshape sequence data to 2D (num_peaks, features_per_peak)
            if k_mod in ["hsqc", "h_nmr", "c_nmr", "mass_spec"]:
                if k_mod == "hsqc":
                    # HSQC: [H1, C13, intensity, ...] -> [C13, H1, intensity, ...]
                    if len(v) % 3 == 0:
                        tensor = tensor.view(-1, 3)
                        tensor = tensor[:, [1, 0, 2]]
                    else:
                        logger.warning(
                            "HSQC data length %d is not divisible by 3, skipping", len(v)
                        )
                        continue
                elif k_mod == "c_nmr":
                    tensor = tensor.view(-1, 1)
                    # tensor = F.pad(tensor, (0, 2), "constant", 0)
                elif k_mod == "h_nmr":
                    tensor = tensor.view(-1, 1)
                    # tensor = F.pad(tensor, (1, 1), "constant", 0)
                elif k_mod == "mass_spec":
                    if len(v) % 2 == 0:
                        tensor = tensor.view(-1, 2)
                        # tensor = F.pad(tensor, (0, 1), "constant", 0)
                    else:
                        logger.warning(
                            "Mass spec data length %d is not divisible by 2, skipping",
                            len(v),
                        )
                        continue
            processed[k_mod] = tensor
        elif k_mod == "mw":
            processed[k_mod] = torch.tensor([v], dtype=torch.float).view(1, 1)
        elif isinstance(v, torch.Tensor):
            processed[k_mod] = v
        else:
            try:
                processed[k_mod] = torch.tensor(v, dtype=torch.float)
            except Exception:
                # leave non-tensor values out
                continue
    return processed


def _run_model(session: ModelSession, processed: Dict[str, torch.Tensor]) -> torch.Tensor:
    """
    Run the MARINA model on preprocessed inputs and return the first prediction vector.
    """
    batch_inputs, _ = collate(
        [(processed, torch.zeros((session.args.out_dim,), dtype=torch.float))]
    )

    with _session_lock, torch.no_grad():
        logger.info("Calling model forward...")
        out = session.model(batch_inputs)
        logger.info(f"Model output shape: {out.shape}")
        preds = out.detach().cpu()
        return preds[0]


def _retrieve_top_k(session: ModelSession, pred_fp: torch.Tensor, k: int):
    """
    Apply sigmoid to the prediction, then retrieve top-k scores and indices.
    """
    try:
        k_int = int(k)
    except (ValueError, TypeError) as e:
        logger.error("Failed to convert k=%r (type: %s) to int: %s", k, type(k), e)
        k_int = 5

    pred = torch.sigmoid(pred_fp)
    logger.info("Pred FP sum: %s", torch.where(pred.squeeze() > 0.5, 1, 0).sum())

    ranker = session.get_rankingset()
    sims, idxs = ranker.retrieve_with_scores(pred.unsqueeze(0), n=k_int)

    # sims and idxs are already top-k, sorted, and aligned
    return sims.tolist(), idxs.tolist(), pred.tolist()


def predict_from_raw(raw_inputs: Dict[str, Any], k: int = 5):
    """Accepts a raw input dict (same spec as SpectralInputLoader.load output) and returns top-k."""
    logger.info("predict_from_raw called with k=%r (type: %s)", k, type(k))
    logger.info("raw_inputs keys: %s", list(raw_inputs.keys()))

    session = get_session()

    processed = _preprocess_raw_inputs(raw_inputs)
    logger.info("Processed inputs: %s", list(processed.keys()))

    pred_fp = _run_model(session, processed)
    vals, idxs, pred_prob = _retrieve_top_k(session, pred_fp, k)

    return vals, idxs, pred_prob
