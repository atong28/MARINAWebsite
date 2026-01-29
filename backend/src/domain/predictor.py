import logging
import threading
from typing import Any, Dict

import torch

from src.domain.data import collate
from src.domain.model_session import ModelSession
from src.services.model_service import ModelService

logger = logging.getLogger(__name__)
_forward_lock = threading.Lock()


def get_session(model_id: str | None = None) -> ModelSession:
    """Get ModelSession for model_id (default model if None) via ModelService."""
    return ModelService.instance().get_session(model_id)


def load_model(
    ckpt_path: str = "data/best.ckpt",
    params_path: str = "data/params.json",
) -> torch.nn.Module:
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
                elif k_mod == "h_nmr":
                    tensor = tensor.view(-1, 1)
                elif k_mod == "mass_spec":
                    if len(v) % 2 == 0:
                        tensor = tensor.view(-1, 2)
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


def _to_spectre_batch(
    processed: Dict[str, torch.Tensor]
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Convert MARINA-style modality dict into SPECTRE's (inputs, type_indicator).

    Type indicators (aligned with SPECTREDataset._pad_and_stack_input):
        0: HSQC
        1: C NMR
        2: H NMR
        3: MW
        4: Mass Spectrometry
    """
    inputs = []
    type_indicators: list[int] = []

    # HSQC is already (N, 3) after preprocessing.
    hsqc = processed.get("hsqc")
    if hsqc is not None and hsqc.numel() > 0:
        inputs.append(hsqc)
        type_indicators.extend([0] * hsqc.shape[0])

    c_nmr = processed.get("c_nmr")
    if c_nmr is not None and c_nmr.numel() > 0:
        # (N,1) -> (N,3) by padding zeros to the right.
        c_nmr_padded = torch.nn.functional.pad(c_nmr, (0, 2), "constant", 0)
        inputs.append(c_nmr_padded)
        type_indicators.extend([1] * c_nmr_padded.shape[0])

    h_nmr = processed.get("h_nmr")
    if h_nmr is not None and h_nmr.numel() > 0:
        # (N,1) -> (N,3) by padding one zero on each side.
        h_nmr_padded = torch.nn.functional.pad(h_nmr, (1, 1), "constant", 0)
        inputs.append(h_nmr_padded)
        type_indicators.extend([2] * h_nmr_padded.shape[0])

    mass_spec = processed.get("mass_spec")
    if mass_spec is not None and mass_spec.numel() > 0:
        # (N,2) -> (N,3) by padding a trailing zero.
        mass_spec_padded = torch.nn.functional.pad(mass_spec, (0, 1), "constant", 0)
        inputs.append(mass_spec_padded)
        type_indicators.extend([4] * mass_spec_padded.shape[0])

    mw = processed.get("mw")
    if mw is not None and mw.numel() > 0:
        # Scalar -> [mw, 0, 0]
        mw_val = float(mw.view(-1)[0].item())
        mw_tensor = torch.tensor([[mw_val, 0.0, 0.0]], dtype=torch.float)
        inputs.append(mw_tensor)
        type_indicators.append(3)

    if not inputs:
        raise ValueError("No valid spectral inputs provided for SPECTRE model.")

    stacked = torch.vstack(inputs)  # (N_total, 3)
    type_indicator = torch.tensor(type_indicators, dtype=torch.long)  # (N_total,)
    # Add batch dimension: (1, N_total, 3) and (1, N_total)
    return stacked.unsqueeze(0), type_indicator.unsqueeze(0)


def _run_model(session: ModelSession, processed: Dict[str, torch.Tensor]) -> torch.Tensor:
    """
    Run the underlying model on preprocessed inputs and return the first prediction vector.
    Handles both MARINA (dict-based) and SPECTRE (stacked tensor + type_indicator).
    """
    with _forward_lock, torch.no_grad():
        logger.info("Calling model forward for type=%s...", getattr(session, "type", "marina"))

        if getattr(session, "type", "marina") == "spectre":
            inputs, type_indicator = _to_spectre_batch(processed)
            out = session.model(inputs, type_indicator)
        else:
            batch_inputs, _ = collate(
                [(processed, torch.zeros((session.args.out_dim,), dtype=torch.float))]
            )
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


def predict_from_raw(
    raw_inputs: Dict[str, Any],
    k: int = 5,
    model_id: str | None = None,
):
    """Accepts a raw input dict (same spec as SpectralInputLoader.load output) and returns top-k."""
    logger.info("predict_from_raw called with k=%r (type: %s)", k, type(k))
    logger.info("raw_inputs keys: %s", list(raw_inputs.keys()))

    session = get_session(model_id)

    processed = _preprocess_raw_inputs(raw_inputs)
    logger.info("Processed inputs: %s", list(processed.keys()))

    pred_fp = _run_model(session, processed)
    vals, idxs, pred_prob = _retrieve_top_k(session, pred_fp, k)

    return vals, idxs, pred_prob
