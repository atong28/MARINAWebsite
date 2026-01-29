"""Ablation analysis endpoint."""

import logging
from typing import Dict, Any

import torch
from fastapi import APIRouter, HTTPException, status, Request

from src.api.middleware.rate_limit import get_limiter
from src.config import MOLECULE_IMG_SIZE
from src.domain.drawing.draw import (
    compute_bit_environments_batch,
    draw_similarity_comparison,
    render_molecule_with_change_overlays,
    pil_image_to_base64,
)
from src.domain.models.analysis_result import AblationRequest, AblationResponse
from src.domain.predictor import predict_from_raw
from src.services.model_manifest import resolve_and_validate_model_id
from src.services.model_service import ModelService

logger = logging.getLogger(__name__)
router = APIRouter()

limiter = get_limiter()


def _build_raw_input(raw: AblationRequest) -> Dict[str, Any]:
    """Convert spectral input to raw dict expected by predictor."""
    raw_data: Dict[str, Any] = {}

    if raw.raw.hsqc:
        raw_data["hsqc"] = raw.raw.hsqc
    if raw.raw.h_nmr:
        raw_data["h_nmr"] = raw.raw.h_nmr
    if raw.raw.c_nmr:
        raw_data["c_nmr"] = raw.raw.c_nmr
    if raw.raw.mass_spec:
        raw_data["mass_spec"] = raw.raw.mass_spec
    if raw.raw.mw is not None:
        raw_data["mw"] = raw.raw.mw

    return raw_data


@router.post("/ablation", response_model=AblationResponse, status_code=status.HTTP_200_OK)
@limiter.limit("30 per minute")
async def run_ablation(request: Request, data: AblationRequest) -> AblationResponse:
    """Run an ablation prediction using modified spectral data."""
    mid, err = resolve_and_validate_model_id(data.model_id)
    if err is not None:
        code, detail = err
        raise HTTPException(status_code=code, detail=detail)

    model_service = ModelService.instance()
    # Auto-load model if not ready
    if not model_service.is_ready(mid):
        logger.info("Model %s not loaded, attempting to load...", mid)
        try:
            model_service.ensure_loaded(mid)
            if not model_service.is_ready(mid):
                raise HTTPException(
                    status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
                    detail=f"Model '{mid}' failed to load. Please check model files and try again.",
                )
            logger.info("Model %s loaded successfully", mid)
        except HTTPException:
            raise
        except Exception as e:
            logger.error("Failed to auto-load model %s: %s", mid, e, exc_info=True)
            raise HTTPException(
                status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
                detail=f"Model '{mid}' is not available and failed to load: {str(e)}",
            )

    raw_data = _build_raw_input(data)
    if not raw_data:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="At least one spectral data field must be provided for ablation.",
        )

    try:
        prediction_output = predict_from_raw(raw_data, k=1, model_id=mid)
    except Exception as exc:
        logger.error("Ablation prediction failed: %s", exc, exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to run ablation prediction.",
        ) from exc

    if not isinstance(prediction_output, tuple) or len(prediction_output) < 3:
        logger.error("Unexpected predictor output for ablation: %s", prediction_output)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Unexpected predictor output.",
        )

    pred_fp = prediction_output[2]
    if pred_fp is None:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Predictor did not return a fingerprint vector.",
        )

    if isinstance(pred_fp, torch.Tensor):
        pred_fp_list = pred_fp.flatten().tolist()
    else:
        pred_fp_list = list(pred_fp)

    fp_loader = model_service.get_fp_loader(mid)
    if fp_loader is None:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Fingerprint loader is not available.",
        )

    # Compute similarity map for the provided molecule
    similarity_map_encoded = None
    try:
        fp_tensor = torch.tensor(pred_fp_list, dtype=torch.float32)
        similarity_img = draw_similarity_comparison(fp_tensor, data.smiles, fp_loader, img_size=MOLECULE_IMG_SIZE)
        similarity_map_encoded = pil_image_to_base64(similarity_img)
    except Exception as exc:
        logger.warning("Failed to generate similarity map for ablation: %s", exc)

    # Determine active bits and compute environments
    threshold = data.bit_threshold
    active_bit_indices = [idx for idx, value in enumerate(pred_fp_list) if value is not None and value > threshold]

    if data.max_bits and len(active_bit_indices) > data.max_bits:
        active_bit_indices = active_bit_indices[: data.max_bits]

    bit_environments: Dict[int, Dict[str, Any]] = {}
    if active_bit_indices:
        try:
            bit_environments = compute_bit_environments_batch(data.smiles, active_bit_indices, fp_loader)
        except Exception as exc:
            logger.warning("Failed to compute bit environments for ablation: %s", exc)

    change_overlay_svg = None
    if data.reference_fp is not None:
        try:
            change_overlay_svg = render_molecule_with_change_overlays(
                data.smiles,
                data.reference_fp,
                pred_fp_list,
                fp_loader,
                threshold=data.bit_threshold,
                img_size=MOLECULE_IMG_SIZE,
            )
        except Exception as exc:
            logger.warning("Failed to generate change overlay SVG for ablation: %s", exc)

    return AblationResponse(
        pred_fp=pred_fp_list,
        active_bit_indices=active_bit_indices,
        similarity_map=similarity_map_encoded,
        bit_environments=bit_environments,
        change_overlay_svg=change_overlay_svg,
    )


