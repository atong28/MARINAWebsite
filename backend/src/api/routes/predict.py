"""
Prediction endpoint for spectral data.
"""
import logging

import torch
from fastapi import APIRouter, HTTPException, Request, status

from src.api.app import run_heavy
from src.api.middleware.rate_limit import get_limiter
from src.config import MAX_TOP_K, MOLECULE_IMG_SIZE, PREDICT_TIMEOUT_S
from src.domain.models.prediction_result import PredictResponse
from src.domain.models.spectral_data import PredictRequest
from src.services.model_manifest import resolve_and_validate_model_id
from src.services.model_service import ModelService
from src.services.molecule_renderer import MoleculeRenderer
from src.services.compute_pool import ComputeOverloadedError, ComputeTimeoutError
from src.services.retrieval_helpers import (
    build_result_cards,
    kept_indices_for_mw_range,
)

logger = logging.getLogger(__name__)
router = APIRouter()
limiter = get_limiter()


@router.post("/predict", response_model=PredictResponse, status_code=status.HTTP_200_OK)
@limiter.limit("30 per minute")
async def predict(request: Request, data: PredictRequest):
    """
    Predict molecular structures from spectral data.
    """
    mid, err = resolve_and_validate_model_id(data.model_id)
    if err is not None:
        code, detail = err
        raise HTTPException(status_code=code, detail=detail)

    model_service = ModelService.instance()
    model_service.ensure_loaded(mid)

    requested_k = data.k
    if requested_k > MAX_TOP_K:
        requested_k = MAX_TOP_K

    try:
        mw_min, mw_max = data.mw_min, data.mw_max
        if mw_min is not None and mw_max is not None and mw_min > mw_max:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="mw_min cannot be greater than mw_max",
            )

        raw_data = {}
        if data.raw.hsqc:
            raw_data["hsqc"] = data.raw.hsqc
        if data.raw.h_nmr:
            raw_data["h_nmr"] = data.raw.h_nmr
        if data.raw.c_nmr:
            raw_data["c_nmr"] = data.raw.c_nmr
        if data.raw.mass_spec:
            raw_data["mass_spec"] = data.raw.mass_spec
        if data.raw.mw:
            raw_data["mw"] = data.raw.mw

        pool = request.app.state.compute_pool
        out = await run_heavy(
            request,
            pool.run(
                "predict",
                {"raw_inputs": raw_data, "k": MAX_TOP_K, "model_id": mid},
                timeout=PREDICT_TIMEOUT_S,
            ),
        )
        if not isinstance(out, tuple) or len(out) < 2:
            logger.error("Unexpected predictor output: %s", out)
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail="Unexpected predictor output",
            )
        pred_prob = out[2] if len(out) > 2 else None
        if pred_prob is None:
            logger.error("Predictor did not return pred_prob as third element")
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail="Predictor did not return a fingerprint vector",
            )

        session = model_service.get_session(mid)
        kept_indices = kept_indices_for_mw_range(session, mw_min, mw_max)
        if not kept_indices:
            return PredictResponse(
                results=[],
                total_count=0,
                offset=0,
                limit=0,
                pred_fp=pred_prob,
            )

        ranker = session.get_filtered_rankingset(mw_min, mw_max)
        pred_tensor = torch.tensor(pred_prob, dtype=torch.float32)
        n = min(requested_k, len(kept_indices))
        sims, local_idxs = ranker.retrieve_with_scores(pred_tensor.unsqueeze(0), n=n)
        sims = sims.squeeze()
        local_idxs = local_idxs.squeeze()
        sim_list = sims.tolist()
        idx_list = local_idxs.tolist()
        if not isinstance(sim_list, list):
            sim_list = [sim_list]
        if not isinstance(idx_list, list):
            idx_list = [idx_list]
        pairs = [
            (int(kept_indices[idx_list[i]]), float(sim_list[i]))
            for i in range(len(sim_list))
        ]
        fp_loader = model_service.get_fp_loader(mid)
        molecule_renderer = MoleculeRenderer.instance()
        # Use the full rankingset tensor for build_result_cards
        rankingset_full = session.get_rankingset().data
        results = build_result_cards(
            session,
            rankingset_full,
            pairs,
            pred_tensor,
            fp_loader,
            molecule_renderer,
            img_size=MOLECULE_IMG_SIZE,
            max_cards=requested_k,
            include_plain_svg=True,
        )

        return PredictResponse(
            results=results,
            total_count=len(results),
            offset=0,
            limit=len(results),
            pred_fp=pred_prob,
        )
    except HTTPException:
        raise
    except ComputeOverloadedError:
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Server is busy. Please retry shortly.",
        )
    except ComputeTimeoutError:
        raise HTTPException(
            status_code=status.HTTP_504_GATEWAY_TIMEOUT,
            detail="Prediction timed out.",
        )
    except Exception as e:
        logger.error("Prediction error: %s", e, exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=str(e),
        )
