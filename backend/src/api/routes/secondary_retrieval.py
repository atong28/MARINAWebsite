"""
Secondary retrieval endpoint using difference fingerprint.
"""
import logging

import torch
from fastapi import APIRouter, HTTPException, Request, status

from src.api.middleware.rate_limit import get_limiter
from src.config import MAX_TOP_K, MOLECULE_IMG_SIZE
from src.domain.models.analysis_result import (
    SecondaryRetrievalRequest,
    SecondaryRetrievalResponse,
)
from src.domain.ranker import RankingSet
from src.services.model_manifest import resolve_and_validate_model_id
from src.services.model_service import ModelService
from src.services.molecule_renderer import MoleculeRenderer
from src.services.retrieval_helpers import build_result_cards

logger = logging.getLogger(__name__)
router = APIRouter()
limiter = get_limiter()


@router.post(
    "/secondary-retrieval",
    response_model=SecondaryRetrievalResponse,
    status_code=status.HTTP_200_OK,
)
@limiter.limit("20 per minute")
async def secondary_retrieval(request: Request, data: SecondaryRetrievalRequest):
    """Secondary retrieval using difference fingerprint (predicted - overlap)."""
    mid, err = resolve_and_validate_model_id(data.model_id)
    if err is not None:
        code, detail = err
        raise HTTPException(status_code=code, detail=detail)

    model_service = ModelService.instance()
    if not model_service.is_ready(mid):
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Model is not ready. Try again in a moment.",
        )

    try:
        k = data.k
        if k > MAX_TOP_K:
            k = MAX_TOP_K

        predicted_fp = data.predicted_fp
        retrieved_fp = data.retrieved_fp
        if not predicted_fp or not retrieved_fp:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Missing required fields: predicted_fp, retrieved_fp",
            )

        predicted_tensor = torch.tensor(predicted_fp, dtype=torch.float32)
        retrieved_tensor = torch.tensor(retrieved_fp, dtype=torch.float32)
        predicted_tensor = torch.clamp(predicted_tensor, 0.0, 1.0)
        retrieved_tensor = torch.clamp(retrieved_tensor, 0.0, 1.0)
        overlap = torch.min(predicted_tensor, retrieved_tensor)
        difference_fp = predicted_tensor - overlap
        difference_fp = torch.clamp(difference_fp, 0.0, 1.0)

        if difference_fp.sum() < 1e-6:
            logger.warning("Difference fingerprint is all zeros")
            return SecondaryRetrievalResponse(
                results=[],
                total_count=0,
                difference_fp=difference_fp.tolist(),
            )

        rankingset = model_service.get_rankingset(mid)
        if rankingset is None:
            raise HTTPException(
                status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
                detail="Rankingset not available",
            )
        ranker = RankingSet(store=rankingset, metric="cosine")
        sims, idxs = ranker.retrieve_with_scores(difference_fp.unsqueeze(0), n=k)
        sims = sims.squeeze()
        idxs = idxs.squeeze()
        sim_list = sims.tolist()
        idx_list = idxs.tolist()
        if not isinstance(sim_list, list):
            sim_list = [sim_list]
        if not isinstance(idx_list, list):
            idx_list = [idx_list]
        pairs = [(int(idx_list[i]), float(sim_list[i])) for i in range(len(sim_list))]
        session = model_service.get_session(mid)
        fp_loader = model_service.get_fp_loader(mid)
        molecule_renderer = MoleculeRenderer.instance()
        results = build_result_cards(
            session,
            rankingset,
            pairs,
            difference_fp,
            fp_loader,
            molecule_renderer,
            img_size=MOLECULE_IMG_SIZE,
            max_cards=k,
            include_plain_svg=False,
        )

        return SecondaryRetrievalResponse(
            results=results,
            total_count=len(results),
            difference_fp=difference_fp.tolist(),
        )
    except HTTPException:
        raise
    except Exception as e:
        logger.error("Secondary retrieval error: %s", e, exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=str(e),
        )
