"""
SMILES search endpoint.
"""
import logging

import torch
from fastapi import APIRouter, HTTPException, Request, status

from src.api.middleware.rate_limit import get_limiter
from src.config import MAX_TOP_K, MOLECULE_IMG_SIZE
from src.domain.ranker import RankingSet, filter_rankingset_rows
from src.domain.models.prediction_result import SmilesSearchResponse
from src.domain.models.spectral_data import SmilesSearchRequest
from src.services.model_manifest import resolve_and_validate_model_id
from src.services.model_service import ModelService
from src.services.molecule_renderer import MoleculeRenderer
from src.services.retrieval_helpers import (
    build_result_cards,
    kept_indices_for_mw_range,
)

logger = logging.getLogger(__name__)
router = APIRouter()
limiter = get_limiter()


@router.post("/smiles-search", response_model=SmilesSearchResponse, status_code=status.HTTP_200_OK)
@limiter.limit("20 per minute")
async def smiles_search(request: Request, data: SmilesSearchRequest):
    """Search for similar molecules using a SMILES string."""
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

    requested_k = data.k
    if requested_k > MAX_TOP_K:
        requested_k = MAX_TOP_K

    try:
        smiles = data.smiles.strip()
        if not smiles:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="SMILES string cannot be empty",
            )

        mw_min, mw_max = data.mw_min, data.mw_max
        if mw_min is not None and mw_max is not None and mw_min > mw_max:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="mw_min cannot be greater than mw_max",
            )

        fp_loader = model_service.get_fp_loader(mid)
        if fp_loader is None:
            raise HTTPException(
                status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
                detail="Model not loaded. Ensure the model is initialized first.",
            )
        fp = fp_loader.build_mfp_for_smiles(smiles)
        if fp is None:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Invalid SMILES or unable to generate fingerprint",
            )

        query_tensor = torch.tensor(fp, dtype=torch.float32)
        session = model_service.get_session(mid)
        rankingset_full = model_service.get_rankingset(mid)
        if rankingset_full is None:
            raise HTTPException(
                status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
                detail="Rankingset not available",
            )

        kept_indices = kept_indices_for_mw_range(session, mw_min, mw_max)
        if not kept_indices:
            return SmilesSearchResponse(
                results=[],
                total_count=0,
                offset=0,
                limit=0,
                query_smiles=smiles,
                query_fp=query_tensor.tolist(),
            )

        filtered_store = filter_rankingset_rows(rankingset_full, kept_indices)
        ranker = RankingSet(store=filtered_store, metric="cosine")
        n = min(requested_k, len(kept_indices))
        sims, local_idxs = ranker.retrieve_with_scores(query_tensor.unsqueeze(0), n=n)
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
        molecule_renderer = MoleculeRenderer.instance()
        results = build_result_cards(
            session,
            rankingset_full,
            pairs,
            query_tensor,
            fp_loader,
            molecule_renderer,
            img_size=MOLECULE_IMG_SIZE,
            max_cards=requested_k,
            include_plain_svg=True,
        )

        return SmilesSearchResponse(
            results=results,
            total_count=len(results),
            offset=0,
            limit=len(results),
            query_smiles=smiles,
            query_fp=query_tensor.tolist(),
        )
    except HTTPException:
        raise
    except Exception as e:
        logger.error("SMILES search error: %s", e, exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=str(e),
        )
