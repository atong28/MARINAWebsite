"""
Secondary retrieval endpoint using difference fingerprint.
"""
import logging
import torch
from fastapi import APIRouter, HTTPException, status, Request

from src.domain.models.analysis_result import SecondaryRetrievalRequest, SecondaryRetrievalResponse
from src.services.model_service import ModelService
from src.services.metadata_service import MetadataService
from src.services.molecule_renderer import MoleculeRenderer
from src.services.result_builder import build_result_card
from src.domain.ranker import RankingSet
from src.config import MOLECULE_IMG_SIZE, MAX_TOP_K
from src.api.middleware.rate_limit import get_limiter

logger = logging.getLogger(__name__)
router = APIRouter()

# Get limiter instance
limiter = get_limiter()


@router.post("/secondary-retrieval", response_model=SecondaryRetrievalResponse, status_code=status.HTTP_200_OK)
@limiter.limit("20 per minute")
async def secondary_retrieval(request: Request, data: SecondaryRetrievalRequest):
    """Secondary retrieval endpoint using difference fingerprint (predicted - overlap)."""
    # Check if model is ready
    model_service = ModelService.instance()
    if not model_service.is_ready():
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Model is not ready. Please try again in a moment."
        )
    
    try:
        # Validate k
        k = data.k
        if k > MAX_TOP_K:
            k = MAX_TOP_K
        
        # Get fingerprints from request
        predicted_fp = data.predicted_fp
        retrieved_fp = data.retrieved_fp
        
        if not predicted_fp or not retrieved_fp:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Missing required fields: predicted_fp, retrieved_fp"
            )
        
        # Convert to tensors
        predicted_tensor = torch.tensor(predicted_fp, dtype=torch.float32)
        retrieved_tensor = torch.tensor(retrieved_fp, dtype=torch.float32)
        
        # Ensure fingerprints are in valid range [0, 1]
        predicted_tensor = torch.clamp(predicted_tensor, 0.0, 1.0)
        retrieved_tensor = torch.clamp(retrieved_tensor, 0.0, 1.0)
        
        # Compute overlap (element-wise minimum for binary/probabilistic fingerprints)
        overlap = torch.min(predicted_tensor, retrieved_tensor)
        
        # Compute difference: predicted_fp - overlap
        difference_fp = predicted_tensor - overlap
        difference_fp = torch.clamp(difference_fp, 0.0, 1.0)
        
        # Check if difference fingerprint is all zeros
        if difference_fp.sum() < 1e-6:
            logger.warning("Difference fingerprint is all zeros, returning empty results")
            return SecondaryRetrievalResponse(
                results=[],
                total_count=0,
                difference_fp=difference_fp.tolist()
            )
        
        # Use difference fingerprint directly (already in [0,1] range)
        # RankingSet will handle normalization internally
        difference_prob = difference_fp
        
        # Get rankingset
        rankingset = model_service.get_rankingset()
        if rankingset is None:
            raise HTTPException(
                status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
                detail="Rankingset not available"
            )
        
        # Use RankingSet to retrieve results
        ranker = RankingSet(store=rankingset, metric="cosine")
        
        # Retrieve top-k indices using difference fingerprint
        idxs = ranker.retrieve_idx(difference_prob.unsqueeze(0), n=k)
        
        # Get similarity scores for the retrieved indices
        sims = ranker._sims(difference_prob.unsqueeze(0))
        sims_sorted, _ = torch.topk(sims.squeeze(), k=k, dim=0)
        
        # Build results
        results = []
        metadata_service = MetadataService.instance()
        molecule_renderer = MoleculeRenderer.instance()
        fp_loader = model_service.get_fp_loader()
        
        for i, idx in enumerate(idxs.squeeze().tolist()):
            entry = metadata_service.get_entry(idx)
            if not entry:
                continue
            
            result_smiles = metadata_service.get_smiles(idx)
            if not result_smiles:
                continue
            
            # Render molecule with difference fingerprint highlighting
            try:
                enhanced_svg = molecule_renderer.render(
                    result_smiles,
                    predicted_fp=difference_prob,
                    fp_loader=fp_loader,
                    img_size=MOLECULE_IMG_SIZE
                )
            except Exception as e:
                logger.warning(f"Failed to render molecule: {e}")
                enhanced_svg = None
            
            # Similarity
            similarity_val = sims_sorted[i] if i < len(sims_sorted) and sims_sorted[i] is not None else 0.0
            try:
                similarity_float = float(similarity_val)
                if not isinstance(similarity_float, (int, float)) or similarity_float != similarity_float:
                    similarity_float = 0.0
            except (ValueError, TypeError):
                similarity_float = 0.0
            
            # Build card
            card = build_result_card(idx, entry, similarity_float, enhanced_svg)
            results.append(card)
        
        # Pagination
        offset = 0
        limit = len(results)
        paged_results = results[offset:offset + limit]
        
        return SecondaryRetrievalResponse(
            results=paged_results,
            total_count=len(results),
            difference_fp=difference_fp.tolist()
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Secondary retrieval error: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=str(e)
        )

