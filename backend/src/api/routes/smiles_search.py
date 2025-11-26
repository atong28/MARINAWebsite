"""
SMILES search endpoint.
"""
import logging
import torch
from fastapi import APIRouter, HTTPException, status, Request

from src.domain.models.spectral_data import SmilesSearchRequest
from src.domain.models.prediction_result import SmilesSearchResponse
from src.services.model_service import ModelService
from src.services.metadata_service import MetadataService
from src.services.molecule_renderer import MoleculeRenderer
from src.services.result_builder import build_result_card
from src.domain.ranker import RankingSet
from src.domain.drawing.draw import compute_bit_environments_batch
from src.domain.fingerprint.fp_loader import EntropyFPLoader
from src.config import MOLECULE_IMG_SIZE, MAX_TOP_K
from src.api.middleware.rate_limit import get_limiter

logger = logging.getLogger(__name__)
router = APIRouter()

# Get limiter instance
limiter = get_limiter()


@router.post("/smiles-search", response_model=SmilesSearchResponse, status_code=status.HTTP_200_OK)
@limiter.limit("20 per minute")
async def smiles_search(request: Request, data: SmilesSearchRequest):
    """Search for similar molecules using a SMILES string."""
    # Check if model is ready
    model_service = ModelService.instance()
    if not model_service.is_ready():
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Model is not ready. Please try again in a moment."
        )
    
    # Validate k
    k = data.k
    if k > MAX_TOP_K:
        k = MAX_TOP_K
    
    try:
        smiles = data.smiles.strip()
        if not smiles:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="SMILES string cannot be empty"
            )
        
        logger.info(f"SMILES search for: {smiles} with k={k}")
        
        # Get fingerprint for SMILES
        fp_loader = model_service.get_fp_loader()
        if fp_loader is None:
            raise HTTPException(
                status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
                detail="Model not loaded. Please ensure the model is initialized first."
            )
        
        fp = fp_loader.build_mfp_for_smiles(smiles)
        if fp is None:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Invalid SMILES string or unable to generate fingerprint"
            )
        
        # Convert to tensor
        pred = torch.tensor(fp, dtype=torch.float)
        
        # Get rankingset
        rankingset = model_service.get_rankingset()
        if rankingset is None:
            raise HTTPException(
                status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
                detail="Rankingset not available"
            )
        
        # Use RankingSet to get results
        ranker = RankingSet(store=rankingset, metric="cosine")
        
        # Get similarity scores and indices together
        sims, idxs = ranker.retrieve_with_scores(pred.unsqueeze(0), n=k)
        sims_sorted, idxs_sorted = torch.topk(sims.squeeze(), k=k, dim=0)
        
        # Build results
        results = []
        metadata_service = MetadataService.instance()
        molecule_renderer = MoleculeRenderer.instance()
        
        for i, idx in enumerate(idxs.squeeze().tolist()):
            entry = metadata_service.get_entry(idx)
            if not entry:
                continue
            
            result_smiles = metadata_service.get_smiles(idx)
            if not result_smiles:
                continue
            
            # Render molecule
            try:
                query_tensor = torch.tensor(fp, dtype=torch.float32)
                enhanced_svg = molecule_renderer.render(
                    result_smiles, 
                    predicted_fp=query_tensor, 
                    fp_loader=fp_loader, 
                    img_size=MOLECULE_IMG_SIZE
                )
            except Exception as e:
                logger.warning(f"Failed to render molecule: {e}")
                enhanced_svg = None
            
            # Get plain SVG
            try:
                plain_svg = molecule_renderer.render(
                    result_smiles, 
                    predicted_fp=None, 
                    fp_loader=None, 
                    img_size=MOLECULE_IMG_SIZE
                )
            except Exception:
                plain_svg = None
            
            # Similarity (clamped to [0.0, 1.0])
            similarity_val = sims_sorted[i] if i < len(sims_sorted) and sims_sorted[i] is not None else 0.0
            try:
                similarity_float = float(similarity_val)
                if not isinstance(similarity_float, (int, float)) or similarity_float != similarity_float:
                    similarity_float = 0.0
            except (ValueError, TypeError):
                similarity_float = 0.0

            # Ensure similarity is always in [0.0, 1.0] before returning
            similarity_float = max(0.0, min(1.0, similarity_float))
            
            # Build card
            card = build_result_card(idx, entry, similarity_float, enhanced_svg)
            if plain_svg:
                card['plain_svg'] = plain_svg
            
            # Get retrieved molecule fingerprint indices
            retrieved_indices = []
            try:
                if rankingset is not None:
                    if rankingset.layout == torch.sparse_csr:
                        row_start = rankingset.crow_indices()[idx]
                        row_end = rankingset.crow_indices()[idx + 1]
                        if row_start < row_end:
                            col_indices = rankingset.col_indices()[row_start:row_end]
                            retrieved_indices = col_indices.cpu().tolist()
                    else:
                        row_tensor = rankingset[idx].cpu()
                        retrieved_indices = torch.nonzero(row_tensor > 0.5, as_tuple=False).squeeze(-1).tolist()
            except Exception as e:
                logger.warning(f"Failed to get retrieved indices: {e}")
            
            card['retrieved_molecule_fp_indices'] = retrieved_indices
            
            # Compute bit environments
            try:
                if retrieved_indices and fp_loader:
                    bit_environments = compute_bit_environments_batch(result_smiles, retrieved_indices, fp_loader)
                    card['bit_environments'] = bit_environments
                else:
                    card['bit_environments'] = {}
            except Exception as e:
                logger.warning(f"Failed to compute bit environments: {e}")
                card['bit_environments'] = {}
            
            results.append(card)
        
        # Pagination
        offset = 0
        limit = len(results)
        paged_results = results[offset:offset + limit]
        
        return SmilesSearchResponse(
            results=paged_results,
            total_count=len(results),
            offset=offset,
            limit=limit,
            query_smiles=smiles
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"SMILES search error: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=str(e)
        )

