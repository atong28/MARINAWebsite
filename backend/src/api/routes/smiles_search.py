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
from src.domain.ranker import RankingSet, filter_rankingset_rows
from src.domain.drawing.draw import compute_bit_environments_batch
from src.domain.fingerprint.fp_loader import EntropyFPLoader
from src.domain.fingerprint.fp_utils import tanimoto_similarity
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
    requested_k = data.k
    if requested_k > MAX_TOP_K:
        requested_k = MAX_TOP_K
    
    try:
        smiles = data.smiles.strip()
        if not smiles:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="SMILES string cannot be empty"
            )
        
        logger.info(f"SMILES search for: {smiles} with k={requested_k}")
        
        # Validate MW filter bounds (one-sided filters allowed)
        mw_min = data.mw_min
        mw_max = data.mw_max
        if mw_min is not None and mw_max is not None and mw_min > mw_max:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="mw_min cannot be greater than mw_max",
            )

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
        pred = torch.tensor(fp, dtype=torch.float32)
        
        # Get full rankingset and prefilter rows by MW range
        rankingset_full = model_service.get_rankingset()
        if rankingset_full is None:
            raise HTTPException(
                status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
                detail="Rankingset not available"
            )
        
        num_rows = rankingset_full.shape[0]
        kept_indices: list[int] = []
        for idx in range(num_rows):
            if metadata_service.within_mw_range(idx, mw_min=mw_min, mw_max=mw_max):
                kept_indices.append(idx)

        query_tensor = torch.tensor(fp, dtype=torch.float32)

        if not kept_indices:
            # No molecules satisfy MW; return empty result set
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
        
        # Build results
        results = []
        metadata_service = MetadataService.instance()
        molecule_renderer = MoleculeRenderer.instance()
        
        def _dense_from_indices(length: int, indices: list[int]) -> torch.Tensor:
            vec = torch.zeros(length, dtype=torch.float32)
            if not indices:
                return vec
            for j in indices:
                if 0 <= j < length:
                    vec[j] = 1.0
            return vec

        for i, local_row in enumerate(local_idxs.tolist()):
            global_idx = kept_indices[local_row]
            entry = metadata_service.get_entry(global_idx)
            if not entry:
                continue
            
            result_smiles = metadata_service.get_smiles(global_idx)
            if not result_smiles:
                continue

            # Render molecule
            try:
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
            similarity_val = sims[i] if i < len(sims) and sims[i] is not None else 0.0
            try:
                cosine_sim = float(similarity_val)
                if not isinstance(cosine_sim, (int, float)) or cosine_sim != cosine_sim:
                    cosine_sim = 0.0
            except (ValueError, TypeError):
                cosine_sim = 0.0

            # Ensure similarity is always in [0.0, 1.0] before returning
            cosine_sim = max(0.0, min(1.0, cosine_sim))

            # Get retrieved molecule fingerprint indices
            retrieved_indices = []
            try:
                rankingset = rankingset_full
                if rankingset is not None:
                    if rankingset.layout == torch.sparse_csr:
                        row_start = rankingset.crow_indices()[global_idx]
                        row_end = rankingset.crow_indices()[global_idx + 1]
                        if row_start < row_end:
                            col_indices = rankingset.col_indices()[row_start:row_end]
                            retrieved_indices = col_indices.cpu().tolist()
                    else:
                        row_tensor = rankingset[global_idx].cpu()
                        retrieved_indices = torch.nonzero(row_tensor > 0.5, as_tuple=False).squeeze(-1).tolist()
            except Exception as e:
                logger.warning(f"Failed to get retrieved indices: {e}")

            # Compute Tanimoto similarity between query fingerprint and retrieved vector
            tanimoto_sim = 0.0
            try:
                if retrieved_indices:
                    retrieved_vec = _dense_from_indices(query_tensor.numel(), retrieved_indices)
                    tanimoto_sim = tanimoto_similarity(query_tensor, retrieved_vec)
            except Exception as e:
                logger.warning(f"Failed to compute Tanimoto similarity for idx {idx}: {e}")

            # Build card with both cosine and Tanimoto similarities
            card = build_result_card(
                idx,
                entry,
                cosine_sim,
                enhanced_svg,
                cosine_similarity=cosine_sim,
                tanimoto_similarity=tanimoto_sim,
            )
            if plain_svg:
                card['plain_svg'] = plain_svg

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

            # Stop once we've collected the requested number of MW-filtered results
            if len(results) >= requested_k:
                break
        
        # Pagination
        offset = 0
        limit = len(results)
        paged_results = results[offset:offset + limit]
        
        return SmilesSearchResponse(
            results=paged_results,
            total_count=len(results),
            offset=offset,
            limit=limit,
            query_smiles=smiles,
            query_fp=query_tensor.tolist(),
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"SMILES search error: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=str(e)
        )

