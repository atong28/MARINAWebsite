"""
Prediction endpoint for spectral data.
"""
import logging
import torch
from fastapi import APIRouter, HTTPException, status, Request

from src.domain.models.spectral_data import PredictRequest
from src.domain.models.prediction_result import PredictResponse
from src.services.model_service import ModelService
from src.services.metadata_service import MetadataService
from src.services.molecule_renderer import MoleculeRenderer
from src.services.result_builder import build_result_card
from src.domain.predictor import predict_from_raw
from src.domain.drawing.draw import compute_bit_environments_batch
from src.config import MOLECULE_IMG_SIZE, MAX_TOP_K
from src.api.middleware.rate_limit import get_limiter

logger = logging.getLogger(__name__)
router = APIRouter()

# Get limiter instance
limiter = get_limiter()


@router.post("/predict", response_model=PredictResponse, status_code=status.HTTP_200_OK)
@limiter.limit("30 per minute")
async def predict(request: Request, data: PredictRequest):
    """
    Predict molecular structures from spectral data.
    The request body is a JSON object with the following fields:
    - k: The number of predictions to return
    - raw: A JSON object with the following fields:
      - hsqc: A list of floats in a flattened 1d list: [delta H, delta C, intensity, delta H, delta C, intensity, ...]
      - h_nmr: A list of floats of H chemical shifts
      - c_nmr: A list of floats of C chemical shifts
      - mass_spec: A list of floats in a flattened 1d list: [m/z, intensity, m/z, intensity, ...]
      - mw: A float representing the molecular weight
    """
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
        # Convert request to raw input format
        raw_data = {}
        if data.raw.hsqc:
            raw_data['hsqc'] = data.raw.hsqc
        if data.raw.h_nmr:
            raw_data['h_nmr'] = data.raw.h_nmr
        if data.raw.c_nmr:
            raw_data['c_nmr'] = data.raw.c_nmr
        if data.raw.mass_spec:
            raw_data['mass_spec'] = data.raw.mass_spec
        if data.raw.mw:
            raw_data['mw'] = data.raw.mw
        
        # Run prediction
        out = predict_from_raw(raw_data, k=k)
        
        if not isinstance(out, tuple) or len(out) < 2:
            logger.error(f"Unexpected predictor output: {out}")
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail="Unexpected predictor output"
            )
        
        scores, indices = out[0], out[1]
        pred_fp = out[2] if len(out) > 2 else None
        
        # Build results
        results = []
        metadata_service = MetadataService.instance()
        molecule_renderer = MoleculeRenderer.instance()
        fp_loader = model_service.get_fp_loader()
        
        for i, idx in enumerate(indices):
            entry = metadata_service.get_entry(idx)
            if not entry:
                continue
            
            result_smiles = metadata_service.get_smiles(idx)
            if not result_smiles:
                continue
            
            # Render molecule
            try:
                pred_tensor = torch.tensor(pred_fp, dtype=torch.float32) if pred_fp is not None else None
                enhanced_svg = molecule_renderer.render(
                    result_smiles, 
                    predicted_fp=pred_tensor, 
                    fp_loader=fp_loader, 
                    img_size=MOLECULE_IMG_SIZE
                )
            except Exception as e:
                logger.warning(f"Failed to render enhanced molecule: {e}")
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
            
            # Similarity
            similarity_val = scores[i] if i < len(scores) and scores[i] is not None else 0.0
            try:
                similarity_float = float(similarity_val)
                if not isinstance(similarity_float, (int, float)) or similarity_float != similarity_float:
                    similarity_float = 0.0
            except (ValueError, TypeError):
                similarity_float = 0.0
            
            # Build card
            card = build_result_card(idx, entry, similarity_float, enhanced_svg)
            if plain_svg:
                card['plain_svg'] = plain_svg
            
            # Get retrieved molecule fingerprint indices
            retrieved_indices = []
            try:
                rankingset = model_service.get_rankingset()
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
        
        return PredictResponse(
            results=paged_results,
            total_count=len(results),
            offset=offset,
            limit=limit,
            pred_fp=pred_fp
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Prediction error: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=str(e)
        )

