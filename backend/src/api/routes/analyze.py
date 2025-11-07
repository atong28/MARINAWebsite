"""
Analysis endpoint for detailed molecular fingerprint analysis.
"""
import logging
import torch
from fastapi import APIRouter, HTTPException, status, Request

from src.domain.models.analysis_result import AnalysisRequest, AnalysisResponse
from src.services.model_service import ModelService
from src.services.molecule_renderer import MoleculeRenderer
from src.domain.drawing.draw import compute_bit_environments_batch
from src.config import MOLECULE_IMG_SIZE
from src.api.middleware.rate_limit import get_limiter

logger = logging.getLogger(__name__)
router = APIRouter()

# Get limiter instance
limiter = get_limiter()


@router.post("/analyze", response_model=AnalysisResponse, status_code=status.HTTP_200_OK)
@limiter.limit("10 per minute")
async def analyze(request: Request, data: AnalysisRequest):
    """Analysis endpoint for detailed molecular fingerprint analysis."""
    # Guard clause: Check if model is ready
    model_service = ModelService.instance()
    if not model_service.is_ready():
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Model is not ready. Please try again in a moment."
        )
    
    # Guard clause: Validate SMILES
    target_smiles = data.smiles.strip()
    if not target_smiles:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="SMILES string cannot be empty"
        )
    
    # Guard clause: Validate retrieved_fp is provided
    retrieved_fp = data.retrieved_fp
    
    if retrieved_fp is None:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="retrieved_fp is required"
        )
    
    if len(retrieved_fp) == 0:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="retrieved_fp array cannot be empty"
        )
    
    # Guard clause: Get fp_loader and validate it exists
    fp_loader = model_service.get_fp_loader()
    if fp_loader is None:
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Fingerprint loader is not available. Model may not be fully initialized."
        )
    
    try:
        retrieved_tensor = torch.tensor(retrieved_fp, dtype=torch.float32)
        retrieved_molecule_fp_indices = []
        try:
            retrieved_molecule_fp_indices = torch.nonzero(retrieved_tensor > 0.5, as_tuple=False).squeeze(-1).tolist()
            if not isinstance(retrieved_molecule_fp_indices, list):
                retrieved_molecule_fp_indices = [retrieved_molecule_fp_indices] if retrieved_molecule_fp_indices is not None else []
        except Exception as e:
            logger.error(f"[Analyze] Failed to extract retrieved fingerprint indices: {e}", exc_info=True)
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail=f"Failed to extract fingerprint indices: {str(e)}"
            )
        
        if not retrieved_molecule_fp_indices:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="No active fingerprint bits found in retrieved_fp"
            )
        
        try:
            bit_environments = compute_bit_environments_batch(
                smiles=target_smiles,
            fp_indices=retrieved_molecule_fp_indices,
                fp_loader=fp_loader
            )
        except Exception as e:
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail=f"Failed to compute bit environments: {str(e)}"
                )
        molecule_renderer = MoleculeRenderer.instance()
        molecule_svg = molecule_renderer.render(
            smiles=target_smiles,
            img_size=MOLECULE_IMG_SIZE,
            bit_environments=bit_environments
        )
        
        if not molecule_svg:
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail="Molecule rendering failed: renderer returned None"
            )
        
        return AnalysisResponse(
            retrieved_molecule_fp_indices=retrieved_molecule_fp_indices,
            molecule_svg=molecule_svg
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"[Analyze] Analysis error: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Analysis failed: {str(e)}"
        )
