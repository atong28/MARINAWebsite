"""
Analysis and custom-card endpoints for detailed molecular fingerprint analysis.
"""
import logging
from typing import List

import torch
from fastapi import APIRouter, HTTPException, status, Request

from src.domain.models.analysis_result import (
    AnalysisRequest,
    AnalysisResponse,
    CustomSmilesCardRequest,
    CustomSmilesCardResponse,
)
from src.services.model_manifest import resolve_and_validate_model_id
from src.services.model_service import ModelService
from src.services.molecule_renderer import MoleculeRenderer
from src.services.result_builder import build_result_card
from src.domain.fingerprint.fp_loader import EntropyFPLoader
from src.domain.fingerprint.fp_utils import tanimoto_similarity
from src.api.app import run_heavy
from src.config import ANALYZE_TIMEOUT_S, CUSTOM_SMILES_TIMEOUT_S, MOLECULE_IMG_SIZE
from src.api.middleware.rate_limit import get_limiter
from src.services.compute_pool import ComputeOverloadedError, ComputeTimeoutError

logger = logging.getLogger(__name__)
router = APIRouter()

# Get limiter instance
limiter = get_limiter()


@router.post("/analyze", response_model=AnalysisResponse, status_code=status.HTTP_200_OK)
@limiter.limit("10 per minute")
async def analyze(request: Request, data: AnalysisRequest):
    """Analysis endpoint for detailed molecular fingerprint analysis."""
    mid, err = resolve_and_validate_model_id(data.model_id)
    if err is not None:
        code, detail = err
        raise HTTPException(status_code=code, detail=detail)

    model_service = ModelService.instance()
    model_service.ensure_loaded(mid)
    
    target_smiles = data.smiles.strip()
    if not target_smiles:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="SMILES string cannot be empty"
        )
    
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
    
    fp_loader = model_service.get_fp_loader(mid)
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
            pool = request.app.state.compute_pool
            bit_environments = await run_heavy(
                request,
                pool.run(
                    "bit_envs",
                    {"smiles": target_smiles, "fp_indices": retrieved_molecule_fp_indices},
                    timeout=ANALYZE_TIMEOUT_S,
                ),
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
    except ComputeOverloadedError:
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Server is busy. Please retry shortly.",
        )
    except ComputeTimeoutError:
        raise HTTPException(
            status_code=status.HTTP_504_GATEWAY_TIMEOUT,
            detail="Analysis timed out.",
        )
    except Exception as e:
        logger.error(f"[Analyze] Analysis error: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Analysis failed: {str(e)}"
        )


@router.post("/custom-smiles-card", response_model=CustomSmilesCardResponse, status_code=status.HTTP_200_OK)
@limiter.limit("10 per minute")
async def custom_smiles_card(request: Request, data: CustomSmilesCardRequest):
    """
    Build a single result card for an arbitrary SMILES using a provided reference
    fingerprint (e.g., predicted or query FP), without running a new retrieval.
    """
    mid, err = resolve_and_validate_model_id(data.model_id)
    if err is not None:
        code, detail = err
        raise HTTPException(status_code=code, detail=detail)

    model_service = ModelService.instance()
    model_service.ensure_loaded(mid)

    smiles = data.smiles.strip()
    if not smiles:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="SMILES string cannot be empty",
        )

    reference_fp: List[float] = data.reference_fp
    if not reference_fp:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="reference_fp cannot be empty",
        )

    fp_loader = model_service.get_fp_loader(mid)
    if fp_loader is None or not isinstance(fp_loader, EntropyFPLoader):
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Fingerprint loader is not available. Model may not be fully initialized.",
        )

    try:
        ref_tensor = torch.tensor(reference_fp, dtype=torch.float32)
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid reference_fp: {str(e)}",
        )

    try:
        # Build special retrieval fingerprint for the custom SMILES
        special_mfp = fp_loader.build_mfp_for_smiles(smiles)
    except Exception as e:
        logger.error(f"[custom-smiles-card] Failed to build fingerprint for SMILES: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Invalid SMILES string or unable to generate fingerprint",
        )

    # Ensure vectors have matching length
    if ref_tensor.numel() != special_mfp.numel():
        n = min(ref_tensor.numel(), special_mfp.numel())
        ref_tensor = ref_tensor.view(-1)[:n]
        special_mfp = special_mfp.view(-1)[:n]
    else:
        ref_tensor = ref_tensor.view(-1)
        special_mfp = special_mfp.view(-1)

    # Cosine similarity in [0, 1]
    try:
        dot = float(torch.dot(ref_tensor, special_mfp))
        ref_norm = float(torch.linalg.norm(ref_tensor) + 1e-8)
        special_norm = float(torch.linalg.norm(special_mfp) + 1e-8)
        cos_raw = dot / (ref_norm * special_norm) if ref_norm > 0 and special_norm > 0 else 0.0
    except Exception as e:
        logger.warning(f"[custom-smiles-card] Failed to compute cosine similarity: {e}")
        cos_raw = 0.0
    cosine_sim = max(0.0, min(1.0, float(cos_raw))) if cos_raw == cos_raw else 0.0

    # Tanimoto similarity
    try:
        tanimoto_sim = tanimoto_similarity(ref_tensor, special_mfp)
    except Exception as e:
        logger.warning(f"[custom-smiles-card] Failed to compute Tanimoto similarity: {e}")
        tanimoto_sim = 0.0

    # Render molecule SVGs
    molecule_renderer = MoleculeRenderer.instance()
    try:
        enhanced_svg = molecule_renderer.render(
            smiles=smiles,
            predicted_fp=ref_tensor,
            fp_loader=fp_loader,
            img_size=MOLECULE_IMG_SIZE,
        )
    except Exception as e:
        logger.warning(f"[custom-smiles-card] Failed to render enhanced molecule: {e}")
        enhanced_svg = None

    try:
        plain_svg = molecule_renderer.render(
            smiles=smiles,
            predicted_fp=None,
            fp_loader=None,
            img_size=MOLECULE_IMG_SIZE,
        )
    except Exception:
        plain_svg = None

    # Build a minimal metadata entry using the provided SMILES
    entry = {
        "canonical_3d_smiles": smiles,
        "canonical_2d_smiles": smiles,
    }

    # Build base card
    card = build_result_card(
        idx=-1,
        entry=entry,
        similarity=cosine_sim,
        svg=enhanced_svg,
        cosine_similarity=cosine_sim,
        tanimoto_similarity=tanimoto_sim,
    )
    if plain_svg:
        card["plain_svg"] = plain_svg

    # Retrieved fingerprint indices and bit environments for the custom SMILES
    try:
        retrieved_tensor = special_mfp
        retrieved_indices = torch.nonzero(retrieved_tensor > 0.5, as_tuple=False).squeeze(-1).tolist()
        if not isinstance(retrieved_indices, list):
            retrieved_indices = [retrieved_indices] if retrieved_indices is not None else []
    except Exception as e:
        logger.warning(f"[custom-smiles-card] Failed to extract fingerprint indices: {e}")
        retrieved_indices = []

    card["retrieved_molecule_fp_indices"] = retrieved_indices

    try:
        if retrieved_indices:
            pool = request.app.state.compute_pool
            bit_envs = await run_heavy(
                request,
                pool.run(
                    "bit_envs",
                    {"smiles": smiles, "fp_indices": retrieved_indices},
                    timeout=CUSTOM_SMILES_TIMEOUT_S,
                ),
            )
            card["bit_environments"] = bit_envs
        else:
            card["bit_environments"] = {}
    except Exception as e:
        logger.warning(f"[custom-smiles-card] Failed to compute bit environments: {e}")
        card["bit_environments"] = {}

    return CustomSmilesCardResponse(result=card)
