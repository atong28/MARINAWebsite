"""Fingerprint utility endpoints."""

import logging

from fastapi import APIRouter, HTTPException, status

from src.domain.models.fingerprints import (
    FingerprintBatchRequest,
    FingerprintBatchResponse,
    FingerprintIndicesRequest,
    FingerprintIndicesResponse,
)
from src.services.model_manifest import resolve_and_validate_model_id
from src.services.model_service import ModelService

logger = logging.getLogger(__name__)
router = APIRouter()


def _get_fp_loader(model_id: str | None):
    mid, err = resolve_and_validate_model_id(model_id)
    if err is not None:
        code, detail = err
        raise HTTPException(status_code=code, detail=detail)
    model_service = ModelService.instance()
    if not model_service.is_ready(mid):
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail=f"Model '{mid}' is not available (not loaded).",
        )
    fp_loader = model_service.get_fp_loader(mid)
    if fp_loader is None:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Fingerprint loader is not available.",
        )
    return fp_loader


@router.post(
    "/fingerprints/indices",
    response_model=FingerprintIndicesResponse,
    status_code=status.HTTP_200_OK,
)
async def fingerprint_indices(payload: FingerprintIndicesRequest) -> FingerprintIndicesResponse:
    """Return entropy fingerprint indices for a single SMILES string."""
    fp_loader = _get_fp_loader(payload.model_id)

    indices = fp_loader.build_fp_indices_for_smiles(payload.smiles)
    if indices is None:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail="Failed to compute fingerprint indices for the provided SMILES.",
        )

    return FingerprintIndicesResponse(smiles=payload.smiles, fp_indices=indices)


@router.post(
    "/fingerprints/batch-indices",
    response_model=FingerprintBatchResponse,
    status_code=status.HTTP_200_OK,
)
async def fingerprint_batch_indices(payload: FingerprintBatchRequest) -> FingerprintBatchResponse:
    """Return entropy fingerprint indices for a list of SMILES strings."""
    fp_loader = _get_fp_loader(payload.model_id)

    fingerprints = fp_loader.build_fp_dict_for_smiles(payload.smiles_list)

    return FingerprintBatchResponse(fingerprints=fingerprints)


