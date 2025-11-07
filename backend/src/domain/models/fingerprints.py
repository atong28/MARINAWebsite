"""Pydantic models for fingerprint endpoints."""

from typing import Dict, List, Optional

from pydantic import BaseModel, Field


class FingerprintIndicesRequest(BaseModel):
    smiles: str = Field(..., min_length=1, description="SMILES string")


class FingerprintIndicesResponse(BaseModel):
    smiles: str = Field(..., description="Input SMILES string")
    fp_indices: Optional[List[int]] = Field(
        None,
        description="Sorted list of active fingerprint indices. None if SMILES could not be processed.",
    )


class FingerprintBatchRequest(BaseModel):
    smiles_list: List[str] = Field(
        ..., min_length=1, description="List of SMILES strings to process"
    )


class FingerprintBatchResponse(BaseModel):
    fingerprints: Dict[str, Optional[List[int]]] = Field(
        ..., description="Mapping of SMILES to sorted fingerprint indices (or None on failure)"
    )

