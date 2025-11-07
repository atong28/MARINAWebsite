from typing import Optional, Dict, Any, List
from pydantic import BaseModel, Field


class AnalysisRequest(BaseModel):
    """Request model for /analyze endpoint"""
    smiles: str = Field(..., min_length=1, description="SMILES string of the molecule to analyze")
    predicted_fp: Optional[List[float]] = Field(None, description="Predicted fingerprint vector")
    retrieved_fp: Optional[List[float]] = Field(None, description="Retrieved molecule fingerprint vector")


class AnalysisResponse(BaseModel):
    """Response model for /analyze endpoint"""
    retrieved_molecule_fp_indices: List[int] = Field(..., description="Indices of retrieved molecule fingerprint bits")
    molecule_svg: str = Field(..., description="Molecule SVG with embedded bit environment overlays")


class SecondaryRetrievalRequest(BaseModel):
    """Request model for /secondary-retrieval endpoint"""
    predicted_fp: List[float] = Field(..., description="Predicted fingerprint vector")
    retrieved_fp: List[float] = Field(..., description="Retrieved molecule fingerprint vector")
    k: int = Field(default=10, ge=1, le=50, description="Number of results to retrieve")


class SecondaryRetrievalResponse(BaseModel):
    """Response model for /secondary-retrieval endpoint"""
    results: list = Field(..., description="List of secondary retrieval results")
    total_count: int = Field(..., ge=0, description="Total number of results")
    difference_fp: Optional[List[float]] = Field(None, description="Difference fingerprint vector (predicted - overlap)")

