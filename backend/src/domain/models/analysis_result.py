from typing import Optional, Dict, Any, List
from pydantic import BaseModel, Field

from .spectral_data import SpectralDataInput
from .prediction_result import ResultCard as ResultCardModel


class AnalysisRequest(BaseModel):
    """Request model for /analyze endpoint"""
    smiles: str = Field(..., min_length=1, description="SMILES string of the molecule to analyze")
    model_id: Optional[str] = Field(None, description="Model to use (default from models.json)")
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
    model_id: Optional[str] = Field(None, description="Model to use (default from models.json)")


class SecondaryRetrievalResponse(BaseModel):
    """Response model for /secondary-retrieval endpoint"""
    results: list = Field(..., description="List of secondary retrieval results")
    total_count: int = Field(..., ge=0, description="Total number of results")
    difference_fp: Optional[List[float]] = Field(None, description="Difference fingerprint vector (predicted - overlap)")


class AblationRequest(BaseModel):
    """Request model for /ablation endpoint"""
    raw: SpectralDataInput = Field(..., description="Spectral data input for ablation analysis")
    smiles: str = Field(..., min_length=1, description="SMILES string of the molecule to compare against")
    model_id: Optional[str] = Field(None, description="Model to use (default from models.json)")
    bit_threshold: float = Field(0.5, ge=0.0, le=1.0, description="Threshold for considering fingerprint bits active")
    max_bits: int = Field(256, ge=1, le=2048, description="Maximum number of active bits to return")
    reference_fp: Optional[List[float]] = Field(None, description="Optional reference fingerprint vector to compare against")


class AblationResponse(BaseModel):
    """Response model for /ablation endpoint"""
    pred_fp: List[float] = Field(..., description="Predicted fingerprint vector from ablation run")
    active_bit_indices: List[int] = Field(..., description="Indices of fingerprint bits above the threshold")
    similarity_map: Optional[str] = Field(None, description="Base64-encoded similarity map image")
    bit_environments: Dict[int, Dict[str, Any]] = Field(default_factory=dict, description="Bit environments for active fingerprint bits")
    change_overlay_svg: Optional[str] = Field(None, description="SVG string showing gained and lost fingerprint bit overlays")


class CustomSmilesCardRequest(BaseModel):
    """Request model for /custom-smiles-card endpoint"""
    smiles: str = Field(..., min_length=1, description="Custom SMILES string to decorate")
    model_id: Optional[str] = Field(None, description="Model to use (default from models.json)")
    reference_fp: List[float] = Field(..., description="Reference fingerprint (e.g., predicted or query fingerprint)")


class CustomSmilesCardResponse(BaseModel):
    """Response model for /custom-smiles-card endpoint"""
    result: ResultCardModel = Field(..., description="Result card for the custom SMILES")

