from typing import List, Optional, Dict, Any
from pydantic import BaseModel, Field


class DatabaseLinks(BaseModel):
    """Links to external databases"""
    coconut: Optional[str] = None
    lotus: Optional[str] = None
    npmrd: Optional[str] = None


class ResultCard(BaseModel):
    """Single result card in prediction response"""
    index: int = Field(..., description="Database index of the molecule")
    smiles: str = Field(..., description="SMILES string of the molecule")
    similarity: float = Field(..., ge=0.0, le=1.0, description="Tanimoto similarity score")
    svg: Optional[str] = Field(None, description="Enhanced SVG representation with highlighting")
    plain_svg: Optional[str] = Field(None, description="Plain SVG representation without highlighting")
    name: Optional[str] = Field(None, description="Primary name of the molecule")
    primary_link: Optional[str] = Field(None, description="Primary database link")
    database_links: DatabaseLinks = Field(default_factory=DatabaseLinks, description="All database links")
    retrieved_molecule_fp_indices: List[int] = Field(default_factory=list, description="Fingerprint bit indices for this molecule")
    bit_environments: Dict[int, Dict[str, Any]] = Field(default_factory=dict, description="Bit environment data for each bit index")
    exact_mass: Optional[float] = Field(None, description="Exact molecular mass in g/mol computed with RDKit")


class PredictResponse(BaseModel):
    """Response model for /predict endpoint"""
    results: List[ResultCard] = Field(..., description="List of result cards")
    total_count: int = Field(..., ge=0, description="Total number of results")
    offset: int = Field(default=0, ge=0, description="Pagination offset")
    limit: int = Field(..., ge=0, description="Pagination limit")
    pred_fp: Optional[List[float]] = Field(None, description="Predicted fingerprint vector")


class SmilesSearchResponse(BaseModel):
    """Response model for /smiles-search endpoint"""
    results: List[ResultCard] = Field(..., description="List of result cards")
    total_count: int = Field(..., ge=0, description="Total number of results")
    offset: int = Field(default=0, ge=0, description="Pagination offset")
    limit: int = Field(..., ge=0, description="Pagination limit")
    query_smiles: str = Field(..., description="The SMILES string that was searched")

