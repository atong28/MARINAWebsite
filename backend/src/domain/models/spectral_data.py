from typing import List, Optional
from pydantic import BaseModel, Field, field_validator

class HSQCData(BaseModel):
    """HSQC NMR data: triplets of [H_shift, C_shift, intensity]"""
    hsqc: List[float] = Field(..., description="Array of HSQC triplets: [H1, C1, I1, H2, C2, I2, ...]")
    
    @field_validator('hsqc')
    def validate_hsqc_length(cls, v):
        if len(v) % 3 != 0:
            raise ValueError('HSQC data must be triplets (length must be multiple of 3)')
        return v


class HNMRData(BaseModel):
    """H NMR data: array of chemical shift values"""
    h_nmr: List[float] = Field(..., description="Array of H NMR chemical shifts in ppm")


class CNMRData(BaseModel):
    """C NMR data: array of chemical shift values"""
    c_nmr: List[float] = Field(..., description="Array of C NMR chemical shifts in ppm")


class MassSpecData(BaseModel):
    """Mass spectrometry data: pairs of [m/z, intensity]"""
    mass_spec: List[float] = Field(..., description="Array of mass spec pairs: [mz1, I1, mz2, I2, ...]")
    
    @field_validator('mass_spec')
    def validate_mass_spec_length(cls, v):
        if len(v) % 2 != 0:
            raise ValueError('Mass spec data must be pairs (length must be multiple of 2)')
        return v


class MolecularWeightData(BaseModel):
    """Molecular weight data: single scalar value"""
    mw: float = Field(..., gt=0, description="Molecular weight in g/mol")


class SpectralDataInput(BaseModel):
    """Combined spectral data input for prediction"""
    hsqc: Optional[List[float]] = None
    h_nmr: Optional[List[float]] = None
    c_nmr: Optional[List[float]] = None
    mass_spec: Optional[List[float]] = None
    mw: Optional[float] = None
    
    @field_validator('hsqc')
    def validate_hsqc(cls, v):
        if v is not None and len(v) % 3 != 0:
            raise ValueError('HSQC data must be triplets (length must be multiple of 3)')
        return v
    
    @field_validator('mass_spec')
    def validate_mass_spec(cls, v):
        if v is not None and len(v) % 2 != 0:
            raise ValueError('Mass spec data must be pairs (length must be multiple of 2)')
        return v
    
    @field_validator('mw')
    def validate_mw(cls, v):
        if v is not None and v <= 0:
            raise ValueError('Molecular weight must be positive')
        return v


class PredictRequest(BaseModel):
    """Request model for /predict endpoint"""
    raw: SpectralDataInput = Field(..., description="Spectral data input")
    k: int = Field(default=10, ge=1, le=50, description="Number of results to retrieve")
    mw_min: Optional[float] = Field(
        None,
        gt=0,
        description="Optional minimum molecular weight filter for retrieval (g/mol)",
    )
    mw_max: Optional[float] = Field(
        None,
        gt=0,
        description="Optional maximum molecular weight filter for retrieval (g/mol)",
    )


class SmilesSearchRequest(BaseModel):
    """Request model for /smiles-search endpoint"""
    smiles: str = Field(..., min_length=1, description="SMILES string to search")
    k: int = Field(default=10, ge=1, le=50, description="Number of results to retrieve")
    mw_min: Optional[float] = Field(
        None,
        gt=0,
        description="Optional minimum molecular weight filter for retrieval (g/mol)",
    )
    mw_max: Optional[float] = Field(
        None,
        gt=0,
        description="Optional maximum molecular weight filter for retrieval (g/mol)",
    )

