"""
Pytest configuration and fixtures for MARINA backend tests.
"""
import pytest
from fastapi.testclient import TestClient

from src.api.app import app


@pytest.fixture
def client():
    """Create a test client for the FastAPI app."""
    return TestClient(app)


@pytest.fixture
def sample_spectral_data():
    """Sample spectral data for testing."""
    return {
        "raw": {
            "hsqc": [7.2, 120.5, 1.0, 6.8, 110.3, 0.8],
            "h_nmr": [7.2, 6.8, 3.5],
            "c_nmr": [120.5, 110.3, 55.2],
            "mass_spec": [180.5, 1000, 150.3, 800],
            "mw": 180.16
        },
        "k": 5
    }


@pytest.fixture
def sample_smiles():
    """Sample SMILES string for testing."""
    return "CCO"  # Ethanol

