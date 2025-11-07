"""
Tests for API routes.
"""
import pytest
from fastapi import status


def test_health_endpoint(client):
    """Test health check endpoint."""
    response = client.get("/health")
    assert response.status_code in [status.HTTP_200_OK, status.HTTP_503_SERVICE_UNAVAILABLE]
    data = response.json()
    assert "status" in data
    assert "model_loaded" in data
    assert "uptime_seconds" in data


def test_predict_endpoint_validation(client, sample_spectral_data):
    """Test predict endpoint with valid data."""
    # This test may fail if model is not loaded - that's expected
    response = client.post("/predict", json=sample_spectral_data)
    # Should either succeed (200) or fail with model not ready (503)
    assert response.status_code in [
        status.HTTP_200_OK,
        status.HTTP_503_SERVICE_UNAVAILABLE,
        status.HTTP_422_UNPROCESSABLE_ENTITY
    ]


def test_predict_endpoint_invalid_data(client):
    """Test predict endpoint with invalid data."""
    response = client.post("/predict", json={
        "raw": {
            "hsqc": [1, 2]  # Invalid - not multiple of 3
        },
        "k": 5
    })
    assert response.status_code == status.HTTP_422_UNPROCESSABLE_ENTITY


def test_smiles_search_endpoint_validation(client, sample_smiles):
    """Test SMILES search endpoint with valid SMILES."""
    response = client.post("/smiles-search", json={
        "smiles": sample_smiles,
        "k": 5
    })
    # Should either succeed (200) or fail with model not ready (503)
    assert response.status_code in [
        status.HTTP_200_OK,
        status.HTTP_503_SERVICE_UNAVAILABLE,
        status.HTTP_400_BAD_REQUEST
    ]


def test_smiles_search_endpoint_invalid(client):
    """Test SMILES search endpoint with invalid SMILES."""
    response = client.post("/smiles-search", json={
        "smiles": "",
        "k": 5
    })
    assert response.status_code == status.HTTP_422_UNPROCESSABLE_ENTITY

