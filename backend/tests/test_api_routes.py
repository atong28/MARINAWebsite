"""
Tests for API routes.
"""
import pytest
from fastapi import status


def test_health_endpoint(client):
    """Test health check endpoint."""
    response = client.get("/api/health")
    assert response.status_code in [status.HTTP_200_OK, status.HTTP_503_SERVICE_UNAVAILABLE]
    data = response.json()
    assert "status" in data
    assert "model_loaded" in data
    assert "uptime_seconds" in data


def test_predict_endpoint_validation(client, sample_spectral_data):
    """Test predict endpoint with valid data."""
    # 200 if model loaded, 503 if not (e.g. no data/marina_best)
    response = client.post("/api/predict", json=sample_spectral_data)
    assert response.status_code in [
        status.HTTP_200_OK,
        status.HTTP_503_SERVICE_UNAVAILABLE,
        status.HTTP_422_UNPROCESSABLE_ENTITY
    ]


def test_predict_endpoint_invalid_data(client):
    """Test predict endpoint with invalid data."""
    response = client.post("/api/predict", json={
        "raw": {
            "hsqc": [1, 2]  # Invalid - not multiple of 3
        },
        "k": 5
    })
    assert response.status_code == status.HTTP_422_UNPROCESSABLE_ENTITY


def test_smiles_search_endpoint_validation(client, sample_smiles):
    """Test SMILES search endpoint with valid SMILES."""
    response = client.post("/api/smiles-search", json={
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
    response = client.post("/api/smiles-search", json={
        "smiles": "",
        "k": 5
    })
    assert response.status_code == status.HTTP_422_UNPROCESSABLE_ENTITY


def test_models_endpoint(client):
    """GET /api/models returns model list and default_model_id."""
    response = client.get("/api/models")
    assert response.status_code == status.HTTP_200_OK
    data = response.json()
    assert "models" in data
    assert "default_model_id" in data
    assert isinstance(data["models"], list)


def test_predict_model_id_spectre_501(client, sample_spectral_data):
    """Requesting model_id=spectre_best returns 501 (not yet supported)."""
    payload = {**sample_spectral_data, "model_id": "spectre_best"}
    response = client.post("/api/predict", json=payload)
    assert response.status_code == 501
    assert "spectre" in response.json().get("detail", "").lower()


def test_predict_model_id_unknown_400(client, sample_spectral_data):
    """Requesting unknown model_id returns 400."""
    payload = {**sample_spectral_data, "model_id": "nonexistent_model"}
    response = client.post("/api/predict", json=payload)
    assert response.status_code == status.HTTP_400_BAD_REQUEST
    assert "unknown" in response.json().get("detail", "").lower()

