import json
import types


def test_health_ok(app_client):
    r = app_client.get('/health')
    assert r.status_code in (200, 202, 503)


def test_openapi(app_client):
    r = app_client.get('/openapi.json')
    assert r.status_code == 200
    data = r.get_json()
    assert data.get('openapi')


def test_predict_requires_raw_and_paginates(app_client, monkeypatch):
    # Missing body
    r = app_client.post('/predict', json={})
    assert r.status_code == 400

    # Stub predictor output and metadata
    import src.predictor as predictor
    def fake_predict_from_raw(raw, k=10):
        return ([0.9, 0.8], [1, 2], [0.0, 1.0, 0.0])
    monkeypatch.setattr(predictor, 'predict_from_raw', fake_predict_from_raw)

    from src.services.metadata_service import MetadataService
    svc = MetadataService.instance()
    monkeypatch.setattr(svc, 'get_entry', lambda idx: {'canonical_3d_smiles': 'N/A', 'canonical_2d_smiles': 'C'})
    monkeypatch.setattr(svc, 'get_smiles', lambda idx: 'C')

    # Paginate
    r = app_client.post('/predict', json={'raw': {}, 'k': 2, 'offset': 1, 'limit': 1})
    assert r.status_code == 200
    data = r.get_json()
    assert data['success'] is True
    assert len(data['results']) == 1
    assert data['offset'] == 1
    assert data['limit'] == 1
    assert data['total_count'] == 2


def test_smiles_search_requires_smiles(app_client):
    r = app_client.post('/smiles-search', json={})
    assert r.status_code == 400

