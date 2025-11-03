import time
import threading
import queue
import pytest


@pytest.mark.load
def test_predict_concurrent(app_client, monkeypatch):
    # Stub predictor to be fast
    import src.predictor as predictor
    def fake_predict_from_raw(raw, k=10):
        return ([0.9] * k, list(range(1, k + 1)), [0.0, 1.0, 0.0])
    monkeypatch.setattr(predictor, 'predict_from_raw', fake_predict_from_raw)

    # Minimal metadata lookups
    from src.services.metadata_service import MetadataService
    svc = MetadataService.instance()
    monkeypatch.setattr(svc, 'get_entry', lambda idx: {'canonical_3d_smiles': 'N/A', 'canonical_2d_smiles': 'C'})
    monkeypatch.setattr(svc, 'get_smiles', lambda idx: 'C')

    N = 30
    errors = queue.Queue()

    def worker():
        try:
            r = app_client.post('/predict', json={'raw': {}, 'k': 5})
            if r.status_code != 200:
                errors.put(r.status_code)
            else:
                data = r.get_json()
                if not data.get('success'):
                    errors.put('no-success')
        except Exception as e:
            errors.put(str(e))

    threads = [threading.Thread(target=worker) for _ in range(N)]
    start = time.time()
    for t in threads: t.start()
    for t in threads: t.join()
    elapsed = time.time() - start

    # No errors and time reasonable for stubbed responses
    assert errors.qsize() == 0
    assert elapsed < 5.0


