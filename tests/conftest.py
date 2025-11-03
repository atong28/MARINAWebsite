import types
import pytest


@pytest.fixture(scope="session")
def app_client():
    # Ensure predictor model is stubbed before importing app
    import src.predictor as predictor
    predictor._model = object()
    predictor._fp_loader = None
    predictor._rankingset_cache = None

    import app as app_module
    # Mark initialized to bypass gating
    app_module._model_initialized = True
    app = app_module.app
    return app.test_client()


