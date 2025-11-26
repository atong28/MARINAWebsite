import types
import pytest


@pytest.fixture(scope="session")
def app_client():
    import app as app_module
    app_module._model_initialized = True
    app = app_module.app
    return app.test_client()


