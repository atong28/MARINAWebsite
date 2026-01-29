import os

# General
DEFAULT_TOP_K = int(os.getenv("DEFAULT_TOP_K", "10"))
MAX_TOP_K = int(os.getenv("MAX_TOP_K", "50"))

# Fingerprints
FP_THRESHOLD = float(os.getenv("FP_THRESHOLD", "0.5"))

# Rendering
MOLECULE_IMG_SIZE = int(os.getenv("MOLECULE_IMG_SIZE", "400"))
HIGHLIGHT_DEBUG = os.getenv("HIGHLIGHT_DEBUG", "false").lower() == "true"

# RDKit
RDKIT_ENABLED = os.getenv("RDKIT_ENABLED", "true").lower() == "true"

# Paths
DATA_DIR = os.getenv("DATA_DIR", "data")
MODELS_JSON_PATH = os.getenv(
    "MODELS_JSON_PATH",
    os.path.join(DATA_DIR, "models.json"),
)
# Per-model root (e.g. data/marina_best). Fallback when models.json missing.
MODEL_ROOT = os.getenv("MODEL_ROOT", os.path.join(DATA_DIR, "marina_best"))
DEFAULT_MODEL_ID = os.getenv("DEFAULT_MODEL_ID", "marina_best")

# Model / fingerprint artifacts (derived from MODEL_ROOT)
CKPT_PATH = os.getenv("CKPT_PATH", os.path.join(MODEL_ROOT, "best.ckpt"))
PARAMS_PATH = os.getenv("PARAMS_PATH", os.path.join(MODEL_ROOT, "params.json"))
RETRIEVAL_PATH = os.getenv("RETRIEVAL_PATH", os.path.join(MODEL_ROOT, "retrieval.pkl"))
METADATA_PATH = os.getenv("METADATA_PATH", os.path.join(MODEL_ROOT, "metadata.json"))

# Server
# Port configuration: Use BACKEND_PORT from root .env file (single source of truth)
PORT = int(os.getenv("BACKEND_PORT", "5000"))
HOST = os.getenv("BACKEND_HOST", "0.0.0.0")

