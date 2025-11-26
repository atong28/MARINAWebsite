import os
from pathlib import Path

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
METADATA_PATH = os.path.join(DATA_DIR, "metadata.json")

# Model / fingerprint artifacts
CKPT_PATH = os.getenv("CKPT_PATH", os.path.join(DATA_DIR, "best.ckpt"))
PARAMS_PATH = os.getenv("PARAMS_PATH", os.path.join(DATA_DIR, "params.json"))
RETRIEVAL_PATH = os.getenv("RETRIEVAL_PATH", os.path.join(DATA_DIR, "retrieval.pkl"))

# Server
# Port configuration: Use BACKEND_PORT from root .env file (single source of truth)
PORT = int(os.getenv("BACKEND_PORT", "5000"))
HOST = os.getenv("BACKEND_HOST", "0.0.0.0")

