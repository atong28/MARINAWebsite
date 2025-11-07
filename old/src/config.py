import os

# General
DEFAULT_TOP_K = int(os.getenv("DEFAULT_TOP_K", "10"))
MAX_TOP_K = int(os.getenv("MAX_TOP_K", "50"))

# Fingerprints
FP_THRESHOLD = float(os.getenv("FP_THRESHOLD", "0.5"))

# Rendering
MOLECULE_IMG_SIZE = int(os.getenv("MOLECULE_IMG_SIZE", "400"))
ENABLE_RENDERING_CACHE = os.getenv("ENABLE_RENDERING_CACHE", "true").lower() == "true"
RENDERING_CACHE_SIZE = int(os.getenv("RENDERING_CACHE_SIZE", "1000"))
HIGHLIGHT_DEBUG = os.getenv("HIGHLIGHT_DEBUG", "false").lower() == "true"

# RDKit
RDKIT_ENABLED = os.getenv("RDKIT_ENABLED", "true").lower() == "true"

# Paths
METADATA_PATH = os.getenv("METADATA_PATH", os.path.join("data", "metadata.json"))


