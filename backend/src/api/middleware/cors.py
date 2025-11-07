from fastapi.middleware.cors import CORSMiddleware


def setup_cors(app):
    """Setup CORS middleware for FastAPI app."""
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],  # In production, specify actual origins
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

