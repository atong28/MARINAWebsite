from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.util import get_remote_address
from slowapi.errors import RateLimitExceeded


# Global limiter instance - will be initialized in setup_rate_limiting
limiter = None


def setup_rate_limiting(app):
    """Setup rate limiting for FastAPI app."""
    global limiter
    limiter = Limiter(key_func=get_remote_address)
    app.state.limiter = limiter
    app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)
    return limiter


def get_limiter():
    """Get the global limiter instance."""
    if limiter is None:
        raise RuntimeError(
            "Rate limiter not initialized. Call setup_rate_limiting(app) before importing routes."
        )
    return limiter

