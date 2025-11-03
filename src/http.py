import logging
from flask import jsonify


logger = logging.getLogger(__name__)


def json_error_response(message: str, status: int = 500):
    return jsonify({'success': False, 'error': message}), status


def json_success(data: dict):
    # Preserve backward compatibility by returning both wrapped and flat data
    # Callers can read either data[...] or the top-level keys
    wrapped = {'success': True, 'data': data}
    # Also merge flat keys (non-conflicting) for compatibility
    for k, v in data.items():
        if k not in wrapped:
            wrapped[k] = v
    return jsonify(wrapped)


def register_error_handlers(app):
    # Decorator to wrap route handlers with try/except and JSON errors
    def route_errors(fn):
        def wrapper(*args, **kwargs):
            try:
                return fn(*args, **kwargs)
            except Exception as e:
                try:
                    logger.exception("Route error")
                except Exception:
                    pass
                return json_error_response(str(e) or 'Internal server error', status=500)
        # Preserve function name for Flask routing
        wrapper.__name__ = getattr(fn, '__name__', 'wrapped')
        return wrapper

    app.route_errors = route_errors
    @app.errorhandler(Exception)
    def handle_unexpected_error(e):
        try:
            logger.exception("Unhandled exception")
        except Exception:
            pass
        return json_error_response(str(e) or 'Internal server error', status=500)


