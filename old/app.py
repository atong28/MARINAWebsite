import os
import logging
from logging.handlers import RotatingFileHandler
import pickle
import traceback
import numpy as np
import torch
import threading
import time
import base64
import io

from flask import Flask, jsonify, request, send_from_directory
try:
    from flask_compress import Compress
    _compress_available = True
except Exception:
    _compress_available = False
try:
    from flask_limiter import Limiter
    from flask_limiter.util import get_remote_address
    _rate_limit_available = True
except Exception:
    _rate_limit_available = False
from flask_cors import CORS

from src.fp_loader import EntropyFPLoader
from src import predictor
from src.ranker import RankingSet
from src.draw import draw_molecule, draw_fingerprint_changes, draw_similarity_comparison, get_fingerprint_differences, compute_bit_environments_batch
from src.services.metadata_service import MetadataService
from src.services.molecule_renderer import MoleculeRenderer
from src.services.result_builder import build_result_card
from src.services.model_service import ModelService
from src.config import MOLECULE_IMG_SIZE, METADATA_PATH, DEFAULT_TOP_K, MAX_TOP_K
from src.validators import validate_k, validate_smiles
from src.spec import get_openapi_spec

# Try to import RDKit, but don't fail if not available
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

# Register JSON error handlers
try:
    from src.http import register_error_handlers
    register_error_handlers(app)
except Exception:
    pass

# Optional request path logging for 404 tracing (debug only)
try:
    from src.config import HIGHLIGHT_DEBUG
    if HIGHLIGHT_DEBUG:
        @app.before_request
        def _log_request_path():
            try:
                logger.info(f"REQ {request.method} {request.path}")
            except Exception:
                pass
except Exception:
    pass

# Enable compression if available
if _compress_available:
    try:
        Compress(app)
    except Exception:
        pass

# Enable simple rate limiting if available
if _rate_limit_available:
    try:
        limiter = Limiter(get_remote_address, app=app, default_limits=["100 per minute"])  # coarse default
    except Exception:
        limiter = None
else:
    limiter = None

# Add security headers to enable clipboard access
@app.after_request
def after_request(response):
    # Enable clipboard access by setting appropriate headers (use only Permissions-Policy)
    response.headers['Permissions-Policy'] = 'clipboard-read=*, clipboard-write=*'
    # Allow clipboard access in iframe if needed
    response.headers['X-Frame-Options'] = 'SAMEORIGIN'
    return response

# Initialize production logging (minimal noise, write to logs/)
os.makedirs("logs", exist_ok=True)

# Root logger: warnings and above only
logging.getLogger().setLevel(logging.WARNING)

# Application logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)

# Rotating file handler for app logs
app_log_path = os.path.join("logs", "app.log")
if not any(isinstance(h, RotatingFileHandler) for h in logger.handlers):
    file_handler = RotatingFileHandler(app_log_path, maxBytes=2_000_000, backupCount=5)
    file_handler.setLevel(logging.WARNING)
    formatter = logging.Formatter(
        fmt='%(asctime)s %(levelname)s %(name)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

# Quiet down werkzeug/flask request logs in production
logging.getLogger("werkzeug").setLevel(logging.WARNING)

# Global model initialization flag
_model_initialized = False
_model_init_error = None
_model_loading = False
_model_init_attempted = False

# Server startup time for debugging
_server_start_time = time.time()

def initialize_model_background():
    """Initialize the model in a background thread to avoid blocking server startup."""
    global _model_initialized, _model_init_error, _model_loading
    
    # If predictor already has a model, mark initialized and return
    try:
        if getattr(predictor, "_model", None) is not None:
            _model_initialized = True
            _model_loading = False
            return
    except Exception:
        pass

    if _model_initialized or _model_loading:
        return
    
    _model_loading = True
    
    def load_model_worker():
        global _model_initialized, _model_init_error, _model_loading
        try:
            logger.info("Starting background model initialization...")
            predictor.load_model()
            _model_initialized = True
            _model_loading = False
            logger.info("Background model initialization completed successfully!")
        except Exception as e:
            _model_init_error = str(e)
            _model_loading = False
            logger.error(f"Background model initialization failed: {e}")
    
    # Start model loading in background thread
    thread = threading.Thread(target=load_model_worker, daemon=True)
    thread.start()
    logger.info("Model initialization started in background thread")

def initialize_model():
    """Initialize the model synchronously (for compatibility)."""
    global _model_initialized, _model_init_error, _model_init_attempted
    
    if _model_initialized:
        return True
    
    if _model_init_attempted:
        logger.info("Model initialization already attempted, skipping...")
        return False
    
    _model_init_attempted = True
    
    try:
        logger.info("Initializing MARINA model synchronously...")
        predictor.load_model()
        _model_initialized = True
        logger.info("Model initialization completed successfully!")
        return True
    except Exception as e:
        _model_init_error = str(e)
        logger.error(f"Model initialization failed: {e}")
        return False

# Initialize model once at startup (import time). Idempotent via predictor._model check.
# Only initialize if not already initialized to prevent multiple initializations
if not _model_initialized and not _model_loading:
    initialize_model()

# Preload metadata cache on startup (non-blocking)
try:
    threading.Thread(target=lambda: MetadataService.instance().preload(), daemon=True).start()
except Exception:
    pass


def pil_image_to_base64(image, format='PNG'):
    """Convert PIL Image to base64 string for web display."""
    buffer = io.BytesIO()
    image.save(buffer, format=format)
    img_str = base64.b64encode(buffer.getvalue()).decode()
    return f"data:image/{format.lower()};base64,{img_str}"



@app.route('/health', methods=['GET'])
def health():
    """Health check endpoint with model status."""
    global _model_initialized, _model_init_error, _model_loading, _model_init_attempted, _server_start_time
    
    # Log health check for debugging
    uptime = time.time() - _server_start_time
    logger.debug(f"Health check - initialized: {_model_initialized}, loading: {_model_loading}, attempted: {_model_init_attempted}, uptime: {uptime:.1f}s")
    
    base_response = {
        'uptime_seconds': round(uptime, 1),
        'server_start_time': _server_start_time
    }
    
    if _model_initialized:
        return jsonify({
            **base_response,
            'status': 'ok',
            'model_loaded': True,
            'message': 'Model is ready for predictions'
        })
    elif _model_loading:
        return jsonify({
            **base_response,
            'status': 'loading',
            'model_loaded': False,
            'message': 'Model is loading in background'
        })
    elif _model_init_error:
        return jsonify({
            **base_response,
            'status': 'error',
            'model_loaded': False,
            'error': _model_init_error,
            'message': 'Model initialization failed'
        }), 503
    elif _model_init_attempted:
        return jsonify({
            **base_response,
            'status': 'initializing',
            'model_loaded': False,
            'message': 'Model initialization in progress'
        }), 202
    else:
        return jsonify({
            **base_response,
            'status': 'initializing',
            'model_loaded': False,
            'message': 'Model initialization not started'
        }), 202


@app.route('/', methods=['GET'])
def index():
  # Serve the simple frontend
  return send_from_directory('static', 'index.html')
@app.route('/openapi.json', methods=['GET'])
def openapi_json():
    base = request.host_url.rstrip('/')
    return jsonify(get_openapi_spec(base))

@app.route('/docs', methods=['GET'])
def docs():
    return send_from_directory('static', 'docs.html')






## Removed /similarity endpoint for MVP (frontend uses top-k scores already)

# Apply route error decorator if available
try:
    _route_errors = app.route_errors
except Exception:
    _route_errors = lambda f: f


@app.route('/predict', methods=['POST'])
@_route_errors
@(limiter.limit("30 per minute") if limiter else (lambda f: f))
def predict():
    
    # Check if model is ready
    global _model_initialized, _model_init_error, _model_loading
    if not _model_initialized:
        if _model_loading:
            return jsonify({'error': 'Model is still loading in background. Please try again in a moment.'}), 503
        elif _model_init_error:
            return jsonify({'error': f'Model not available: {_model_init_error}'}), 503
        else:
            return jsonify({'error': 'Model initialization not started. Please try again in a moment.'}), 503
    
    try:
        payload = request.get_json()
    except Exception as e:
        return jsonify({'error': 'expected JSON body'}), 400

    if payload is None:
        return jsonify({'error': 'expected JSON body'}), 400

    k, k_err = validate_k(payload.get('k', DEFAULT_TOP_K))
    if k_err:
        return jsonify({'error': k_err}), 400
    if k <= 0:
        return jsonify({'error': 'k must be positive integer'}), 400
    if k > MAX_TOP_K:
        k = MAX_TOP_K

    try:
        # Model is already loaded from startup
        
        if 'raw' in payload:
            raw_data = payload['raw']
            out = predictor.predict_from_raw(raw_data, k=k)
            
            if isinstance(out, tuple) and len(out) >= 2:
                scores, indices = out[0], out[1]
                pred_fp = out[2] if len(out) > 2 else None
                
                # Get SMILES and molecular structures for all results
                results = []
                for i, idx in enumerate(indices):
                    entry = MetadataService.instance().get_entry(idx)
                    if entry:
                        result_smiles = MetadataService.instance().get_smiles(idx)
                        if result_smiles:
                            try:
                                fp_loader = ModelService.instance().get_fp_loader()
                                pred_tensor = torch.tensor(pred_fp, dtype=torch.float32) if pred_fp is not None else None
                                enhanced_svg = MoleculeRenderer.instance().render(
                                    result_smiles, predicted_fp=pred_tensor, fp_loader=fp_loader, img_size=MOLECULE_IMG_SIZE
                                )
                            except Exception:
                                enhanced_svg = None
                            similarity_val = scores[i] if scores[i] is not None else 0.0
                            try:
                                similarity_float = float(similarity_val)
                                if not isinstance(similarity_float, (int, float)) or similarity_float != similarity_float:
                                    similarity_float = 0.0
                            except (ValueError, TypeError):
                                similarity_float = 0.0
                            card = build_result_card(idx, entry, similarity_float, enhanced_svg)
                            # Also provide plain (non-highlight) depiction for analysis header
                            try:
                                plain_svg = MoleculeRenderer.instance().render(result_smiles, predicted_fp=None, fp_loader=None, img_size=MOLECULE_IMG_SIZE)
                                if plain_svg:
                                    card['plain_svg'] = plain_svg
                            except Exception:
                                pass
                            # Attach highlight flags in debug mode
                            try:
                                from src.config import HIGHLIGHT_DEBUG
                                if HIGHLIGHT_DEBUG:
                                    dbg = MoleculeRenderer.instance().last_debug()
                                    if dbg and isinstance(dbg, dict):
                                        card['highlighted'] = bool(dbg.get('highlighted'))
                                        card['highlight_debug'] = {k: dbg.get(k) for k in ('branch','reason')}
                            except Exception:
                                pass
                            # Attach retrieved molecule bit indices for instant analysis indices display
                            retrieved_indices = []
                            try:
                                from src.predictor import _rankingset_cache as _rs
                                if _rs is not None:
                                    if _rs.layout == torch.sparse_csr:
                                        row_start = _rs.crow_indices()[idx]
                                        row_end = _rs.crow_indices()[idx + 1]
                                        if row_start < row_end:
                                            col_indices = _rs.col_indices()[row_start:row_end]
                                            retrieved_indices = col_indices.cpu().tolist()
                                        else:
                                            retrieved_indices = []
                                    else:
                                        row_tensor = _rs[idx].cpu()
                                        retrieved_indices = torch.nonzero(row_tensor > 0.5, as_tuple=False).squeeze(-1).tolist()
                            except Exception:
                                pass
                            card['retrieved_molecule_fp_indices'] = retrieved_indices
                            
                            # Compute bit environments for all retrieved indices
                            try:
                                if retrieved_indices and fp_loader:
                                    bit_environments = compute_bit_environments_batch(result_smiles, retrieved_indices, fp_loader)
                                    card['bit_environments'] = bit_environments
                                else:
                                    card['bit_environments'] = {}
                            except Exception as e:
                                logger.warning(f"Failed to compute bit environments for result {idx}: {e}")
                                card['bit_environments'] = {}
                            
                            results.append(card)

                resp = {'results': results, 'total_count': len(results), 'offset': int(payload.get('offset', 0)), 'limit': int(payload.get('limit', len(results)))}
                # Apply pagination slice on results
                try:
                    offset = max(0, int(payload.get('offset', 0)))
                    limit = int(payload.get('limit', len(results)))
                    if limit < 0:
                        limit = len(results)
                    resp['results'] = results[offset: offset + limit]
                except Exception:
                    pass
                if pred_fp is not None:
                    resp['pred_fp'] = pred_fp
                from src.http import json_success
                return json_success(resp)
            else:
                logger.error(f"Unexpected predictor output: {out}")
                return jsonify({'error': 'unexpected predictor output'}), 500
        else:
            return jsonify({'error': "provide 'raw' in request body"}), 400
    except Exception as e:
        logger.error(f"Prediction error: {e}", exc_info=True)
        return jsonify({'error': str(e)}), 500



@app.route('/smiles-search', methods=['POST'])
@(limiter.limit("20 per minute") if limiter else (lambda f: f))
def smiles_search():
    """Search for similar molecules using a SMILES string."""
    global _model_initialized, _model_init_error, _model_loading
    if not _model_initialized:
        if _model_loading:
            return jsonify({'error': 'Model is still loading in background. Please try again in a moment.'}), 503
        elif _model_init_error:
            return jsonify({'error': f'Model not available: {_model_init_error}'}), 503
        else:
            return jsonify({'error': 'Model initialization not started. Please try again in a moment.'}), 503
    
    try:
        payload = request.get_json()
        if not payload or 'smiles' not in payload:
            return jsonify({'error': 'SMILES string required'}), 400
        
        smiles, smiles_err = validate_smiles(payload.get('smiles'))
        if smiles_err:
            return jsonify({'error': smiles_err}), 400
        k_raw = payload.get('k', DEFAULT_TOP_K)
        k, k_err = validate_k(k_raw)
        if k_err:
            return jsonify({'error': k_err}), 400
        if k <= 0:
            return jsonify({'error': 'k must be positive integer'}), 400
        if k > MAX_TOP_K:
            k = MAX_TOP_K

        if not smiles:
            return jsonify({'error': 'SMILES string cannot be empty'}), 400
        
        logger.info(f"SMILES search for: {smiles} with k={k}")
        
        # Use the existing fp_loader from predictor (reuse to save memory)
        # Get fingerprint for SMILES using the existing fp_loader
        try:
            # Access the global fp_loader from predictor module
            from src.predictor import _fp_loader
            if _fp_loader is None:
                return jsonify({'error': 'Model not loaded. Please ensure the model is initialized first.'}), 503
            
            fp = _fp_loader.build_mfp_for_smiles(smiles)
            if fp is None:
                return jsonify({'error': 'Invalid SMILES string or unable to generate fingerprint'}), 400
        except Exception as e:
            logger.error(f"Error converting SMILES to fingerprint: {e}")
            return jsonify({'error': 'Failed to process SMILES string. Please check the SMILES format.'}), 400
        
        # Use the same ranking approach as the predictor
        
        # Convert to tensor (don't apply sigmoid - rankingset is L2-normalized binary fingerprints)
        pred = torch.tensor(fp, dtype=torch.float)
        
        # Use RankingSet to get results (reuse cached store to save memory)
        from src.predictor import _rankingset_cache
        if _rankingset_cache is None:
            # If cache is not available, load it (this should rarely happen)
            store = _fp_loader.load_rankingset("RankingEntropy")
        else:
            store = _rankingset_cache
        ranker = RankingSet(store=store, metric="cosine")
        
        # Get top-k results
        idxs = ranker.retrieve_idx(pred.unsqueeze(0), n=k)
        
        # Get similarity scores for the retrieved indices
        sims = ranker._sims(pred.unsqueeze(0))  # (N, 1)
        sims_sorted, _ = torch.topk(sims.squeeze(), k=k, dim=0)
        
        # Get SMILES for all results
        results = []
        for i, idx in enumerate(idxs.squeeze().tolist()):
            entry = MetadataService.instance().get_entry(idx)
            if entry:
                result_smiles = MetadataService.instance().get_smiles(idx)
                if result_smiles:
                    # Render molecule (enhanced if possible)
                    try:
                        fp_loader = ModelService.instance().get_fp_loader()
                        query_tensor = torch.tensor(fp, dtype=torch.float32)
                        enhanced_svg = MoleculeRenderer.instance().render(
                            result_smiles, predicted_fp=query_tensor, fp_loader=fp_loader, img_size=MOLECULE_IMG_SIZE
                        )
                    except Exception:
                        enhanced_svg = None
                    # Similarity value
                    similarity_val = sims_sorted[i] if sims_sorted[i] is not None else 0.0
                    try:
                        similarity_float = float(similarity_val)
                        if not isinstance(similarity_float, (int, float)) or similarity_float != similarity_float:  # NaN guard
                            similarity_float = 0.0
                    except (ValueError, TypeError):
                        similarity_float = 0.0
                    # Build card
                    card = build_result_card(idx, entry, similarity_float, enhanced_svg)
                    try:
                        plain_svg = MoleculeRenderer.instance().render(result_smiles, predicted_fp=None, fp_loader=None, img_size=MOLECULE_IMG_SIZE)
                        if plain_svg:
                            card['plain_svg'] = plain_svg
                    except Exception:
                        pass
                    try:
                        from src.config import HIGHLIGHT_DEBUG
                        if HIGHLIGHT_DEBUG:
                            dbg = MoleculeRenderer.instance().last_debug()
                            if dbg and isinstance(dbg, dict):
                                card['highlighted'] = bool(dbg.get('highlighted'))
                                card['highlight_debug'] = {k: dbg.get(k) for k in ('branch','reason')}
                    except Exception:
                        pass
                    
                    # Attach retrieved molecule bit indices and compute bit environments
                    retrieved_indices = []
                    try:
                        if store is not None:
                            if store.layout == torch.sparse_csr:
                                row_start = store.crow_indices()[idx]
                                row_end = store.crow_indices()[idx + 1]
                                if row_start < row_end:
                                    col_indices = store.col_indices()[row_start:row_end]
                                    retrieved_indices = col_indices.cpu().tolist()
                                else:
                                    retrieved_indices = []
                            else:
                                row_tensor = store[idx].cpu()
                                retrieved_indices = torch.nonzero(row_tensor > 0.5, as_tuple=False).squeeze(-1).tolist()
                    except Exception:
                        pass
                    card['retrieved_molecule_fp_indices'] = retrieved_indices
                    
                    # Compute bit environments for all retrieved indices
                    try:
                        if retrieved_indices and fp_loader:
                            bit_environments = compute_bit_environments_batch(result_smiles, retrieved_indices, fp_loader)
                            card['bit_environments'] = bit_environments
                        else:
                            card['bit_environments'] = {}
                    except Exception as e:
                        logger.warning(f"Failed to compute bit environments for SMILES search result {idx}: {e}")
                        card['bit_environments'] = {}
                    
                    results.append(card)

        # Pagination
        offset = max(0, int(payload.get('offset', 0)))
        limit = int(payload.get('limit', len(results)))
        if limit < 0:
            limit = len(results)
        paged = results[offset: offset + limit]
        from src.http import json_success
        return json_success({'results': paged, 'total_count': len(results), 'offset': offset, 'limit': limit, 'query_smiles': smiles})
        
    except Exception as e:
        logger.error(f"SMILES search error: {e}", exc_info=True)
        return jsonify({'error': str(e)}), 500


@app.route('/analyze', methods=['POST'])
@(limiter.limit("10 per minute") if limiter else (lambda f: f))
def analyze():
    """Analysis endpoint for detailed molecular fingerprint analysis."""
    global _model_initialized, _model_init_error, _model_loading
    if not _model_initialized:
        if _model_loading:
            return jsonify({'error': 'Model is still loading in background. Please try again in a moment.'}), 503
        elif _model_init_error:
            return jsonify({'error': f'Model not available: {_model_init_error}'}), 503
        else:
            return jsonify({'error': 'Model initialization not started. Please try again in a moment.'}), 503
    
    try:
        payload = request.get_json()
        if not payload:
            return jsonify({'error': 'JSON body required'}), 400
        
        # Extract required fields
        current_data = payload.get('current_data')
        original_data = payload.get('original_data')
        target_smiles = payload.get('target_smiles')
        target_index = payload.get('target_index')
        
        if not all([current_data, original_data, target_smiles, target_index is not None]):
            return jsonify({'error': 'Missing required fields: current_data, original_data, target_smiles, target_index'}), 400
        
        logger.info(f"Analysis request for SMILES: {target_smiles}, index: {target_index}")
        
        # Get fingerprints for both datasets
        try:
            # Current prediction
            current_out = predictor.predict_from_raw(current_data, k=1)
            if isinstance(current_out, tuple) and len(current_out) >= 2:
                current_pred_fp = current_out[2] if len(current_out) > 2 else None
            else:
                return jsonify({'error': 'Failed to generate current fingerprint'}), 500
            
            # Original prediction
            original_out = predictor.predict_from_raw(original_data, k=1)
            if isinstance(original_out, tuple) and len(original_out) >= 2:
                original_pred_fp = original_out[2] if len(original_out) > 2 else None
            else:
                return jsonify({'error': 'Failed to generate original fingerprint'}), 500
            
            if current_pred_fp is None or original_pred_fp is None:
                return jsonify({'error': 'Failed to generate fingerprints for analysis'}), 500
            
            # Get retrieval molecule fingerprint from rankingset
            retrieved_molecule_fp = None
            retrieved_molecule_fp_indices = None
            try:
                from src.predictor import _rankingset_cache
                if _rankingset_cache is not None:
                    # Extract fingerprint row for target_index
                    if _rankingset_cache.layout == torch.sparse_csr:
                        # For sparse CSR, extract row
                        row_start = _rankingset_cache.crow_indices()[target_index]
                        row_end = _rankingset_cache.crow_indices()[target_index + 1]
                        if row_start < row_end:
                            col_indices = _rankingset_cache.col_indices()[row_start:row_end]
                            values = _rankingset_cache.values()[row_start:row_end]
                            # Build dense fingerprint
                            retrieved_fp_dense = torch.zeros(_rankingset_cache.size(1), dtype=torch.float32)
                            retrieved_fp_dense[col_indices] = values
                            retrieved_molecule_fp = retrieved_fp_dense.tolist()
                            # Extract indices (bit positions that are set)
                            retrieved_molecule_fp_indices = col_indices.cpu().tolist()
                        else:
                            retrieved_molecule_fp = [0.0] * _rankingset_cache.size(1)
                            retrieved_molecule_fp_indices = []
                    else:
                        # For dense
                        retrieved_fp_tensor = _rankingset_cache[target_index].cpu()
                        retrieved_molecule_fp = retrieved_fp_tensor.tolist()
                        # Extract indices where value is above threshold (0.5 for binary-like)
                        retrieved_molecule_fp_indices = torch.nonzero(retrieved_fp_tensor > 0.5, as_tuple=False).squeeze(-1).tolist()
            except Exception as e:
                logger.warning(f"Failed to get retrieval molecule fingerprint: {e}")
                retrieved_molecule_fp = None
                retrieved_molecule_fp_indices = None
            
            # Extract predicted fingerprint indices
            predicted_fp_indices = None
            if current_pred_fp:
                try:
                    pred_tensor = torch.tensor(current_pred_fp, dtype=torch.float32)
                    # For probabilistic fingerprints, use threshold of 0.5 or find max values
                    # Check if it's binary-like (all values are 0 or 1) or probabilistic
                    max_val = pred_tensor.max().item()
                    if max_val <= 1.0:
                        # Use 0.5 threshold for probabilistic
                        predicted_fp_indices = torch.nonzero(pred_tensor > 0.5, as_tuple=False).squeeze(-1).tolist()
                    else:
                        # Use max value threshold
                        threshold = max_val * 0.5
                        predicted_fp_indices = torch.nonzero(pred_tensor > threshold, as_tuple=False).squeeze(-1).tolist()
                except Exception as e:
                    logger.warning(f"Failed to extract predicted fingerprint indices: {e}")
                    predicted_fp_indices = None
            
        except Exception as e:
            logger.error(f"Error generating fingerprints: {e}")
            return jsonify({'error': f'Failed to generate fingerprints: {str(e)}'}), 500
        
        # Generate visualizations and prepare response
        result = {
            'target_smiles': target_smiles,
            'target_index': target_index,
            'predicted_fp': current_pred_fp,
            'retrieved_molecule_fp': retrieved_molecule_fp,
            'predicted_fp_indices': predicted_fp_indices,
            'retrieved_molecule_fp_indices': retrieved_molecule_fp_indices
        }
        
        try:
            from src.predictor import _fp_loader
            if _fp_loader is not None and RDKIT_AVAILABLE:
                # Generate similarity highlighting for current prediction
                current_tensor = torch.tensor(current_pred_fp, dtype=torch.float)
                original_tensor = torch.tensor(original_pred_fp, dtype=torch.float)
                
                # Similarity visualization (current prediction)
                similarity_img = draw_similarity_comparison(
                    predicted_fp=current_tensor,
                    retrieval_smiles=target_smiles,
                    fp_loader=_fp_loader,
                    need_to_clean_H=False,
                    img_size=MOLECULE_IMG_SIZE
                )
                result['similarity_visualization'] = pil_image_to_base64(similarity_img)
                
                # Change visualization (difference between original and current)
                change_img = draw_fingerprint_changes(
                    original_fp=original_tensor,
                    new_fp=current_tensor,
                    retrieval_smiles=target_smiles,
                    fp_loader=_fp_loader,
                    need_to_clean_H=False,
                    img_size=MOLECULE_IMG_SIZE
                )
                result['change_visualization'] = pil_image_to_base64(change_img)
                
                # Get fingerprint differences
                fingerprint_diffs = get_fingerprint_differences(
                    original_fp=original_tensor,
                    new_fp=current_tensor,
                    retrieval_smiles=target_smiles,
                    fp_loader=_fp_loader,
                    need_to_clean_H=False
                )
                
                
                result['fingerprint_differences'] = fingerprint_diffs
            else:
                result['similarity_visualization'] = None
                result['change_visualization'] = None
                result['fingerprint_differences'] = None
                
        except Exception as e:
            logger.warning(f"Failed to generate analysis visualizations: {e}")
            result['similarity_visualization'] = None
            result['change_visualization'] = None
        
        # Secondary retrieval computed within analyze (avoid extra API call)
        try:
            # Only proceed if both fingerprints are present
            if current_pred_fp is not None and retrieved_molecule_fp is not None:
                pred_tensor = torch.tensor(current_pred_fp, dtype=torch.float32)
                # Build retrieved tensor from list
                retrieved_tensor = torch.tensor(retrieved_molecule_fp, dtype=torch.float32)
                
                # Clamp to [0,1]
                pred_tensor = torch.clamp(pred_tensor, 0.0, 1.0)
                retrieved_tensor = torch.clamp(retrieved_tensor, 0.0, 1.0)
                
                # Compute overlap and difference
                overlap = torch.min(pred_tensor, retrieved_tensor)
                difference_fp = torch.clamp(pred_tensor - overlap, 0.0, 1.0)
                
                if difference_fp.sum() < 1e-6:
                    # Equal fingerprints (no remaining bits)
                    result['secondary_skipped'] = True
                    result['secondary_message'] = 'Predicted fingerprint and selected molecule fingerprint are identical. No remaining bits to search.'
                else:
                    # Run retrieval on the difference fingerprint
                    from src.predictor import _rankingset_cache, _fp_loader, _args
                    if _rankingset_cache is None:
                        if _fp_loader is None:
                            raise RuntimeError('Model not loaded')
                        _rankingset_cache = _fp_loader.load_rankingset(_args.fp_type)
                    ranker = RankingSet(store=_rankingset_cache, metric="cosine")
                    idxs = ranker.retrieve_idx(difference_fp.unsqueeze(0), n=10)
                    sims = ranker._sims(difference_fp.unsqueeze(0))
                    sims_sorted, _ = torch.topk(sims.squeeze(), k=min(10, sims.numel()), dim=0)
                    
                    # Load metadata
                    meta_path = os.path.join('data', 'metadata.json')
                    if not os.path.exists(meta_path):
                        raise FileNotFoundError('metadata.json not found')
                    import json
                    with open(meta_path, 'r') as f:
                        meta = json.load(f)
                    
                    secondary_results = []
                    for i, idx in enumerate(idxs.squeeze().tolist()):
                        entry = MetadataService.instance().get_entry(idx)
                        if not entry:
                            continue
                        result_smiles = MetadataService.instance().get_smiles(idx)
                        # Render molecule using services
                        try:
                            fp_loader = ModelService.instance().get_fp_loader()
                            enhanced_svg = MoleculeRenderer.instance().render(
                                result_smiles, predicted_fp=difference_fp, fp_loader=fp_loader, img_size=MOLECULE_IMG_SIZE
                            )
                        except Exception:
                            enhanced_svg = None
                        # Similarity
                        similarity_val = sims_sorted[i] if sims_sorted.numel() > i else 0.0
                        try:
                            similarity_float = float(similarity_val)
                            if not isinstance(similarity_float, (int, float)) or similarity_float != similarity_float:
                                similarity_float = 0.0
                        except (ValueError, TypeError):
                            similarity_float = 0.0
                        # Append card
                        secondary_results.append(build_result_card(idx, entry, similarity_float, enhanced_svg))
                    result['secondary_results'] = secondary_results
                    result['difference_fp'] = difference_fp.tolist()
        except Exception as e:
            logger.warning(f"Failed to compute secondary retrieval in analyze: {e}")

        from src.http import json_success
        return json_success(result)
        
    except Exception as e:
        logger.error(f"Analysis error: {e}", exc_info=True)
        return jsonify({'error': str(e)}), 500


@app.route('/secondary-retrieval', methods=['POST'])
@(limiter.limit("20 per minute") if limiter else (lambda f: f))
def secondary_retrieval():
    """Secondary retrieval endpoint using difference fingerprint (predicted - overlap)."""
    global _model_initialized, _model_init_error, _model_loading
    if not _model_initialized:
        if _model_loading:
            return jsonify({'error': 'Model is still loading in background. Please try again in a moment.'}), 503
        elif _model_init_error:
            return jsonify({'error': f'Model not available: {_model_init_error}'}), 503
        else:
            return jsonify({'error': 'Model initialization not started. Please try again in a moment.'}), 503
    
    try:
        payload = request.get_json()
        if not payload:
            return jsonify({'error': 'JSON body required'}), 400
        
        predicted_fp = payload.get('predicted_fp')
        retrieved_fp = payload.get('retrieved_fp')
        k, k_err = validate_k(payload.get('k', DEFAULT_TOP_K))
        if k_err:
            return jsonify({'error': k_err}), 400
        if k <= 0:
            return jsonify({'error': 'k must be positive integer'}), 400
        if k > MAX_TOP_K:
            k = MAX_TOP_K

        if not predicted_fp or not retrieved_fp:
            return jsonify({'error': 'Missing required fields: predicted_fp, retrieved_fp'}), 400
        
        # Convert to tensors
        predicted_tensor = torch.tensor(predicted_fp, dtype=torch.float32)
        retrieved_tensor = torch.tensor(retrieved_fp, dtype=torch.float32)
        
        # Ensure fingerprints are in valid range [0, 1]
        predicted_tensor = torch.clamp(predicted_tensor, 0.0, 1.0)
        retrieved_tensor = torch.clamp(retrieved_tensor, 0.0, 1.0)
        
        # Compute overlap (element-wise minimum for binary/probabilistic fingerprints)
        overlap = torch.min(predicted_tensor, retrieved_tensor)
        
        # Compute difference: predicted_fp - overlap
        difference_fp = predicted_tensor - overlap
        difference_fp = torch.clamp(difference_fp, 0.0, 1.0)
        
        # Check if difference fingerprint is all zeros
        if difference_fp.sum() < 1e-6:
            logger.warning("Difference fingerprint is all zeros, returning empty results")
            return jsonify({'results': [], 'difference_fp': difference_fp.tolist()})
        
        # Use difference fingerprint directly (already in [0,1] range)
        # RankingSet will handle normalization internally
        difference_prob = difference_fp
        
        # Use the existing RankingSet to retrieve results
        from src.predictor import _rankingset_cache, _fp_loader
        if _rankingset_cache is None:
            if _fp_loader is None:
                return jsonify({'error': 'Model not loaded'}), 503
            from src.predictor import _args
            _rankingset_cache = _fp_loader.load_rankingset(_args.fp_type)
        
        ranker = RankingSet(store=_rankingset_cache, metric="cosine")
        
        # Retrieve top-k indices using difference fingerprint
        idxs = ranker.retrieve_idx(difference_prob.unsqueeze(0), n=k)
        
        # Get similarity scores for the retrieved indices
        sims = ranker._sims(difference_prob.unsqueeze(0))
        sims_sorted, _ = torch.topk(sims.squeeze(), k=k, dim=0)
        
        # Build results using services
        results = []
        for i, idx in enumerate(idxs.squeeze().tolist()):
            entry = MetadataService.instance().get_entry(idx)
            if not entry:
                continue
            result_smiles = MetadataService.instance().get_smiles(idx)
            # Render with services
            try:
                fp_loader = ModelService.instance().get_fp_loader()
                enhanced_svg = MoleculeRenderer.instance().render(
                    result_smiles, predicted_fp=difference_prob, fp_loader=fp_loader, img_size=MOLECULE_IMG_SIZE
                )
            except Exception:
                enhanced_svg = None
            # Similarity
            similarity_val = sims_sorted[i] if sims_sorted[i] is not None else 0.0
            try:
                similarity_float = float(similarity_val)
                if not isinstance(similarity_float, (int, float)) or similarity_float != similarity_float:
                    similarity_float = 0.0
            except (ValueError, TypeError):
                similarity_float = 0.0
            results.append(build_result_card(idx, entry, similarity_float, enhanced_svg))
        
        # Pagination
        offset = max(0, int(payload.get('offset', 0)))
        limit = int(payload.get('limit', len(results)))
        if limit < 0:
            limit = len(results)
        paged = results[offset: offset + limit]
        from src.http import json_success
        return json_success({'results': paged, 'total_count': len(results), 'offset': offset, 'limit': limit, 'difference_fp': difference_fp.tolist()})
        
    except Exception as e:
        logger.error(f"Secondary retrieval error: {e}", exc_info=True)
        return jsonify({'error': str(e)}), 500
