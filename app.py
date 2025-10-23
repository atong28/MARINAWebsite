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
from flask_cors import CORS

from src.fp_loader import EntropyFPLoader
from src import predictor
from src.ranker import RankingSet
from src.draw import draw_molecule, draw_fingerprint_changes, draw_similarity_comparison, get_fingerprint_differences

# Try to import RDKit, but don't fail if not available
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

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
    logger.info("Starting model initialization at import time...")
    initialize_model()
else:
    logger.info(f"Model initialization skipped - initialized: {_model_initialized}, loading: {_model_loading}")


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






## Removed /similarity endpoint for MVP (frontend uses top-k scores already)

# Minimal prediction endpoint (CPU backend)
@app.route('/predict', methods=['POST'])
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
        payload = request.get_json(force=True)
    except Exception as e:
        return jsonify({'error': 'expected JSON body'}), 400

    if payload is None:
        return jsonify({'error': 'expected JSON body'}), 400

    k = int(payload.get('k', 5))
    if k <= 0:
        return jsonify({'error': 'k must be positive integer'}), 400

    try:
        # Model is already loaded from startup
        
        if 'raw' in payload:
            raw_data = payload['raw']
            logger.info(f"Predicting from raw data with k={k}")
            logger.info(f"Raw data keys: {list(raw_data.keys())}")
            
            logger.info(f"About to call predict_from_raw with k={k} (type: {type(k)})")
            logger.info(f"raw_data type: {type(raw_data)}, keys: {list(raw_data.keys())}")
            out = predictor.predict_from_raw(raw_data, k=k)
            logger.info(f"Prediction completed successfully, output type: {type(out)}")
            
            if isinstance(out, tuple) and len(out) >= 2:
                scores, indices = out[0], out[1]
                pred_fp = out[2] if len(out) > 2 else None
                logger.info(f"Prediction completed: {len(scores)} scores and {len(indices)} indices")
                
                # Load metadata to get SMILES for all results
                meta_path = os.path.join('data', 'metadata.json')
                if not os.path.exists(meta_path):
                    return jsonify({'error': 'metadata.json not found'}), 404
                
                import json
                with open(meta_path, 'r') as f:
                    meta = json.load(f)
                
                # Get SMILES and molecular structures for all results
                results = []
                for i, idx in enumerate(indices):
                    entry = meta.get(str(idx))
                    if entry:
                        # Use canonical_3d_smiles for display (stereochemistry enabled)
                        result_smiles = entry['canonical_3d_smiles'] if entry['canonical_3d_smiles'] != 'N/A' else entry['canonical_2d_smiles']
                        if result_smiles:
                            # Generate enhanced molecular structure with similarity highlighting
                            try:
                                if RDKIT_AVAILABLE and pred_fp is not None:
                                    # Use the enhanced draw_molecule function
                                    from src.predictor import _fp_loader
                                    if _fp_loader is not None:
                                        pred_tensor = torch.tensor(pred_fp, dtype=torch.float)
                                        enhanced_img = draw_molecule(
                                            predicted_fp=pred_tensor,
                                            retrieval_smiles=result_smiles,
                                            fp_loader=_fp_loader,
                                            need_to_clean_H=False,
                                            img_size=400
                                        )
                                        # Convert to base64 for web display
                                        enhanced_svg = pil_image_to_base64(enhanced_img)
                                    else:
                                        # Fallback to basic SVG
                                        mol = Chem.MolFromSmiles(result_smiles)
                                        if mol:
                                            enhanced_svg = Draw.MolToSVG(mol, width=200, height=200)
                                        else:
                                            enhanced_svg = None
                                else:
                                    # Fallback to basic SVG
                                    if RDKIT_AVAILABLE:
                                        mol = Chem.MolFromSmiles(result_smiles)
                                        if mol:
                                            enhanced_svg = Draw.MolToSVG(mol, width=200, height=200)
                                        else:
                                            enhanced_svg = None
                                    else:
                                        enhanced_svg = None
                            except Exception as e:
                                logger.warning(f"Failed to render enhanced molecule for SMILES {result_smiles}: {e}")
                                # Fallback to basic SVG
                                try:
                                    if RDKIT_AVAILABLE:
                                        mol = Chem.MolFromSmiles(result_smiles)
                                        if mol:
                                            enhanced_svg = Draw.MolToSVG(mol, width=200, height=200)
                                        else:
                                            enhanced_svg = None
                                    else:
                                        enhanced_svg = None
                                except Exception as e2:
                                    logger.warning(f"Failed to render fallback molecule for SMILES {result_smiles}: {e2}")
                                    enhanced_svg = None
                            
                            # Ensure similarity is a valid number
                            similarity_val = scores[i] if scores[i] is not None else 0.0
                            try:
                                similarity_float = float(similarity_val)
                                if not isinstance(similarity_float, (int, float)) or similarity_float != similarity_float:  # Check for NaN
                                    similarity_float = 0.0
                            except (ValueError, TypeError):
                                similarity_float = 0.0
                            
                            # Extract database information
                            coconut = entry.get('coconut')
                            lotus = entry.get('lotus')
                            npmrd = entry.get('npmrd')
                            
                            # Determine primary name and link
                            primary_name = None
                            primary_link = None
                            
                            if coconut:
                                primary_name = coconut.get('name')
                                coconut_id = coconut.get('coconut_id')
                                if coconut_id:
                                    primary_link = f"https://coconut.naturalproducts.net/compounds/{coconut_id}"
                            elif lotus:
                                primary_name = lotus.get('name')
                                lotus_id = lotus.get('lotus_id')
                                if lotus_id:
                                    primary_link = f"https://lotus.naturalproducts.net/compound/lotus_id/{lotus_id}"
                            elif npmrd:
                                primary_name = npmrd.get('name')
                                npmrd_id = npmrd.get('npmrd_id')
                                if npmrd_id:
                                    primary_link = f"https://np-mrd.org/natural_products/{npmrd_id}"
                            
                            # Prepare database links
                            database_links = {}
                            if coconut and coconut.get('coconut_id'):
                                database_links['coconut'] = f"https://coconut.naturalproducts.net/compounds/{coconut['coconut_id']}"
                            if lotus and lotus.get('lotus_id'):
                                database_links['lotus'] = f"https://lotus.naturalproducts.net/compounds/{lotus['lotus_id']}"
                            if npmrd and npmrd.get('npmrd_id'):
                                database_links['npmrd'] = f"https://np-mrd.org/natural_products/{npmrd['npmrd_id']}"
                            
                            results.append({
                                'index': idx,
                                'smiles': result_smiles,
                                'similarity': similarity_float,
                                'svg': enhanced_svg,
                                'name': primary_name,
                                'primary_link': primary_link,
                                'database_links': database_links,
                                'np_pathway': None,  # Placeholder for future use
                                'np_superclass': None,  # Placeholder for future use
                                'np_class': None  # Placeholder for future use
                            })
                
                resp = {'results': results}
                if pred_fp is not None:
                    resp['pred_fp'] = pred_fp
                return jsonify(resp)
            else:
                logger.error(f"Unexpected predictor output: {out}")
                return jsonify({'error': 'unexpected predictor output'}), 500
        else:
            return jsonify({'error': "provide 'raw' in request body"}), 400
    except Exception as e:
        logger.error(f"Prediction error: {e}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500



@app.route('/smiles-search', methods=['POST'])
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
        
        smiles = payload['smiles'].strip()
        k = payload.get('k', 10)
        
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
        
        # Load metadata to get SMILES for all results
        meta_path = os.path.join('data', 'metadata.json')
        if not os.path.exists(meta_path):
            return jsonify({'error': 'metadata.json not found'}), 404
        
        import json
        with open(meta_path, 'r') as f:
            meta = json.load(f)
        
        # Get SMILES for all results
        results = []
        for i, idx in enumerate(idxs.squeeze().tolist()):
            entry = meta.get(str(idx))
            if entry:
                # Use canonical_3d_smiles for display (stereochemistry enabled)
                result_smiles = entry['canonical_3d_smiles'] if entry['canonical_3d_smiles'] != 'N/A' else entry['canonical_2d_smiles']
                if result_smiles:
                    # Generate enhanced molecular structure with similarity highlighting
                    try:
                        if RDKIT_AVAILABLE:
                            # Use the enhanced draw_molecule function
                            from src.predictor import _fp_loader
                            if _fp_loader is not None:
                                # Convert query fingerprint to tensor
                                query_tensor = torch.tensor(fp, dtype=torch.float)
                                enhanced_img = draw_molecule(
                                    predicted_fp=query_tensor,
                                    retrieval_smiles=result_smiles,
                                    fp_loader=_fp_loader,
                                    need_to_clean_H=False,
                                    img_size=400
                                )
                                # Convert to base64 for web display
                                enhanced_svg = pil_image_to_base64(enhanced_img)
                            else:
                                # Fallback to basic SVG
                                mol = Chem.MolFromSmiles(result_smiles)
                                if mol:
                                    enhanced_svg = Draw.MolToSVG(mol, width=200, height=200)
                                else:
                                    enhanced_svg = None
                        else:
                            enhanced_svg = None
                    except Exception as e:
                        logger.warning(f"Failed to render enhanced molecule for SMILES {result_smiles}: {e}")
                        # Fallback to basic SVG
                        try:
                            if RDKIT_AVAILABLE:
                                mol = Chem.MolFromSmiles(result_smiles)
                                if mol:
                                    enhanced_svg = Draw.MolToSVG(mol, width=200, height=200)
                                else:
                                    enhanced_svg = None
                            else:
                                enhanced_svg = None
                        except Exception as e2:
                            logger.warning(f"Failed to render fallback molecule for SMILES {result_smiles}: {e2}")
                            enhanced_svg = None
                    
                    # Ensure similarity is a valid number
                    similarity_val = sims_sorted[i] if sims_sorted[i] is not None else 0.0
                    try:
                        similarity_float = float(similarity_val)
                        if not isinstance(similarity_float, (int, float)) or similarity_float != similarity_float:  # Check for NaN
                            similarity_float = 0.0
                    except (ValueError, TypeError):
                        similarity_float = 0.0
                    
                    # Extract database information
                    coconut = entry.get('coconut')
                    lotus = entry.get('lotus')
                    npmrd = entry.get('npmrd')
                    
                    # Determine primary name and link
                    primary_name = None
                    primary_link = None
                    
                    if coconut:
                        primary_name = coconut.get('name')
                        coconut_id = coconut.get('coconut_id')
                        if coconut_id:
                            primary_link = f"https://coconut.naturalproducts.net/compounds/{coconut_id}"
                    elif lotus:
                        primary_name = lotus.get('name')
                        lotus_id = lotus.get('lotus_id')
                        if lotus_id:
                            primary_link = f"https://lotus.naturalproducts.net/compounds/{lotus_id}"
                    elif npmrd:
                        primary_name = npmrd.get('name')
                        npmrd_id = npmrd.get('npmrd_id')
                        if npmrd_id:
                            primary_link = f"https://np-mrd.org/natural_products/{npmrd_id}"
                    
                    # Prepare database links
                    database_links = {}
                    if coconut and coconut.get('coconut_id'):
                        database_links['coconut'] = f"https://coconut.naturalproducts.net/compounds/{coconut['coconut_id']}"
                    if lotus and lotus.get('lotus_id'):
                        database_links['lotus'] = f"https://lotus.naturalproducts.net/compounds/{lotus['lotus_id']}"
                    if npmrd and npmrd.get('npmrd_id'):
                        database_links['npmrd'] = f"https://np-mrd.org/natural_products/{npmrd['npmrd_id']}"
                    
                    results.append({
                        'index': idx,
                        'smiles': result_smiles,
                        'similarity': similarity_float,
                        'svg': enhanced_svg,
                        'name': primary_name,
                        'primary_link': primary_link,
                        'database_links': database_links,
                        'np_pathway': None,
                        'np_superclass': None,
                        'np_class': None
                    })
        
        return jsonify({
            'results': results,
            'query_smiles': smiles
        })
        
    except Exception as e:
        logger.error(f"SMILES search error: {e}", exc_info=True)
        return jsonify({'error': str(e)}), 500


@app.route('/analyze', methods=['POST'])
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
            
        except Exception as e:
            logger.error(f"Error generating fingerprints: {e}")
            return jsonify({'error': f'Failed to generate fingerprints: {str(e)}'}), 500
        
        # Generate visualizations
        result = {
            'target_smiles': target_smiles,
            'target_index': target_index
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
                    img_size=400
                )
                result['similarity_visualization'] = pil_image_to_base64(similarity_img)
                
                # Change visualization (difference between original and current)
                change_img = draw_fingerprint_changes(
                    original_fp=original_tensor,
                    new_fp=current_tensor,
                    retrieval_smiles=target_smiles,
                    fp_loader=_fp_loader,
                    need_to_clean_H=False,
                    img_size=400
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
                
                # Debug logging for fingerprint differences
                logger.info(f"=== FINGERPRINT DIFFERENCES DEBUG ===")
                logger.info(f"Target SMILES: {target_smiles}")
                logger.info(f"Original FP shape: {original_tensor.shape}, sum: {original_tensor.sum()}")
                logger.info(f"Current FP shape: {current_tensor.shape}, sum: {current_tensor.sum()}")
                logger.info(f"FP difference sum: {(current_tensor - original_tensor).sum()}")
                logger.info(f"Added bits count: {len(fingerprint_diffs.get('added', []))}")
                logger.info(f"Removed bits count: {len(fingerprint_diffs.get('removed', []))}")
                logger.info(f"Unchanged bits count: {len(fingerprint_diffs.get('unchanged', []))}")
                
                if fingerprint_diffs.get('added'):
                    logger.info(f"Added bits: {[bit['bit_id'] for bit in fingerprint_diffs['added'][:5]]}")
                if fingerprint_diffs.get('removed'):
                    logger.info(f"Removed bits: {[bit['bit_id'] for bit in fingerprint_diffs['removed'][:5]]}")
                logger.info(f"=== END FINGERPRINT DIFFERENCES DEBUG ===")
                
                result['fingerprint_differences'] = fingerprint_diffs
            else:
                result['similarity_visualization'] = None
                result['change_visualization'] = None
                result['fingerprint_differences'] = None
                
        except Exception as e:
            logger.warning(f"Failed to generate analysis visualizations: {e}")
            result['similarity_visualization'] = None
            result['change_visualization'] = None
        
        return jsonify(result)
        
    except Exception as e:
        logger.error(f"Analysis error: {e}", exc_info=True)
        return jsonify({'error': str(e)}), 500