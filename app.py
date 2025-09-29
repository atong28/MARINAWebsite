import os
import logging
from logging.handlers import RotatingFileHandler
import pickle
import traceback
import numpy as np
import torch
import threading
import time

from flask import Flask, jsonify, request, send_from_directory

from src.fp_loader import EntropyFPLoader
from src import predictor
from src.ranker import RankingSet

# Try to import RDKit, but don't fail if not available
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

app = Flask(__name__)

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
    global _model_initialized, _model_init_error
    
    if _model_initialized:
        return True
    
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
initialize_model()


@app.route('/health', methods=['GET'])
def health():
    """Health check endpoint with model status."""
    global _model_initialized, _model_init_error, _model_loading
    
    if _model_initialized:
        return jsonify({
            'status': 'ok',
            'model_loaded': True,
            'message': 'Model is ready for predictions'
        })
    elif _model_loading:
        return jsonify({
            'status': 'loading',
            'model_loaded': False,
            'message': 'Model is loading in background'
        })
    elif _model_init_error:
        return jsonify({
            'status': 'error',
            'model_loaded': False,
            'error': _model_init_error,
            'message': 'Model initialization failed'
        }), 503
    else:
        return jsonify({
            'status': 'initializing',
            'model_loaded': False,
            'message': 'Model initialization not started'
        }), 202


@app.route('/', methods=['GET'])
def index():
  # Serve the simple frontend
  return send_from_directory('static', 'index.html')


@app.route('/meta/<int:idx>', methods=['GET'])
def meta(idx: int):
  # Return metadata (smiles) for a retrieval index from rankingset_meta.pkl
  meta_path = os.path.join('data', 'rankingset_meta.pkl')
  if not os.path.exists(meta_path):
    return jsonify({'error': 'rankingset_meta.pkl not found; mount data/ containing it'}), 404
  with open(meta_path, 'rb') as f:
    meta = pickle.load(f)
  entry = meta.get(idx)
  if entry is None:
    return jsonify({'error': f'index {idx} not found in metadata'}), 404
  return jsonify({'smiles': entry.get('smiles')})


@app.route('/similarity', methods=['POST'])
def similarity():
  """Compute Tanimoto similarity between a predicted floating fingerprint and a SMILES-derived fingerprint.

  Request JSON: {'pred_fp': [...], 'smiles': '...'}
  """
  data = request.get_json(force=True)
  pred_fp = data.get('pred_fp')
  smiles = data.get('smiles')
  if pred_fp is None or smiles is None:
    return jsonify({'error': 'provide pred_fp and smiles in JSON body'}), 400

  # Try to compute fingerprint from SMILES using RDKit if available
  if not RDKIT_AVAILABLE:
    return jsonify({'error': 'RDKit required to compute SMILES fingerprint on server. Install rdkit or compute client-side.'}), 501

  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    return jsonify({'error': 'invalid smiles'}), 400
  # use Morgan fingerprint (radius=2) as bit vector
  fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=len(pred_fp))
  # convert bitvector to list of 0/1
  smi_fp = list(fp.ToBitString())
  smi_fp = [int(x) for x in smi_fp]

  # pred_fp may be float; binarize by threshold 0.5
  pred_bin = [1 if float(x) > 0.5 else 0 for x in pred_fp]

  # compute tanimoto
  a = np.array(pred_bin)
  b = np.array(smi_fp)
  inter = int(np.sum(a & b))
  union = int(np.sum((a | b)))
  tanimoto = inter / union if union > 0 else 0.0
  return jsonify({'tanimoto': float(tanimoto), 'intersection': inter, 'union': union})

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
        
        if 'index' in payload:
            idx = int(payload['index'])
            out = predictor.predict_from_index(idx, k=k, input_types=payload.get('input_types'))
            # predictor may return (scores, indices) or (scores, indices, pred_fp)
            if isinstance(out, tuple) and len(out) >= 2:
                scores, indices = out[0], out[1]
                pred_fp = out[2] if len(out) > 2 else None
                resp = {'topk_scores': scores, 'topk_indices': indices}
                if pred_fp is not None:
                    resp['pred_fp'] = pred_fp
                return jsonify(resp)
            else:
                return jsonify({'error': 'unexpected predictor output'}), 500
        elif 'raw' in payload:
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
                logger.info(f"Returning {len(scores)} scores and {len(indices)} indices")
                resp = {'topk_scores': scores, 'topk_indices': indices}
                if pred_fp is not None:
                    resp['pred_fp'] = pred_fp
                return jsonify(resp)
            else:
                logger.error(f"Unexpected predictor output: {out}")
                return jsonify({'error': 'unexpected predictor output'}), 500
        else:
            return jsonify({'error': "provide 'index' or 'raw' in request body"}), 400
    except Exception as e:
        logger.error(f"Prediction error: {e}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


@app.route('/render-molecule', methods=['POST'])
def render_molecule():
    """Render molecular structure from SMILES string."""
    try:
        payload = request.get_json()
        if not payload or 'smiles' not in payload:
            return jsonify({'error': 'SMILES string required'}), 400
        
        smiles = payload['smiles'].strip()
        if not smiles:
            return jsonify({'error': 'SMILES string cannot be empty'}), 400
        
        if not RDKIT_AVAILABLE:
            return jsonify({'error': 'RDKit not available for molecular rendering'}), 501
        
        # Create molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string'}), 400
        
        # Generate SVG representation
        from rdkit.Chem import Draw
        from rdkit.Chem.Draw import rdMolDraw2D
        
        # Create drawer with smaller dimensions to fit in result cards
        drawer = rdMolDraw2D.MolDraw2DSVG(200, 200)
        drawer.SetFontSize(10)
        
        # Draw molecule
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        
        # Get SVG string
        svg_string = drawer.GetDrawingText()
        
        return jsonify({
            'smiles': smiles,
            'svg': svg_string,
            'molecular_weight': Chem.rdMolDescriptors.CalcExactMolWt(mol),
            'formula': Chem.rdMolDescriptors.CalcMolFormula(mol)
        })
        
    except Exception as e:
        logger.error(f"Render molecule error: {e}", exc_info=True)
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
            
            fp = _fp_loader.build_mfp_for_new_SMILES(smiles)
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
        
        logger.info(f"SMILES search completed: {len(sims_sorted)} results")
        
        return jsonify({
            'topk_scores': sims_sorted.tolist(),
            'topk_indices': idxs.squeeze().tolist(),
            'smiles': smiles
        })
        
    except Exception as e:
        logger.error(f"SMILES search error: {e}", exc_info=True)
        return jsonify({'error': str(e)}), 500