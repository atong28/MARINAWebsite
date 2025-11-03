import threading
from typing import Optional
from functools import lru_cache
from src.config import MOLECULE_IMG_SIZE, ENABLE_RENDERING_CACHE, RENDERING_CACHE_SIZE

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    RDKIT_AVAILABLE = True
except Exception:
    RDKIT_AVAILABLE = False

from src.draw import draw_molecule
import base64
import io


class MoleculeRenderer:
    _instance = None
    _lock = threading.Lock()

    @classmethod
    def instance(cls) -> "MoleculeRenderer":
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = MoleculeRenderer()
        return cls._instance

    def render(self, smiles: str, predicted_fp=None, fp_loader=None, img_size: int = MOLECULE_IMG_SIZE) -> Optional[str]:
        if ENABLE_RENDERING_CACHE:
            return _render_cached(smiles, predicted_fp, fp_loader, img_size)
        return _render_impl(smiles, predicted_fp, fp_loader, img_size)


def _to_key(smiles: str, predicted_fp, img_size: int) -> str:
    fp_key = None
    try:
        if predicted_fp is not None:
            # Use a cheap hash (length + sum) to avoid large keys
            fp_key = f"{len(predicted_fp)}:{round(float(sum(predicted_fp)), 3)}"
    except Exception:
        fp_key = "x"
    return f"{smiles}|{img_size}|{fp_key}"


def _render_impl(smiles: str, predicted_fp=None, fp_loader=None, img_size: int = MOLECULE_IMG_SIZE) -> Optional[str]:
    try:
        if RDKIT_AVAILABLE and predicted_fp is not None and fp_loader is not None:
            enhanced_img = draw_molecule(
                predicted_fp=predicted_fp,
                retrieval_smiles=smiles,
                fp_loader=fp_loader,
                need_to_clean_H=False,
                img_size=img_size
            )
            buffer = io.BytesIO()
            enhanced_img.save(buffer, format='PNG')
            img_str = base64.b64encode(buffer.getvalue()).decode()
            return f"data:image/png;base64,{img_str}"
        elif RDKIT_AVAILABLE:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return Draw.MolToSVG(mol, width=img_size//2, height=img_size//2)
            return None
        else:
            return None
    except Exception:
        return None


@lru_cache(maxsize=RENDERING_CACHE_SIZE)
def _render_cached(smiles: str, predicted_fp_tuple_hashable, fp_loader_id, img_size: int) -> Optional[str]:
    # predicted_fp_tuple_hashable is either None or a tuple (len,sum) used only for keying;
    # we still need to pass the real fp to renderer; cache benefits mostly from repeated smiles.
    # Here we only render basic SVG for cache to avoid fp mismatch. For enhanced, skip cache.
    try:
        if RDKIT_AVAILABLE:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return Draw.MolToSVG(mol, width=img_size//2, height=img_size//2)
    except Exception:
        pass
    return None

def clear_cache():
    try:
        _render_cached.cache_clear()
    except Exception:
        pass


