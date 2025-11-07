import threading
from typing import Optional, Dict
import base64
import io

from rdkit import Chem
from rdkit.Chem import Draw

from src.domain.drawing.draw import draw_molecule, render_molecule_with_overlays
from src.config import MOLECULE_IMG_SIZE, HIGHLIGHT_DEBUG

import torch


class MoleculeRenderer:
    """Service for rendering molecules. No caching - render on demand."""
    _instance: Optional['MoleculeRenderer'] = None
    _lock = threading.Lock()

    @classmethod
    def instance(cls) -> "MoleculeRenderer":
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = MoleculeRenderer()
        return cls._instance

    def render(
        self, 
        smiles: str, 
        predicted_fp=None, 
        fp_loader=None, 
        img_size: int = MOLECULE_IMG_SIZE,
        bit_environments: Optional[Dict[int, dict]] = None
    ) -> Optional[str]:
        """
        Render a molecule. If bit_environments are provided, render SVG with embedded overlays.
        If predicted_fp and fp_loader are provided, render enhanced PNG.
        Otherwise render basic SVG.
        """
        # If bit environments are provided, render SVG with embedded overlays
        if bit_environments is not None:
            try:
                svg = render_molecule_with_overlays(smiles, bit_environments, img_size)
                return svg
            except Exception as e:
                if HIGHLIGHT_DEBUG:
                    import logging
                    logging.getLogger(__name__).warning(f"Overlay rendering failed: {e}")
                # Fall through to other rendering methods
        
        # If we can enhance (have fingerprint + loader), render enhanced
        if predicted_fp is not None and fp_loader is not None:
            try:
                img = self._render_enhanced(smiles, predicted_fp, fp_loader, img_size)
                return img
            except Exception as e:
                # Fall back to basic rendering
                if HIGHLIGHT_DEBUG:
                    import logging
                    logging.getLogger(__name__).warning(f"Enhanced rendering failed: {e}")
        
        # Basic rendering
        try:
            img = self._render_basic(smiles, img_size)
            return img
        except Exception as e:
            if HIGHLIGHT_DEBUG:
                import logging
                logging.getLogger(__name__).warning(f"Basic rendering failed: {e}")
            return None

    def _render_enhanced(
        self, 
        smiles: str, 
        predicted_fp, 
        fp_loader, 
        img_size: int
    ) -> Optional[str]:
        """Render molecule with fingerprint highlighting."""
        
        if isinstance(predicted_fp, list):
            predicted_fp = torch.tensor(predicted_fp, dtype=torch.float32)
        
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

    def _render_basic(self, smiles: str, img_size: int) -> Optional[str]:
        """Render basic molecule SVG without highlighting."""
        
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Draw.MolToSVG(mol, width=img_size//2, height=img_size//2)
        return None

