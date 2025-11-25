import io
import base64
from typing import Iterable, Optional, Sequence, Dict, List
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor
from rdkit.Chem.Draw import SimilarityMaps, rdMolDraw2D
from rdkit.Geometry import Point2D
from src.domain.fingerprint.fp_loader import EntropyFPLoader
from src.domain.fingerprint.fp_utils import get_bitinfos, _mk_rdkit
import torch
from functools import lru_cache
import xml.etree.ElementTree as ET

def compute_cos_sim(fp1: torch.Tensor, fp2: torch.Tensor, eps: float = 1e-12) -> float:
    """
    Robust cosine similarity between 1D tensors.
    Returns 0.0 if either vector has (near-)zero norm.
    """
    fp1 = fp1.flatten().to(dtype=torch.float32, copy=False)
    fp2 = fp2.flatten().to(dtype=torch.float32, copy=False)

    n1 = torch.linalg.norm(fp1).item()
    n2 = torch.linalg.norm(fp2).item()
    if n1 < eps or n2 < eps:
        return 0.0

    dot = torch.dot(fp1, fp2).item()
    return float(dot / (n1 * n2))

def draw_high_res_similarity_map(
    mol: Chem.Mol,
    weights: Sequence[float],
    draw2d: rdMolDraw2D.MolDraw2D,
    colorMap: Optional[Iterable] = None,
    sigma: Optional[float] = None,
    contourLines: int = 10,
    gridResolution: float = 0.05,
) -> rdMolDraw2D.MolDraw2D:
    """
    High-res similarity map (cleaned version of RDKit's SimilarityMaps.GetSimilarityMapFromWeights).
    - Only requires an external MolDraw2D (so you control size/format).
    - Keeps background; overlays Gaussian contours colored by `weights`.
    """
    if mol is None or mol.GetNumAtoms() < 2:
        raise ValueError("Molecule must have at least 2 atoms.")
    if draw2d is None:
        raise ValueError("draw2d must be provided.")

    # Prepare 2D coords if missing
    mprep = rdMolDraw2D.PrepareMolForDrawing(mol, addChiralHs=False)
    if not mprep.GetNumConformers():
        rdDepictor.Compute2DCoords(mprep)

    # Heuristic sigma based on first bond length (in 2D)
    if sigma is None:
        conf = mprep.GetConformer()
        if mprep.GetNumBonds() > 0:
            b = mprep.GetBondWithIdx(0)
            p1 = conf.GetAtomPosition(b.GetBeginAtomIdx())
            p2 = conf.GetAtomPosition(b.GetEndAtomIdx())
        else:
            p1 = conf.GetAtomPosition(0)
            p2 = conf.GetAtomPosition(1)
        sigma = round(0.3 * (p1 - p2).Length(), 2)

    # Atom locations
    conf = mprep.GetConformer()
    locs = [Point2D(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y) for i in range(mprep.GetNumAtoms())]
    sigmas = [sigma] * mprep.GetNumAtoms()

    # Configure contour params
    draw2d.ClearDrawing()
    params = Draw.ContourParams()
    params.fillGrid = True
    params.gridResolution = gridResolution
    params.extraGridPadding = 0.5

    # Optional custom colormap
    if colorMap is not None:
        try:
            # try matplotlib-like API
            import matplotlib
            if isinstance(colorMap, str):
                cm = matplotlib.colormaps[colorMap]
                clrs = [tuple(x) for x in cm([0.0, 0.5, 1.0])]
            else:
                # assume iterable of 3 RGB tuples
                clrs = list(colorMap)
            params.setColourMap(clrs)
        except Exception:
            # fall back silently to RDKit default if matplotlib not available / invalid spec
            pass

    # Draw contours + molecule
    Draw.ContourAndDrawGaussians(draw2d, locs, weights, sigmas, nContours=contourLines, params=params)
    draw2d.drawOptions().clearBackground = False
    draw2d.DrawMolecule(mprep)
    return draw2d


def draw_molecule(
    predicted_fp: torch.Tensor,
    retrieval_smiles: str,
    fp_loader: EntropyFPLoader,
    need_to_clean_H: bool = False,
    img_size: int = 1200,  # 3x larger resolution (400 * 3)
):
    """
    Visualize a retrieved molecule with atom-level similarity highlighting
    based on the predicted fingerprint vector.

    Args:
        predicted_fp: torch.Tensor — predicted fingerprint (1D float tensor)
        retrieval_smiles: str — SMILES of the retrieved molecule
        need_to_clean_H: bool — whether to remove explicit hydrogens
        img_size: int — output image size in pixels

    Returns:
        PIL.Image.Image — image with highlighted similarity map
    """
    def _show_png(data: bytes) -> Image.Image:
        return Image.open(io.BytesIO(data))

    # ---- Initialize loader ----
    retrieval_mol = Chem.MolFromSmiles(retrieval_smiles)
    if retrieval_mol is None:
        raise ValueError(f"Invalid SMILES: {retrieval_smiles}")

    if need_to_clean_H:
        retrieval_mol = Chem.RemoveHs(retrieval_mol)

    # ---- Compute per-atom contributions ----
    atom_to_bit_infos, _ = get_bitinfos(retrieval_smiles, fp_loader.max_radius)
    retrieval_fp = fp_loader.build_mfp_from_bitinfo(atom_to_bit_infos)
    pred = predicted_fp.detach().to(dtype=torch.float32).cpu().flatten()
    base_similarity = compute_cos_sim(retrieval_fp, pred)

    weights = []
    for atom_id in range(retrieval_mol.GetNumAtoms()):
        masked_fp = fp_loader.build_mfp_from_bitinfo(atom_to_bit_infos, [atom_id])
        masked_sim = compute_cos_sim(pred, masked_fp)
        weights.append(base_similarity - masked_sim)

    # Standardize weights and add robustness: if near-flat, fall back to boolean contribution
    weights, _ = SimilarityMaps.GetStandardizedWeights(weights)
    try:
        import numpy as _np
        if _np.std(weights) < 1e-6:
            bool_weights = [1.0 if atom_to_bit_infos.get(aid) else 0.0 for aid in range(retrieval_mol.GetNumAtoms())]
            weights = bool_weights
    except Exception:
        pass

    # ---- Draw similarity map ----
    drawer = Draw.MolDraw2DCairo(img_size, img_size)
    # Increase contour lines and grid density for stronger visual contrast
    draw_high_res_similarity_map(
        retrieval_mol, weights, draw2d=drawer, contourLines=5, gridResolution=0.04
    )
    drawer.FinishDrawing()
    img = _show_png(drawer.GetDrawingText())
    return img


def draw_fingerprint_changes(
    original_fp: torch.Tensor,
    new_fp: torch.Tensor,
    retrieval_smiles: str,
    fp_loader: EntropyFPLoader,
    need_to_clean_H: bool = False,
    img_size: int = 1200,  # 3x larger resolution (400 * 3)
) -> Image.Image:
    """
    Visualize fingerprint changes between original and new predictions.
    
    Args:
        original_fp: torch.Tensor — original fingerprint (1D float tensor)
        new_fp: torch.Tensor — new fingerprint (1D float tensor)
        retrieval_smiles: str — SMILES of the molecule to highlight
        fp_loader: EntropyFPLoader — fingerprint loader instance
        need_to_clean_H: bool — whether to remove explicit hydrogens
        img_size: int — output image size in pixels
    
    Returns:
        PIL.Image.Image — image with highlighted changes (red for removed, green for added)
    """
    def _show_png(data: bytes) -> Image.Image:
        return Image.open(io.BytesIO(data))

    # Initialize molecule
    retrieval_mol = Chem.MolFromSmiles(retrieval_smiles)
    if retrieval_mol is None:
        raise ValueError(f"Invalid SMILES: {retrieval_smiles}")

    if need_to_clean_H:
        retrieval_mol = Chem.RemoveHs(retrieval_mol)

    # Get atom-to-bit mappings
    atom_to_bit_infos, _ = get_bitinfos(retrieval_smiles, fp_loader.max_radius)
    
    # Compute fingerprint changes - convert to binary for proper comparison
    original_binary = (original_fp.cpu() > 0.5).float()
    new_binary = (new_fp.cpu() > 0.5).float()
    fp_diff = new_binary - original_binary  # -1: removed, 0: no change, 1: added
    
    # For each atom, compute its contribution to the fingerprint changes
    num_atoms = retrieval_mol.GetNumAtoms()
    change_weights = []
    
    for atom_id in range(num_atoms):
        # Compute the contribution of this atom to the fingerprint change
        atom_contribution = 0.0
        for bit_info in atom_to_bit_infos.get(atom_id, []):
            bit_idx = bit_info[0]
            
            # Check if we need to map through the fingerprint loader's mapping
            if hasattr(fp_loader, 'bitinfo_to_fp_index_map'):
                fp_index = fp_loader.bitinfo_to_fp_index_map.get(bit_info)
                if fp_index is None:
                    continue  # This bit is not in the selected fingerprint
                bit_idx = fp_index
            
            if bit_idx < len(fp_diff):
                # Weight positive changes (added bits) more heavily
                if fp_diff[bit_idx].item() > 0:  # Added bit
                    atom_contribution += 2.0  # Green highlighting
                elif fp_diff[bit_idx].item() < 0:  # Removed bit
                    atom_contribution += 1.0  # Red highlighting
        
        change_weights.append(atom_contribution)
    
    # Normalize weights to [0, 1] range
    if change_weights:
        max_weight = max(change_weights)
        if max_weight > 0:
            change_weights = [w / max_weight for w in change_weights]
    
    # Ensure we have some weights to work with
    if not change_weights:
        change_weights = [0.0] * num_atoms
    
    # Create custom colormap: white to red for changes
    custom_colors = [(1.0, 1.0, 1.0), (1.0, 0.8, 0.8), (1.0, 0.0, 0.0)]  # White to light red to red
    
    # Render the molecule with change highlighting
    drawer = Draw.MolDraw2DCairo(img_size, img_size)
    
    # Use a more targeted approach to avoid uniform overlay
    # Only highlight atoms with significant changes
    significant_atoms = []
    significant_colors = []
    for i, weight in enumerate(change_weights):
        if weight > 0.1:  # Only highlight atoms with significant changes
            significant_atoms.append(i)
            # Color intensity based on weight
            intensity = min(weight, 1.0)
            significant_colors.append((1.0, 1.0 - intensity, 1.0 - intensity))  # White to red gradient
    
    if significant_atoms:
        # Draw molecule with atom highlighting
        # Convert colors to RDKit format (RGBA tuples)
        rdkit_colors = {}
        for i, atom_idx in enumerate(significant_atoms):
            color = significant_colors[i]
            rdkit_colors[atom_idx] = color  # RDKit expects RGB tuple
        
        drawer.DrawMolecule(retrieval_mol, highlightAtoms=significant_atoms, highlightAtomColors=rdkit_colors)
    else:
        # Draw molecule without highlighting if no significant changes
        drawer.DrawMolecule(retrieval_mol)
    
    drawer.FinishDrawing()
    img = _show_png(drawer.GetDrawingText())
    return img


def draw_substructure(smiles: str, img_size: int = 200) -> Optional[Image.Image]:
    """Draw a substructure from SMILES string"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        drawer = Draw.MolDraw2DCairo(img_size, img_size)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        
        img_data = drawer.GetDrawingText()
        img = Image.open(io.BytesIO(img_data))
        return img
    except Exception as e:
        print(f"Error drawing substructure {smiles}: {e}")
        return None


def pil_image_to_base64(img: Image.Image) -> str:
    """Convert PIL Image to base64 string"""
    try:
        buffer = io.BytesIO()
        img.save(buffer, format='PNG')
        img_str = base64.b64encode(buffer.getvalue()).decode()
        return f"data:image/png;base64,{img_str}"
    except Exception as e:
        print(f"Error converting image to base64: {e}")
        return None


def draw_similarity_comparison(
    predicted_fp: torch.Tensor,
    retrieval_smiles: str,
    fp_loader: EntropyFPLoader,
    need_to_clean_H: bool = False,
    img_size: int = 1200,  # 3x larger resolution (400 * 3)
) -> Image.Image:
    """
    Visualize similarity highlighting for the analysis page.
    This is a wrapper around the existing draw_molecule function.
    
    Args:
        predicted_fp: torch.Tensor — predicted fingerprint (1D float tensor)
        retrieval_smiles: str — SMILES of the retrieved molecule
        fp_loader: EntropyFPLoader — fingerprint loader instance
        need_to_clean_H: bool — whether to remove explicit hydrogens
        img_size: int — output image size in pixels
    
    Returns:
        PIL.Image.Image — image with similarity highlighting
    """
    return draw_molecule(
        predicted_fp=predicted_fp,
        retrieval_smiles=retrieval_smiles,
        fp_loader=fp_loader,
        need_to_clean_H=need_to_clean_H,
        img_size=img_size
    )


def get_fingerprint_differences(
    original_fp: torch.Tensor,
    new_fp: torch.Tensor,
    retrieval_smiles: str,
    fp_loader: EntropyFPLoader,
    need_to_clean_H: bool = False,
) -> Dict[str, List[Dict]]:
    """
    Get detailed fingerprint differences with associated substructures.
    
    Args:
        original_fp: torch.Tensor — original fingerprint
        new_fp: torch.Tensor — new fingerprint
        retrieval_smiles: str — SMILES of the molecule
        fp_loader: EntropyFPLoader — fingerprint loader instance
        need_to_clean_H: bool — whether to remove explicit hydrogens
    
    Returns:
        Dict with 'added', 'removed', and 'unchanged' lists containing bit info
    """
    # Initialize molecule
    retrieval_mol = Chem.MolFromSmiles(retrieval_smiles)
    if retrieval_mol is None:
        return {'added': [], 'removed': [], 'unchanged': []}

    if need_to_clean_H:
        retrieval_mol = Chem.RemoveHs(retrieval_mol)

    # Get atom-to-bit mappings
    atom_to_bit_infos, _ = get_bitinfos(retrieval_smiles, fp_loader.max_radius)
    
    # Compute fingerprint changes - convert to binary for proper comparison
    original_binary = (original_fp.cpu() > 0.5).float()
    new_binary = (new_fp.cpu() > 0.5).float()
    fp_diff = new_binary - original_binary  # -1: removed, 0: no change, 1: added
    
    # Categorize changes
    added_bits = []
    removed_bits = []
    unchanged_bits = []
    
    # Get all unique bit infos from the molecule
    all_bit_infos = set()
    for atom_idx, bit_infos in atom_to_bit_infos.items():
        for bit_info in bit_infos:
            all_bit_infos.add(bit_info)
    
    # Check each bit info against the fingerprint differences
    for bit_info in all_bit_infos:
        try:
            bit_id = bit_info[0]
            
            # Check if we need to map through the fingerprint loader's mapping
            if hasattr(fp_loader, 'bitinfo_to_fp_index_map'):
                fp_index = fp_loader.bitinfo_to_fp_index_map.get(bit_info)
                if fp_index is None:
                    continue  # This bit is not in the selected fingerprint
                bit_id = fp_index
            
            if bit_id < len(fp_diff):
                change = fp_diff[bit_id].item()
                
                bit_data = {
                    'bit_id': bit_id,
                    'atom_symbol': bit_info[1],
                    'fragment_smiles': bit_info[2],
                    'radius': bit_info[3],
                    'change': change
                }
                
                if change > 0:  # Added bit
                    added_bits.append(bit_data)
                elif change < 0:  # Removed bit
                    removed_bits.append(bit_data)
                else:  # No change
                    unchanged_bits.append(bit_data)
        except Exception as e:
            print(f"Error processing bit_info {bit_info}: {e}")
            continue
    
    # Add substructure images to the bit data
    for bit_data in added_bits + removed_bits + unchanged_bits:
        try:
            # Try to render the substructure
            substructure_img = draw_substructure(bit_data['fragment_smiles'], img_size=200)
            if substructure_img:
                bit_data['substructure_image'] = pil_image_to_base64(substructure_img)
        except Exception as e:
            print(f"Error rendering substructure {bit_data['fragment_smiles']}: {e}")
            bit_data['substructure_image'] = None
    
    return {
        'added': added_bits,
        'removed': removed_bits,
        'unchanged': unchanged_bits
    }


@lru_cache(maxsize=2000)
def render_bit_preview_base64(smiles: str, bit_index: int, radius: int, img_size: int = 200) -> Optional[str]:
    """
    Render a small substructure preview image for a given fingerprint bit index.
    Attempts to locate a BitInfo whose mapped fingerprint index equals bit_index.
    Falls back to matching by raw bit_id when mapping is identity.
    Returns a base64 data URL string or None.
    """
    try:
        atom_to_bit_infos, all_bit_infos = get_bitinfos(smiles, radius)
        if not all_bit_infos:
            return None
        # Try exact raw bit_id match first
        candidates = [bi for bi in all_bit_infos if bi[0] == bit_index]
        # If not found, there may be a mapping in the loader (handled by wrapper below)
        target_frag_smiles = None
        if candidates:
            target_frag_smiles = candidates[0][2]
        else:
            # Heuristic: pick a fragment whose string hash mod big prime equals the index (very rare fallback)
            for bi in all_bit_infos:
                if (abs(hash(bi[2])) % 1000003) == (bit_index % 1000003):
                    target_frag_smiles = bi[2]
                    break
        if not target_frag_smiles:
            return None
        img = draw_substructure(target_frag_smiles, img_size=img_size)
        if not img:
            return None
        return pil_image_to_base64(img)
    except Exception:
        return None


def render_bit_preview(smiles: str, bit_index: int, fp_loader: EntropyFPLoader, img_size: int = 200) -> Optional[str]:
    """
    Wrapper that uses the loader's mapping (if available) to find BitInfo candidates for a fp index.
    Returns base64 data URL or None.
    """
    try:
        # If loader exposes inverse mapping, use it
        inverse = None
        if hasattr(fp_loader, 'bitinfo_to_fp_index_map'):
            try:
                mapping = fp_loader.bitinfo_to_fp_index_map
                inverse = {}
                for bi, idx in mapping.items():
                    inverse.setdefault(idx, []).append(bi)
            except Exception:
                inverse = None
        if inverse and bit_index in inverse and inverse[bit_index]:
            frag = inverse[bit_index][0][2]
            img = draw_substructure(frag, img_size=img_size)
            if img:
                return pil_image_to_base64(img)
        # Fallback to base function using raw bit ids
        return render_bit_preview_base64(smiles, bit_index, getattr(fp_loader, 'max_radius', 2), img_size)
    except Exception:
        return None


def compute_bit_environment(smiles: str, bit_index: int, radius: int):
    """
    Return dict with atom indices, bond index pairs, and normalized atom 2D coordinates for a given bit index.
    Coordinates normalized to [0,1] based on drawing conformer.
    """
    mol = Chem.MolFromSmiles(smiles)
    Chem.RemoveStereochemistry(mol)
    if mol is None:
        return None
    # Ensure 2D coords
    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    # Try to locate environment by bit id using RDKit bit info
    gen, ao = _mk_rdkit(radius)
    _ = gen.GetSparseFingerprint(mol, additionalOutput=ao)
    info = ao.GetBitInfoMap()
    atom_envs = info.get(bit_index)
    if not atom_envs:
        return None
    # Use first occurrence
    atom_idx, curr_radius = atom_envs[0]
    env_bonds = Chem.FindAtomEnvironmentOfRadiusN(mol, curr_radius, atom_idx)
    atom_set = set()
    for bidx in env_bonds:
        b = mol.GetBondWithIdx(bidx)
        atom_set.add(b.GetBeginAtomIdx())
        atom_set.add(b.GetEndAtomIdx())
    atom_list = sorted(atom_set) if atom_set else [atom_idx]
    # Compute normalized coords
    xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
    minx, maxx = min(xs), max(xs)
    miny, maxy = min(ys), max(ys)
    spanx = max(maxx - minx, 1e-6)
    spany = max(maxy - miny, 1e-6)
    atom_coords = []
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        atom_coords.append({
            'id': i,
            'x': (p.x - minx) / spanx,
            'y': (p.y - miny) / spany
        })
    bond_pairs = []
    for bidx in env_bonds:
        b = mol.GetBondWithIdx(bidx)
        bond_pairs.append([b.GetBeginAtomIdx(), b.GetEndAtomIdx()])
    return {
        'atoms': atom_list,
        'bonds': bond_pairs,
        'coords': {'atoms': atom_coords}
    }


def render_molecule_with_overlays(
    smiles: str,
    bit_environments: Dict[int, dict],
    img_size: int = 400,
    bit_colors: Optional[Dict[int, Dict[str, str]]] = None,
    default_color: str = 'rgba(16,185,129,0.7)',
    default_fill_color: str = 'rgba(16,185,129,0.25)'
) -> str:
    """
    Render molecule SVG with embedded bit environment overlays.
    
    Args:
        smiles: SMILES string of the molecule
        bit_environments: Dict mapping bit_index -> environment dict with atoms, bonds, coords
        img_size: Size of the SVG in pixels
    
    Returns:
        SVG string with embedded overlay elements (initially hidden)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(img_size, img_size)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    base_svg = drawer.GetDrawingText()
    atom_positions = {}
    for i in range(mol.GetNumAtoms()):
        pt = drawer.GetDrawCoords(i)  # RDGeom.Point2D in canvas coords
        atom_positions[i] = (pt.x, pt.y)

    root = ET.fromstring(base_svg)
    svg_elem = root if root.tag.endswith('svg') else root.find('.//{http://www.w3.org/2000/svg}svg')
    if svg_elem is None:
        svg_elem = root if root.tag == 'svg' else root.find('.//svg')
    
    if svg_elem is None:
        raise ValueError("Could not find SVG element in generated SVG")
    
    overlays_group = ET.Element('g', {
        'id': 'bit-overlays'
    })
    
    color_map = bit_colors or {}

    for bit_index, env in bit_environments.items():
        if not env or 'atoms' not in env or 'bonds' not in env:
            continue
        
        colors = color_map.get(bit_index, {})
        stroke_color = colors.get('stroke', default_color)
        fill_color = colors.get('fill', default_fill_color)
        
        bit_group = ET.Element('g', {
            'id': f'bit-overlay-{bit_index}',
            'class': 'bit-overlay',
            'data-bit-index': str(bit_index),
            'style': 'display:none'  # Initially hidden
        })
        
        if 'bonds' in env and isinstance(env['bonds'], list):
            for bond_idx, (atom1, atom2) in enumerate(env['bonds']):
                if atom1 not in atom_positions or atom2 not in atom_positions:
                    continue
                pos1 = atom_positions[atom1]
                pos2 = atom_positions[atom2]
                line = ET.Element('line', {
                    'id': f'bit-{bit_index}-bond-{bond_idx}',
                    'class': 'bit-bond',
                    'x1': str(pos1[0]),
                    'y1': str(pos1[1]),
                    'x2': str(pos2[0]),
                    'y2': str(pos2[1]),
                    'stroke': stroke_color,
                    'stroke-width': '4',
                    'stroke-linecap': 'round'
                })
                bit_group.append(line)
        
        if 'atoms' in env and isinstance(env['atoms'], list):
            for atom_idx, atom_id in enumerate(env['atoms']):
                if atom_id not in atom_positions:
                    continue
                pos = atom_positions[atom_id]
                circle = ET.Element('circle', {
                    'id': f'bit-{bit_index}-atom-{atom_idx}',
                    'class': 'bit-atom',
                    'cx': str(pos[0]),
                    'cy': str(pos[1]),
                    'r': '10',
                    'fill': fill_color,
                    'stroke': stroke_color,
                    'stroke-width': '2'
                })
                bit_group.append(circle)
        
        overlays_group.append(bit_group)
    
    svg_elem.append(overlays_group)
    ET.register_namespace('', 'http://www.w3.org/2000/svg')
    return ET.tostring(root, encoding='unicode')


def render_molecule_with_change_overlays(
    smiles: str,
    original_fp: Sequence[float],
    new_fp: Sequence[float],
    fp_loader,
    threshold: float = 0.5,
    img_size: int = 400
) -> str:
    """Render molecule overlays highlighting gained (green) and lost (red) fingerprint bits."""

    if fp_loader is None:
        raise ValueError("Fingerprint loader is required to render change overlays")

    original_bits = {idx for idx, value in enumerate(original_fp) if value is not None and value > threshold}
    new_bits = {idx for idx, value in enumerate(new_fp) if value is not None and value > threshold}

    gained_bits = sorted(new_bits - original_bits)
    lost_bits = sorted(original_bits - new_bits)

    if not gained_bits and not lost_bits:
        # No changes; return empty overlays
        return render_molecule_with_overlays(smiles, {}, img_size=img_size)

    target_bits = list({*gained_bits, *lost_bits})

    bit_environments = compute_bit_environments_batch(smiles, target_bits, fp_loader) if target_bits else {}

    gain_color = {
        'stroke': 'rgba(34,197,94,0.7)',
        'fill': 'rgba(34,197,94,0.25)'
    }
    loss_color = {
        'stroke': 'rgba(239,68,68,0.7)',
        'fill': 'rgba(239,68,68,0.25)'
    }

    bit_colors: Dict[int, Dict[str, str]] = {}
    for idx in gained_bits:
        bit_colors[idx] = gain_color
    for idx in lost_bits:
        bit_colors[idx] = loss_color

    return render_molecule_with_overlays(
        smiles,
        bit_environments,
        img_size=img_size,
        bit_colors=bit_colors,
    )


def compute_bit_environments_batch(smiles: str, fp_indices: List[int], fp_loader) -> Dict[int, dict]:
    """
    Compute bit environments for multiple fingerprint indices in batch.
    Returns dict mapping fp_index -> env_data, only including bits with valid environments (omits null).
    
    Args:
        smiles: SMILES string of the molecule
        fp_indices: List of fingerprint indices to compute environments for
        fp_loader: EntropyFPLoader instance with mapping information
        
    Returns:
        Dict[int, dict]: Mapping from fp_index to environment data (only valid envs included)
    """
    if not fp_indices or not smiles:
        return {}
    
    result = {}
    max_radius = getattr(fp_loader, 'max_radius', 2)
    
    # Build mapping: fp_index -> (raw_bit_id, radius) for fast lookup
    fp_to_raw = {}
    
    # Strategy 1: Use fp_index_to_bitinfo_map if available (fastest)
    if hasattr(fp_loader, 'fp_index_to_bitinfo_map') and fp_loader.fp_index_to_bitinfo_map:
        for fp_idx in fp_indices:
            bi = fp_loader.fp_index_to_bitinfo_map.get(int(fp_idx))
            if bi is not None and isinstance(bi, tuple) and len(bi) >= 4:
                raw_bit_id = int(bi[0])
                r = int(bi[3])
                fp_to_raw[fp_idx] = (raw_bit_id, r)
    
    # Strategy 2: For unmapped indices, derive from SMILES bitinfos
    unmapped = [idx for idx in fp_indices if idx not in fp_to_raw]
    if unmapped and hasattr(fp_loader, 'bitinfo_to_fp_index_map') and fp_loader.bitinfo_to_fp_index_map:
        atom_to_bitinfos, all_bit_infos = get_bitinfos(smiles, max_radius)
        for cand in (all_bit_infos or []):
            mapped_idx = fp_loader.bitinfo_to_fp_index_map.get(cand)
            if mapped_idx is not None and mapped_idx in unmapped:
                raw_bit_id = int(cand[0])
                r = int(cand[3])
                fp_to_raw[mapped_idx] = (raw_bit_id, r)
                unmapped.remove(mapped_idx)
    
    # Strategy 3: Fallback - treat remaining indices as raw bit IDs
    for fp_idx in unmapped:
        fp_to_raw[fp_idx] = (int(fp_idx), max_radius)
    
    # Compute environments for all mapped bits
    for fp_idx, (raw_bit_id, radius) in fp_to_raw.items():
        env = compute_bit_environment(smiles, raw_bit_id, radius)
        if env is not None:
            result[fp_idx] = env
        else:
            # Try radius sweep as fallback
            for r in range(max_radius, -1, -1):
                env = compute_bit_environment(smiles, raw_bit_id, r)
                if env is not None:
                    result[fp_idx] = env
                    break
    
    return result
