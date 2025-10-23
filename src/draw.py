import io
import base64
from typing import Iterable, Optional, Sequence, Dict, List
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor
from rdkit.Chem.Draw import SimilarityMaps, rdMolDraw2D
from rdkit.Geometry import Point2D
from src.fp_loader import EntropyFPLoader
from src.fp_utils import get_bitinfos
import torch

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
    base_similarity = compute_cos_sim(retrieval_fp, predicted_fp.cpu())

    weights = []
    for atom_id in range(retrieval_mol.GetNumAtoms()):
        masked_fp = fp_loader.build_mfp_from_bitinfo(atom_to_bit_infos, [atom_id])
        masked_sim = compute_cos_sim(predicted_fp.cpu(), masked_fp)
        weights.append(base_similarity - masked_sim)

    weights, _ = SimilarityMaps.GetStandardizedWeights(weights)

    # ---- Draw similarity map ----
    drawer = Draw.MolDraw2DCairo(img_size, img_size)
    draw_high_res_similarity_map(
        retrieval_mol, weights, draw2d=drawer, contourLines=0, gridResolution=0.06
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
    
    # Debug: Check if the bit IDs from get_bitinfos match the fingerprint indices
    print(f"Fingerprint length: {len(fp_diff)}")
    print(f"Fingerprint indices range: 0 to {len(fp_diff)-1}")
    if atom_to_bit_infos:
        sample_atom_bits = list(atom_to_bit_infos.values())[0][:3] if list(atom_to_bit_infos.values())[0] else []
        print(f"Sample bit IDs from get_bitinfos: {[bit_info[0] for bit_info in sample_atom_bits]}")
        max_bit_id = max([bit_info[0] for atom_bits in atom_to_bit_infos.values() for bit_info in atom_bits]) if any(atom_to_bit_infos.values()) else -1
        print(f"Max bit ID from get_bitinfos: {max_bit_id}")
        print(f"Are bit IDs within fingerprint range? {max_bit_id < len(fp_diff) if max_bit_id >= 0 else 'No bits found'}")
    
    # Debug logging
    print(f"=== CHANGE HIGHLIGHTING DEBUG ===")
    print(f"Original FP shape: {original_fp.shape}, sum: {original_fp.sum()}")
    print(f"New FP shape: {new_fp.shape}, sum: {new_fp.sum()}")
    print(f"Original binary sum: {original_binary.sum()}")
    print(f"New binary sum: {new_binary.sum()}")
    print(f"FP diff sum: {fp_diff.sum()}")
    print(f"FP diff non-zero count: {(fp_diff != 0).sum()}")
    print(f"Retrieval SMILES: {retrieval_smiles}")
    print(f"Max radius: {fp_loader.max_radius}")
    print(f"Number of atoms: {retrieval_mol.GetNumAtoms()}")
    print(f"Atom-to-bit mappings count: {len(atom_to_bit_infos)}")
    
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
    
    # Debug: Check if any atoms have associated bits
    print(f"Atom-to-bit mappings sample: {dict(list(atom_to_bit_infos.items())[:3])}")
    print(f"Total bits mapped to atoms: {sum(len(bits) for bits in atom_to_bit_infos.values())}")
    print(f"Sample changed bit indices: {[i for i, val in enumerate(fp_diff) if val != 0][:10]}")
    
    print(f"Change weights: {change_weights}")
    print(f"Max change weight: {max(change_weights) if change_weights else 0}")
    print(f"Non-zero change weights: {[w for w in change_weights if w > 0]}")
    print(f"=== END CHANGE HIGHLIGHTING DEBUG ===")
    
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
    
    # Debug logging
    print(f"=== FINGERPRINT DIFFERENCES DEBUG ===")
    print(f"Original FP shape: {original_fp.shape}, sum: {original_fp.sum()}")
    print(f"New FP shape: {new_fp.shape}, sum: {new_fp.sum()}")
    print(f"Original binary sum: {original_binary.sum()}")
    print(f"New binary sum: {new_binary.sum()}")
    print(f"FP diff sum: {fp_diff.sum()}")
    print(f"FP diff non-zero count: {(fp_diff != 0).sum()}")
    print(f"Retrieval SMILES: {retrieval_smiles}")
    print(f"Max radius: {fp_loader.max_radius}")
    print(f"Number of atoms: {retrieval_mol.GetNumAtoms()}")
    print(f"Atom-to-bit mappings count: {len(atom_to_bit_infos)}")
    
    # Categorize changes
    added_bits = []
    removed_bits = []
    unchanged_bits = []
    
    # Get all unique bit infos from the molecule
    all_bit_infos = set()
    for atom_idx, bit_infos in atom_to_bit_infos.items():
        for bit_info in bit_infos:
            all_bit_infos.add(bit_info)
    
    # Debug: Check the all_bit_infos
    print(f"Total bit infos from molecule: {len(all_bit_infos)}")
    print(f"Sample bit infos: {list(all_bit_infos)[:3]}")
    
    # Check if we need to map bit IDs through the fingerprint loader's mapping
    print(f"FP loader bitinfo_to_fp_index_map available: {hasattr(fp_loader, 'bitinfo_to_fp_index_map')}")
    if hasattr(fp_loader, 'bitinfo_to_fp_index_map'):
        print(f"FP loader mapping size: {len(fp_loader.bitinfo_to_fp_index_map)}")
        print(f"Sample FP loader mapping: {dict(list(fp_loader.bitinfo_to_fp_index_map.items())[:3])}")
    
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
    
    print(f"Added bits count: {len(added_bits)}")
    print(f"Removed bits count: {len(removed_bits)}")
    print(f"Unchanged bits count: {len(unchanged_bits)}")
    if added_bits:
        print(f"Sample added bits: {[bit['bit_id'] for bit in added_bits[:5]]}")
    if removed_bits:
        print(f"Sample removed bits: {[bit['bit_id'] for bit in removed_bits[:5]]}")
    print(f"=== END FINGERPRINT DIFFERENCES DEBUG ===")
    
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
