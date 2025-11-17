from typing import Dict, Optional

from rdkit import Chem
from rdkit.Chem import Descriptors

from src.domain.models.prediction_result import DatabaseLinks


def build_result_card(
    idx: int, 
    entry: Dict, 
    similarity: float, 
    svg: Optional[str]
) -> Dict:
    """Build a result card dictionary from metadata entry."""
    # Choose display SMILES
    smiles = entry.get('canonical_3d_smiles') if entry.get('canonical_3d_smiles') != 'N/A' else entry.get('canonical_2d_smiles')
    
    # Compute exact mass using RDKit
    exact_mass = None
    if smiles:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                exact_mass = Descriptors.ExactMolWt(mol)
        except Exception as e:
            # Handle errors gracefully - exact_mass remains None
            print(f"Error computing exact mass for smiles: {smiles}")
            print(f"Error: {e}")
            pass
    
    # Primary name/link
    primary_name = None
    primary_link = None
    coconut = entry.get('coconut')
    lotus = entry.get('lotus')
    npmrd = entry.get('npmrd')
    
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

    # Database links
    database_links = DatabaseLinks()
    if coconut and coconut.get('coconut_id'):
        database_links.coconut = f"https://coconut.naturalproducts.net/compounds/{coconut['coconut_id']}"
    if lotus and lotus.get('lotus_id'):
        database_links.lotus = f"https://lotus.naturalproducts.net/compounds/{lotus['lotus_id']}"
    if npmrd and npmrd.get('npmrd_id'):
        database_links.npmrd = f"https://np-mrd.org/natural_products/{npmrd['npmrd_id']}"

    return {
        'index': idx,
        'smiles': smiles,
        'similarity': float(similarity) if similarity is not None else 0.0,
        'svg': svg,
        'name': primary_name,
        'primary_link': primary_link,
        'database_links': database_links.dict(),
        'np_pathway': None,
        'np_superclass': None,
        'np_class': None,
        'retrieved_molecule_fp_indices': [],
        'bit_environments': {},
        'exact_mass': exact_mass
    }

