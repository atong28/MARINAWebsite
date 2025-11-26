from typing import Dict, Optional
import logging
from rdkit import Chem
from rdkit.Chem import Descriptors

from src.domain.models.prediction_result import DatabaseLinks

logger = logging.getLogger(__name__)

def build_result_card(
    idx: int, 
    entry: Dict, 
    similarity: float, 
    svg: Optional[str],
    *,
    cosine_similarity: Optional[float] = None,
    tanimoto_similarity: Optional[float] = None,
) -> Dict:
    """Build a result card dictionary from metadata entry."""
    # Choose display SMILES
    logger.info(f"Building result card for index: {idx}")
    logger.info(f"Entry: {entry}")
    smiles = entry.get('canonical_3d_smiles')
    if smiles is None or smiles == 'N/A':
        smiles = entry.get('canonical_2d_smiles')
    
    # Compute exact mass using RDKit
    exact_mass = None
    if smiles:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                exact_mass = Descriptors.ExactMolWt(mol)
        except Exception as e:
            pass
    if exact_mass is None:
        logger.error(f"Exact mass is None for smiles: {smiles}")
    
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

    primary_similarity = similarity if similarity is not None else 0.0
    if cosine_similarity is not None:
        primary_similarity = cosine_similarity

    return {
        'index': idx,
        'smiles': smiles,
        'similarity': float(primary_similarity),
        'cosine_similarity': float(cosine_similarity) if cosine_similarity is not None else float(primary_similarity),
        'tanimoto_similarity': float(tanimoto_similarity) if tanimoto_similarity is not None else None,
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

