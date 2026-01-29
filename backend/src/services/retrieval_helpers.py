"""
Shared helpers for predict, smiles_search, and secondary_retrieval.
"""
import logging
from typing import Any, List, Optional, Tuple

import torch

from src.config import MOLECULE_IMG_SIZE
from src.domain.drawing.draw import compute_bit_environments_batch
from src.domain.fingerprint.fp_utils import tanimoto_similarity
from src.services.result_builder import build_result_card

logger = logging.getLogger(__name__)


def dense_from_indices(length: int, indices: List[int]) -> torch.Tensor:
    """Build a dense 0/1 tensor from fp column indices."""
    vec = torch.zeros(length, dtype=torch.float32)
    if not indices:
        return vec
    for j in indices:
        if 0 <= j < length:
            vec[j] = 1.0
    return vec


def kept_indices_for_mw_range(
    session: Any,
    mw_min: Optional[float],
    mw_max: Optional[float],
) -> List[int]:
    """
    Return retrieval indices in [mw_min, mw_max]. Uses session's precomputed
    MW index. No filter -> all indices.
    """
    return session.indices_in_mw_range(mw_min, mw_max)


def _retrieved_indices_from_rankingset(
    rankingset: torch.Tensor, global_idx: int
) -> List[int]:
    if rankingset.layout == torch.sparse_csr:
        crow = rankingset.crow_indices()
        col = rankingset.col_indices()
        start = int(crow[global_idx])
        end = int(crow[global_idx + 1])
        if start >= end:
            return []
        return col[start:end].cpu().tolist()
    row = rankingset[global_idx].cpu()
    nz = torch.nonzero(row > 0.5, as_tuple=False).squeeze(-1)
    return nz.tolist() if nz.numel() else []


def _clamp_similarity(val: Any) -> float:
    try:
        s = float(val)
        if s != s:  # NaN
            return 0.0
        return max(0.0, min(1.0, s))
    except (ValueError, TypeError):
        return 0.0


def build_result_cards(
    session: Any,
    rankingset: torch.Tensor,
    pairs: List[Tuple[int, float]],
    pred_tensor: torch.Tensor,
    fp_loader: Any,
    molecule_renderer: Any,
    img_size: int = MOLECULE_IMG_SIZE,
    *,
    max_cards: Optional[int] = None,
    include_plain_svg: bool = True,
) -> List[dict]:
    """
    Build result cards from (global_idx, similarity) pairs. Uses session for
    metadata (get_entry, get_smiles), rankingset for retrieved fp indices.
    """
    cards: List[dict] = []
    for i, (global_idx, sim_val) in enumerate(pairs):
        if max_cards is not None and len(cards) >= max_cards:
            break
        entry = session.get_entry(global_idx)
        if not entry:
            continue
        result_smiles = session.get_smiles(global_idx)
        if not result_smiles:
            continue

        try:
            enhanced_svg = molecule_renderer.render(
                result_smiles,
                predicted_fp=pred_tensor,
                fp_loader=fp_loader,
                img_size=img_size,
            )
        except Exception as e:
            logger.warning("Failed to render enhanced molecule: %s", e)
            enhanced_svg = None

        plain_svg = None
        if include_plain_svg:
            try:
                plain_svg = molecule_renderer.render(
                    result_smiles,
                    predicted_fp=None,
                    fp_loader=None,
                    img_size=img_size,
                )
            except Exception:
                pass

        cosine_sim = _clamp_similarity(sim_val)
        retrieved_indices: List[int] = []
        try:
            retrieved_indices = _retrieved_indices_from_rankingset(
                rankingset, global_idx
            )
        except Exception as e:
            logger.warning("Failed to get retrieved indices for idx %s: %s", global_idx, e)

        tanimoto_sim = 0.0
        try:
            if pred_tensor is not None and retrieved_indices:
                vec = dense_from_indices(pred_tensor.numel(), retrieved_indices)
                tanimoto_sim = tanimoto_similarity(pred_tensor, vec)
        except Exception as e:
            logger.warning("Failed to compute Tanimoto for idx %s: %s", global_idx, e)

        card = build_result_card(
            global_idx,
            entry,
            cosine_sim,
            enhanced_svg,
            cosine_similarity=cosine_sim,
            tanimoto_similarity=tanimoto_sim,
        )
        if plain_svg:
            card["plain_svg"] = plain_svg
        card["retrieved_molecule_fp_indices"] = retrieved_indices

        try:
            if retrieved_indices and fp_loader:
                card["bit_environments"] = compute_bit_environments_batch(
                    result_smiles, retrieved_indices, fp_loader
                )
            else:
                card["bit_environments"] = {}
        except Exception as e:
            logger.warning("Failed to compute bit environments: %s", e)
            card["bit_environments"] = {}

        cards.append(card)
    return cards
