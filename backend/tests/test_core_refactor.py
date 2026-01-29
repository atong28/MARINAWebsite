import pytest
import torch

from src.domain.ranker import RankingSet
from src.domain.predictor import _preprocess_raw_inputs


def test_ranking_set_retrieve_with_scores_matches_retrieve_idx():
    store = torch.eye(4, dtype=torch.float32)
    ranker = RankingSet(store=store, metric="cosine")

    query = torch.tensor([1.0, 0.0, 0.0, 0.0], dtype=torch.float32)
    sims, idxs = ranker.retrieve_with_scores(query, n=2)

    # Top hit should be index 0 with similarity 1.0
    assert sims[0] == pytest.approx(1.0)
    assert idxs[0] == 0


def test_preprocess_raw_inputs_hsqc_shape_and_order():
    raw = {"hsqc": [1.0, 2.0, 0.5, 3.0, 4.0, 0.8]}
    processed = _preprocess_raw_inputs(raw)
    hsqc = processed["hsqc"]

    assert hsqc.shape == (2, 3)
    # Check that columns are [C, H, intensity]
    assert torch.allclose(hsqc[0], torch.tensor([2.0, 1.0, 0.5]))
    assert torch.allclose(hsqc[1], torch.tensor([4.0, 3.0, 0.8]))


