"""Tests for retrieval helpers."""
import pytest
import torch

from src.services.retrieval_helpers import dense_from_indices, kept_indices_for_mw_range


def test_dense_from_indices():
    vec = dense_from_indices(10, [0, 2, 4])
    assert vec.shape == (10,)
    assert vec[0].item() == 1.0
    assert vec[2].item() == 1.0
    assert vec[4].item() == 1.0
    assert vec[1].item() == 0.0


def test_dense_from_indices_empty():
    vec = dense_from_indices(5, [])
    assert vec.shape == (5,)
    assert vec.sum().item() == 0.0


def test_kept_indices_for_mw_range_no_filter():
    class MockSession:
        def indices_in_mw_range(self, mw_min, mw_max):
            return [0, 1, 2]

    session = MockSession()
    assert kept_indices_for_mw_range(session, None, None) == [0, 1, 2]


def test_kept_indices_for_mw_range_with_filter():
    class MockSession:
        def indices_in_mw_range(self, mw_min, mw_max):
            return [1, 3]

    session = MockSession()
    assert kept_indices_for_mw_range(session, 100.0, 500.0) == [1, 3]
