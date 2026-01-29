"""Tests for model registry."""
import pytest

from src.services.model_registry import get


def test_get_returns_none_when_not_loaded():
    assert get("nonexistent_model") is None
