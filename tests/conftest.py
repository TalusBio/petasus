"""Test fixtures"""

from pathlib import Path

import pytest


@pytest.fixture
def data_path():
    """The path to test data."""
    return Path(__file__).parent / "data"


@pytest.fixture(autouse=True)
def chdir(monkeypatch, tmp_path):
    """Run tests in test dir."""
    monkeypatch.chdir(tmp_path)
