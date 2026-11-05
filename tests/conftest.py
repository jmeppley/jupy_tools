"""
Shared pytest configuration and fixtures for the jupy_tools test suite.

Run the whole suite from the repository root with:

    pytest
"""
import os

import pytest

# Never try to open a display when plotting in tests
os.environ.setdefault("MPLBACKEND", "Agg")

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TESTS_DIR)


@pytest.fixture(autouse=True)
def _run_from_repo_root(monkeypatch):
    """ Some tests use paths relative to the repo root (eg 'tests/test_fs.py') """
    monkeypatch.chdir(REPO_ROOT)


@pytest.fixture
def close_figures():
    """ close any matplotlib figures a test creates """
    yield
    from matplotlib import pyplot as plt
    plt.close("all")


@pytest.fixture
def ax(close_figures):
    """ a fresh matplotlib Axes to draw on """
    from matplotlib import pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=100)
    return ax

