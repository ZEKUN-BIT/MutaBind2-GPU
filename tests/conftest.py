"""Run planning tests without a local native-tool installation."""
from pathlib import Path

import pytest


@pytest.fixture(autouse=True)
def example_configuration(monkeypatch):
    config = Path(__file__).resolve().parents[1] / "config/tools.example.json"
    monkeypatch.setenv("MUTABIND2_CONFIG", str(config))
