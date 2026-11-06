"""Run operational regression tests without installed native tools or production scratch."""
import pytest
from mutabind2_tool.scan import scheduler

@pytest.fixture(autouse=True)
def isolated_pipeline_environment(tmp_path, monkeypatch):
    monkeypatch.setenv('MUTABIND2_WORKSPACE_ROOT', str(tmp_path / 'scratch'))
    monkeypatch.setattr(scheduler, 'protocol', lambda *args: ({'test_fixture': True}, 'isolated-test-protocol'))
