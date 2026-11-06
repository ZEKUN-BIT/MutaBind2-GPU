from pathlib import Path
import json
import pytest
from mutabind2_tool.scan import scheduler


def test_mixed_case_checkout_uses_home_cache(tmp_path, monkeypatch):
    monkeypatch.delenv('MUTABIND2_WORKSPACE_ROOT', raising=False)
    monkeypatch.setattr(scheduler, 'TOOL_ROOT', tmp_path / 'MixedCase')
    monkeypatch.setattr(Path, 'home', classmethod(lambda cls: tmp_path / 'home'))
    output = tmp_path / 'result'
    root = scheduler.workspace_root(output)
    assert root.parent == tmp_path / 'home/.cache/mutabind2/workspaces'
    assert json.loads((root / 'owner.json').read_text())['scan'] == str(output)
    assert scheduler.workspace_root(output) == root


def test_invalid_explicit_scratch_is_not_silently_redirected(tmp_path):
    with pytest.raises(ValueError, match='lowercase ASCII'):
        scheduler.workspace_root(tmp_path / 'result', tmp_path / 'MixedCase')
