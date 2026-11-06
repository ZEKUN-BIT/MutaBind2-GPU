import argparse
import importlib.util
import json
from pathlib import Path
import sys

import pytest

HERE = Path(__file__).resolve().parent
INSTALL = HERE.parents[1]
sys.path[:0] = [str(HERE), str(INSTALL / 'native_build/scan_pipeline_20261030'), str(INSTALL)]
import memory_admission as memory
import recover
from test_pipeline import factory, events

spec = importlib.util.spec_from_file_location('memory_continue', HERE / 'continue_scan.py')
controller = importlib.util.module_from_spec(spec)
spec.loader.exec_module(controller)


def setup_surface(tmp_path, monkeypatch, available, run):
    (tmp_path / 'status.json').write_text(json.dumps({'stage': 'surface', 'protocol_sha256': 'unchanged'}))
    samples = iter(available)
    sleeps = []
    monkeypatch.setattr(memory.surface, '_available', lambda: next(samples) * memory.surface.GIB)
    monkeypatch.setattr(memory.surface, 'run_deck', run)
    monkeypatch.setattr(memory.time, 'sleep', lambda seconds: sleeps.append(seconds))
    memory.install()
    return sleeps


def test_surface_waits_for_memory_and_preserves_original_call_and_state(tmp_path, monkeypatch):
    calls = []
    def run(*args, **kwargs):
        calls.append((args, kwargs))
        return {'native_result': 'unchanged'}
    sleeps = setup_surface(tmp_path, monkeypatch, [29, 25, 30], run)
    result = memory.surface.run_deck('charmm', tmp_path, 'stem', deck='original deck')
    assert result == {'native_result': 'unchanged'}
    assert calls == [(('charmm', tmp_path, 'stem'), {'deck': 'original deck'})]
    assert sleeps == [1, 1]
    state = json.loads((tmp_path / 'status.json').read_text())
    assert state['stage'] == 'surface' and state['protocol_sha256'] == 'unchanged'
    assert state['surface_memory_admission']['stage'] == 'admitted'


def test_surface_retries_admission_race_before_native_launch(tmp_path, monkeypatch):
    calls = []
    def run(*args, **kwargs):
        calls.append(1)
        if len(calls) == 1:
            raise RuntimeError(memory.ADMISSION_ERROR)
        return 'success'
    sleeps = setup_surface(tmp_path, monkeypatch, [31, 29, 31], run)
    assert memory.surface.run_deck('charmm', tmp_path) == 'success'
    assert calls == [1, 1] and sleeps == [1]


def test_surface_native_errors_are_not_retried_or_hidden(tmp_path, monkeypatch):
    calls = []
    def run(*args, **kwargs):
        calls.append(1)
        raise RuntimeError('Surface CHARMM exit code 1')
    sleeps = setup_surface(tmp_path, monkeypatch, [35], run)
    with pytest.raises(RuntimeError, match='exit code 1'):
        memory.surface.run_deck('charmm', tmp_path)
    assert calls == [1] and not sleeps


def test_memory_wait_is_interruptible_and_never_starts_native_work(tmp_path, monkeypatch):
    calls = []
    setup_surface(tmp_path, monkeypatch, [25], lambda *a: calls.append(1))
    def interrupt(seconds):
        raise InterruptedError('stop')
    monkeypatch.setattr(memory.time, 'sleep', interrupt)
    with pytest.raises(InterruptedError):
        memory.surface.run_deck('charmm', tmp_path)
    assert not calls
    assert json.loads((tmp_path / 'status.json').read_text())['stage'] == 'waiting_surface_memory'


def test_guard_installation_is_idempotent(tmp_path, monkeypatch):
    setup_surface(tmp_path, monkeypatch, [35], lambda *a: 'done')
    original = memory.surface.run_deck
    memory.install()
    assert memory.surface.run_deck is original


def test_prefetch_failure_finishes_and_saves_current_batch(factory):
    worker = factory()
    worker.__class__ = controller.MemoryPipeline
    (worker.output / 'fail_batch').write_text('2')
    with pytest.raises(RuntimeError, match='batch 1 completed and saved'):
        worker.run()
    assert worker.pending()[0]['index'] == 2
    assert len(controller.core.read_table(worker.output / 'completed_predictions.tsv')) == 1
    steps = events(worker)
    assert any(e['batch'] == 1 and e['event'] == 'energy_end' for e in steps)
    assert not any(e['batch'] == 2 and e['event'] == 'energy_start' for e in steps)
    assert all(p.poll() is not None for p in worker.children.values())


def test_fixed_pipeline_still_overlaps_preparation_and_serializes_energy(factory):
    worker = factory()
    worker.__class__ = controller.MemoryPipeline
    worker.run()
    steps = events(worker)
    lookup = {(e['batch'], e['event']): e['time'] for e in steps}
    assert lookup[2, 'prep_start'] < lookup[1, 'energy_end']
    assert lookup[2, 'energy_start'] >= lookup[1, 'energy_end']
    assert len(controller.core.read_table(worker.output / 'completed_predictions.tsv')) == 3


def test_fixed_pipeline_preserves_stop_marker(factory):
    worker = factory()
    worker.__class__ = controller.MemoryPipeline
    worker.marker.write_text('user stop\n')
    worker.run()
    assert worker.marker.read_text() == 'user stop\n'
    assert not worker.children


def test_recovery_rejects_protocol_change_without_touching_artifacts(tmp_path):
    parameters = argparse.Namespace(protocol_sha256='expected', cap_mode='neutral')
    with pytest.raises(ValueError, match='same protocol'):
        recover.validate_prepared(tmp_path, tmp_path / 'source', tmp_path / 'baseline',
                                  parameters, {'protocol_sha256': 'changed'})
    assert not list(tmp_path.iterdir())


def test_recovery_accepts_pre_surface_checkpoint_and_rejects_changed_geometry(tmp_path, monkeypatch):
    workspace = tmp_path / 'workspace'
    for name in ('job', 'jobout', 'sinputfiles'):
        (workspace / name).mkdir(parents=True)
    source = tmp_path / 'batch.input'
    source.write_text('immutable batch\n')
    (workspace / 'job/job.input').write_bytes(source.read_bytes())
    (tmp_path / 'structure.pdb').write_text('same structure\n')
    (workspace / 'job/structure.pdb').write_text('same structure\n')
    table = [{'Result_Id': '1', 'PDBid': 'fragments', 'Mutation_cleaned': 'CE43A'}]
    recover.recovery.core.write_table(workspace / 'job/job.input.cleaned', list(table[0]), table)
    coordinate = workspace / 'jobout/fragments_ce43a_min.crd'
    coordinate.write_text('frozen minimized coordinates\n')
    (workspace / 'jobout/min_fragments_ce43a.log').write_text('ENERGY: 100 0\n')
    monkeypatch.setattr(recover.recovery.task, 'inspect_job',
                        lambda *a: ([dict(Result_Id='1', PDBfile='structure.pdb')], None))
    monkeypatch.setattr(recover.recovery.task, 'validate_baseline', lambda *a: {})
    monkeypatch.setattr(recover.recovery.task, '_expected_residues', lambda *a: {})
    monkeypatch.setattr(recover.recovery.task, 'audit_geometry', lambda *a: {})
    monkeypatch.setattr(recover.recovery.task, 'read_foldx_difference', lambda *a: 0)
    baseline = tmp_path / 'baseline'
    parameters = argparse.Namespace(protocol_sha256='same', cap_mode='neutral')
    state = dict(protocol_sha256='same', cap_mode='neutral', baseline=str(baseline),
        timings_seconds={name: 1 for name in ('foldx', 'mutant_vmd', 'mutant_charmm',
            'energy_decks', 'surface_decks', 'minimize_fragments_ce43a', 'geometry_fragments_ce43a')})
    assert recover.freeze(workspace, source, baseline, parameters, state) == table
    assert 'surface' not in state['timings_seconds'] and not state.get('pb_runs')
    coordinate.write_text('changed coordinates\n')
    with pytest.raises(ValueError, match='Frozen recovery artifact changed'):
        recover.validate_prepared(workspace, source, baseline, parameters, state)
