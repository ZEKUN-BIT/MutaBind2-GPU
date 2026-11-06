"""Exercise provenance merging and a real, isolated frozen-geometry handoff."""
import argparse
import importlib.util
import json
from pathlib import Path
import shutil

import pytest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('cells_continuation',HERE/'continue_scan.py')
cells=importlib.util.module_from_spec(spec)
spec.loader.exec_module(cells)


@pytest.fixture
def scans(tmp_path,monkeypatch):
    for key in ('OLD','WARP','CELLS','NEW'):
        path=tmp_path/key.lower();path.mkdir()
        monkeypatch.setattr(cells,key,path)
        cells.core.atomic_json(path/'protocol.json',dict(sha256=key))
    rows=[dict(Result_Id=str(i),Chain='L',Mutation='V130'+mut,Status='pending')
          for i,mut in enumerate('AEKR',1)]
    cells.core.atomic_json(cells.OLD/'scan.json',dict(mutation_identity=[
        [r['Result_Id'],r['Chain'],r['Mutation']] for r in rows]))
    cells.core.write_table(cells.OLD/'predictions.tsv',list(rows[0]),rows)
    for index,scan in enumerate((cells.OLD,cells.WARP,cells.CELLS,cells.NEW)):
        row=dict(rows[index],Status='complete')
        cells.core.write_table(scan/'completed_predictions.tsv',list(row),[row])
    return rows


def test_merge_keeps_original_order_and_four_protocols(scans):
    rows,state=cells.collect()
    assert [r['Result_Id'] for r in rows]==['1','2','3','4']
    assert [r['Protocol_sha256'] for r in rows]==['OLD','WARP','CELLS','NEW']
    assert state['stage']=='complete' and state['completed']==4


def test_duplicate_completed_id_is_rejected(scans):
    row=dict(scans[0],Status='complete')
    cells.core.write_table(cells.NEW/'completed_predictions.tsv',list(row),[row])
    with pytest.raises(ValueError,match='Duplicate'):
        cells.collect()


def test_wrong_mutation_identity_is_rejected(scans):
    row=dict(scans[3],Status='complete',Mutation='V130W')
    cells.core.write_table(cells.NEW/'completed_predictions.tsv',list(row),[row])
    with pytest.raises(ValueError,match='mismatched'):
        cells.collect()


def test_missing_protocol_is_rejected(scans):
    (cells.NEW/'protocol.json').unlink()
    with pytest.raises(ValueError,match='Invalid'):
        cells.collect()


def test_real_dry_handoff_retains_wt_and_all_prepared_geometries(tmp_path,monkeypatch):
    from resume_prepared import validate
    pipeline=cells.pipeline
    current=cells.read(cells.CELLS/'pipeline_status.json')
    donor=Path(current['preparation']['workspace'])
    donor_state=cells.read(donor/'status.json')
    assert donor_state['stage']=='waiting_pipeline_energy' and not donor_state.get('pb_runs')
    original_rows,_=cells.collect()
    batches=cells.read(cells.CELLS/'batches.json')['batches']
    finishing=next(b['result_ids'] for b in batches if b['index']==current['energy']['batch'])
    rows=[dict(r,Status='complete',Source_scan=str(cells.CELLS))
          if r['Result_Id'] in finishing else r for r in original_rows]
    actual_read=cells.read
    def assume_boundary(path):
        if path==cells.CELLS/'batch_status.json':
            return dict(stage='paused')
        return actual_read(path)
    monkeypatch.setattr(cells,'read',assume_boundary)
    monkeypatch.setattr(cells,'collect',lambda:(rows,{}))
    evidence=tmp_path/'evidence';evidence.mkdir()
    (evidence/'pb_flags').symlink_to(HERE/'pb_flags',target_is_directory=True)
    (evidence/'vdw').symlink_to(HERE/'vdw',target_is_directory=True)
    monkeypatch.setattr(cells,'HERE',evidence)
    monkeypatch.setattr(cells,'NEW',tmp_path/'dry_scan')
    monkeypatch.setattr(cells,'BACKEND',HERE/'pb_flags/bin/charmm')
    monkeypatch.setattr(cells,'VDW_BACKEND',HERE/'vdw/production/bin/charmm_gpu')
    manifest,protocol=cells.prepare()
    baseline=cells.NEW/'baseline'
    original=cells.CELLS/'baseline'
    saved=cells.read(original/'status.json')
    assert all(cells.digest(baseline/'wt_cache'/name)==value
               for name,value in saved['wt_cache_sha256'].items())
    assert (donor/'job/job.input').read_bytes()==(cells.NEW/'batches/batch_0001/batch.input').read_bytes()
    cloned=tmp_path/'frozen_preparation'
    shutil.copytree(donor,cloned,symlinks=True)
    hashes={str(p.relative_to(cloned)):cells.digest(p) for folder in ('job','jobout','sinputfiles')
            for p in (cloned/folder).rglob('*') if p.is_file()}
    pipeline.core.atomic_json(cloned/'gpu_migration.json',dict(new_protocol=protocol,scientific_files=hashes))
    donor_state.update(protocol_sha256=protocol,baseline=str(baseline))
    parameters=pipeline.core.task_arguments(manifest,argparse.Namespace(threads=8,no_contributions=False),
                                            cells.NEW/'provean',protocol)
    table=validate(cloned,cells.NEW/'batches/batch_0001/batch.input',baseline,
                   argparse.Namespace(**parameters),donor_state)
    assert len(table)==19
    # The same checkpoint must be rejected if a donated coordinate is altered.
    coordinate=cloned/'jobout'/f"{table[0]['PDBid'].lower()}_{table[0]['Mutation_cleaned'].lower()}_min.crd"
    with coordinate.open('a') as stream:
        stream.write('\nchanged checkpoint\n')
    with pytest.raises(ValueError,match='Donated artifact changed'):
        validate(cloned,cells.NEW/'batches/batch_0001/batch.input',baseline,
                 argparse.Namespace(**parameters),donor_state)
    report=dict(passed=True,frozen_wt_files=len(saved['wt_cache_sha256']),
                frozen_mutants=len(table),corrupted_checkpoint_rejected=True)
    (HERE/'transition_test_results.json').write_text(json.dumps(report,indent=2)+'\n')


def test_resume_pending_set_keeps_own_completed_ids(scans):
    rows,_=cells.collect()
    pending=[r for r in rows if r['Status']!='complete' or r.get('Source_scan')==str(cells.NEW)]
    assert [r['Result_Id'] for r in pending]==['4']
