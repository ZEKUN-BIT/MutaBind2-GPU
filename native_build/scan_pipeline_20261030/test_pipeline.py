import argparse
import csv
import fcntl
import json
from pathlib import Path
import sys
import threading
import time

import pytest

HERE = Path(__file__).resolve().parent
INSTALL = HERE.parents[1]
sys.path[:0] = [str(HERE), str(INSTALL)]
import pipeline
from pipeline import Pipeline, core, read
from pipeline_worker import wait_for_energy
from recover_energy import validate_prepared
from mutabind2_tool.scan import plan
from mutabind2_tool.core.runtime import ALLOWED_KEYS, load_config

FAKE_WORKER = '''
import argparse,csv,fcntl,json,os,sys,time
from pathlib import Path
sys.path.insert(0,WORKER_DIRECTORY)
from pipeline_worker import wait_for_energy
from mutabind2_tool.scan import task
p=argparse.ArgumentParser();p.add_argument('module')
for name in ['source','workspace','parameters','config','baseline']:p.add_argument('--'+name,type=Path)
a=p.parse_args();w=a.workspace;w.mkdir();(w/'job').mkdir();(w/'jobout').mkdir()
root=a.source.parents[2];batch=int(a.source.parent.name.rsplit('_',1)[1])
parameters=json.loads(a.parameters.read_text())
def event(name):
    fd=os.open(root/'test_events.jsonl',os.O_WRONLY|os.O_CREAT|os.O_APPEND,0o600)
    os.write(fd,(json.dumps(dict(event=name,batch=batch,time=time.monotonic()))+'\\n').encode());os.close(fd)
with a.source.open() as stream:rows=list(csv.DictReader(stream,delimiter='\\t'))
state=dict(stage='preparing',protocol_sha256=parameters['protocol_sha256'],pid=os.getpid(),
           result_ids=[r['Result_Id'] for r in rows])
task.atomic_json(w/'status.json',state);event('prep_start');time.sleep(.15)
if (root/'fail_batch').exists() and int((root/'fail_batch').read_text())==batch:
    state.update(stage='failed',error='injected preparation failure');task.atomic_json(w/'status.json',state);sys.exit(1)
event('prep_end');wait_for_energy(w,state,'mutant_fake')
with (root/'test_energy.lock').open('a') as lock:
    fcntl.flock(lock,fcntl.LOCK_EX)
    state.update(stage='mutant_pb_energy_fake',pb_runs={'mutant_fake':{'outcome':'running'}})
    task.atomic_json(w/'status.json',state);event('energy_start');time.sleep(.6);event('energy_end')
    values=[dict(r,SunDDG=str(int(r['Result_Id'])/10),Confidence='High',Interface='no') for r in rows]
    with (w/'results.tsv').open('w') as stream:
        writer=csv.DictWriter(stream,list(values[0]),delimiter='\\t',lineterminator='\\n');writer.writeheader();writer.writerows(values)
    (w/'job/job.input.cleaned.outdata').write_text('fixed scientific features\\n')
    state.update(stage='complete');state['pb_runs']['mutant_fake']['outcome']='complete';task.atomic_json(w/'status.json',state)
'''


@pytest.fixture
def factory(tmp_path, monkeypatch):
    monkeypatch.setattr(pipeline.task, 'validate_baseline', lambda *a, **k: {})
    monkeypatch.setattr(pipeline.energy, '_available', lambda: 35 * 1024**3)
    serial = 0
    def make():
        nonlocal serial
        serial += 1
        output = tmp_path / f'scan_{serial}'
        mutations = tmp_path / f'mutations_{serial}.tsv'
        mutations.write_text('chain\tmutation\nI\tL45P\nI\tL45A\nI\tL45W\n')
        args = argparse.Namespace(output=output, pdb=INSTALL/'examples/1cse.pdb', partner1='E', partner2='I',
            mutations=mutations, sites=None, mutated_partner=None, fragment_mode='auto', fragments=None,
            gap_threshold=2.0, cap_mode='charged', surface_backend='local', provean_mode='fragment',
            provean_cache=None, max_per_batch=1, workers=2, foldx_workers=2, pb_executable=None,
            vdw_executable=None, backend='gpu', config=INSTALL/'config/tools.example.json', no_annotations=False,
            resume=False, threads=8, no_contributions=False, combined=None, handoff=False,
            keep_workspace=False, poll=.01)
        manifest, _ = plan.build(args)
        config = load_config(args.config)
        core.atomic_json(output/'config.snapshot.json', {k: config.get(k) for k in ALLOWED_KEYS})
        (output/'provean').mkdir();(output/'provean/scores.tsv').write_text('fixture scores\n')
        identity, protocol = core.protocol(manifest, config, output/'provean')
        core.atomic_json(output/'protocol.json', dict(sha256=protocol, identity=identity))
        (output/'baseline').mkdir()
        result = Pipeline(args)
        fake = tmp_path / f'worker_{serial}.py'
        fake.write_text(FAKE_WORKER.replace('WORKER_DIRECTORY', repr(str(HERE))))
        result.worker = fake
        return result
    return make


def finish(worker, handle):
    deadline = time.monotonic() + 15
    while worker.alive(handle):
        assert time.monotonic() < deadline
        time.sleep(.01)
    assert read(handle['workspace']/'status.json')['stage'] == 'complete'
    core._archive_batch(handle['workspace'], worker.destination(handle['batch']), worker.root,
                        worker.expected(handle['batch']), worker.protocol)


def events(worker):
    return [json.loads(line) for line in (worker.output/'test_events.jsonl').read_text().splitlines()]


def test_pipeline_overlaps_preparation_but_serializes_energy_and_matches_serial_results(factory):
    serial = factory()
    for batch in serial.batches:
        handle = serial.start(batch);serial.release(handle);finish(serial, handle)
    expected = core.aggregate(serial.output, serial.manifest)
    parallel = factory();parallel.run()
    assert core.read_table(parallel.output/'predictions.tsv') == expected
    steps = events(parallel)
    energy = [e for e in steps if e['event'].startswith('energy_')]
    assert [e['event'] for e in energy] == ['energy_start', 'energy_end'] * 3
    lookup = {(e['batch'],e['event']):e['time'] for e in steps}
    assert lookup[2,'prep_start'] < lookup[1,'energy_end']
    assert lookup[2,'energy_start'] >= lookup[1,'energy_end']
    assert read(parallel.output/'batch_status.json')['stage'] == 'complete'
    assert all((parallel.destination(b)/'scientific_intermediates.tar.gz').exists() for b in parallel.batches)


def test_prepared_worker_is_adopted_after_controller_restart_without_repeating_preparation(factory):
    first = factory();current=first.start(first.batches[0]);first.release(current)
    future=first.start(first.batches[1])
    deadline=time.monotonic()+10
    while read(future['workspace']/'status.json').get('stage')!='waiting_pipeline_energy':
        assert time.monotonic()<deadline;time.sleep(.01)
    finish(first,current)
    resumed=Pipeline(first.args);resumed.worker=first.worker
    try:
        resumed.run()
    finally:
        first.stop_children()
    assert len([e for e in events(resumed) if e['batch']==2 and e['event']=='prep_start'])==1
    assert len(core.read_table(resumed.output/'completed_predictions.tsv'))==3


def test_failed_preparation_stops_pipeline_without_releasing_its_energy(factory):
    worker=factory();(worker.output/'fail_batch').write_text('2')
    with pytest.raises(RuntimeError,match='preparation failed'):
        worker.run()
    assert read(worker.output/'batch_status.json')['stage']=='failed'
    assert not any(e['event']=='energy_start' and e['batch']==2 for e in events(worker))
    assert all(p.poll() is not None for p in worker.children.values())


def test_stop_marker_and_existing_scan_lock_prevent_new_calculation(factory):
    worker=factory();worker.marker.write_text('user stop\n');worker.run()
    assert read(worker.output/'batch_status.json')['stage']=='paused'
    assert not worker.children
    other=factory()
    with (other.output/'scan.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        with pytest.raises(BlockingIOError):other.run()
    assert not other.children


def test_wrong_protocol_gate_is_rejected(tmp_path):
    (tmp_path/'pipeline_energy_ready.json').write_text(json.dumps({'protocol_sha256':'wrong'}))
    with pytest.raises(ValueError,match='different protocol'):
        wait_for_energy(tmp_path,{'protocol_sha256':'expected'},'energy')


def test_low_memory_defers_prefetch(factory,monkeypatch):
    worker=factory();workspace=worker.root/'batch_0001_attempt_001';workspace.mkdir()
    core.atomic_json(workspace/'status.json',{'pb_runs':{'mutant_fake':{}}})
    monkeypatch.setattr(pipeline.energy,'_available',lambda:29*1024**3)
    assert worker.maybe_prepare({'workspace':workspace},worker.batches[1]) is None
    assert not worker.children


def test_clean_handoff_archives_current_batch_before_releasing_prefetched_batch(factory):
    worker=factory();original=worker.start(worker.batches[0]);worker.release(original)
    core.update(worker.output,stage='running',batch=1,active_pid=original['pid'],workspace=str(original['workspace']))
    lock=(worker.output/'scan.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX)
    errors=[]
    def old_scheduler():
        try:
            finish(worker,original)
            core.aggregate(worker.output,worker.manifest)
            core.update(worker.output,stage='paused',active_pid=None)
        except BaseException as exc:errors.append(exc)
        finally:lock.close()
    thread=threading.Thread(target=old_scheduler);thread.start()
    worker.args.handoff=True
    try:worker.run()
    finally:thread.join(15);worker.stop_children()
    assert not errors
    assert not worker.marker.exists()
    steps=events(worker);lookup={(e['batch'],e['event']):e['time'] for e in steps}
    assert lookup[2,'prep_start']<lookup[1,'energy_end']
    assert lookup[2,'energy_start']>=lookup[1,'energy_end']
    assert read(worker.output/'batch_status.json')['stage']=='complete'


def test_recover_handoff_after_old_scheduler_already_paused(factory):
    worker=factory();worker.marker.write_text('Pipeline handoff old\n')
    core.atomic_json(worker.output/'pipeline_handoff.json',{'marker_text':worker.marker.read_text()})
    core.update(worker.output,stage='paused',batch=1,completed=0)
    worker.args.handoff=True;worker.run()
    assert not worker.marker.exists()
    assert read(worker.output/'batch_status.json')['stage']=='complete'


def test_energy_recovery_refuses_changed_protocol_before_touching_workspace(tmp_path):
    parameters = argparse.Namespace(protocol_sha256='expected', cap_mode='neutral')
    with pytest.raises(ValueError, match='same protocol and baseline'):
        validate_prepared(tmp_path, tmp_path/'source', tmp_path/'baseline', parameters,
                          {'protocol_sha256': 'changed'})


def test_energy_recovery_refuses_changed_staged_input(tmp_path):
    (tmp_path/'job').mkdir()
    (tmp_path/'job/job.input').write_text('changed')
    source=tmp_path/'source';source.write_text('original')
    baseline=tmp_path/'baseline'
    parameters=argparse.Namespace(protocol_sha256='expected', cap_mode='neutral')
    state=dict(protocol_sha256='expected',cap_mode='neutral',baseline=str(baseline),pb_runs={'a':{}})
    with pytest.raises(ValueError,match='immutable batch'):
        validate_prepared(tmp_path,source,baseline,parameters,state)
