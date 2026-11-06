"""Continue remaining IDs with validated GPU lists/MSTEP and the exact frozen WT."""
from pathlib import Path
import argparse
import fcntl
import hashlib
import json
import os
import shlex
import shutil
import signal
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
INSTALL=HERE.parents[1]
PIPE=INSTALL/'native_build/scan_pipeline_20261030'
sys.path[:0]=[str(PIPE),str(INSTALL)]
import pipeline
from mutabind2_tool.scan.provean import prepare_scan_cache

OLD=INSTALL/'projects/7mkj_128_255_20261030'
WARP=INSTALL/'projects/7mkj_128_255_20261030_optimized'
CELLS=INSTALL/'projects/7mkj_128_255_20261030_cells'
NEW=INSTALL/'projects/7mkj_128_255_20261030_gpu'
COMBINED=INSTALL/'projects/7mkj_128_255_20261030_combined'
BACKEND=INSTALL/'bin/charmm-pb-gpu-mstep'
VDW_BACKEND=INSTALL/'bin/charmm-vdw-gpu-lists'
read=pipeline.read
digest=pipeline.digest
core=pipeline.core


def collect():
    original=read(OLD/'scan.json')
    expected={rid:(chain,mutation) for rid,chain,mutation in original['mutation_identity']}
    rows=core.read_table(OLD/'predictions.tsv')
    completed={}
    for scan in (OLD,WARP,CELLS,NEW):
        path=scan/'completed_predictions.tsv'
        if not path.exists():
            continue
        protocol=read(scan/'protocol.json').get('sha256')
        for row in core.read_table(path):
            rid=row['Result_Id']
            if rid in completed or expected.get(rid)!=(row['Chain'],row['Mutation']):
                raise ValueError('Duplicate or mismatched result: '+rid)
            if row['Status']!='complete' or not protocol:
                raise ValueError('Invalid completed prediction')
            completed[rid]=dict(row,Source_scan=str(scan),Protocol_sha256=protocol)
    merged=[completed.get(r['Result_Id'],dict(r,Source_scan='',Protocol_sha256='')) for r in rows]
    state=read(NEW/'batch_status.json') or read(CELLS/'batch_status.json')
    if state.get('stage')=='running' and state.get('workspace'):
        native=read(Path(state['workspace'])/'status.json')
        state.update(native_stage=native.get('stage'),current_mutation=native.get('pb_current_mutation'),
            current_energies_complete=sum(k.startswith('mutant_') and v.get('outcome')=='complete'
                for k,v in native.get('pb_runs',{}).items()))
    state.update(total=len(expected),completed=len(completed),pending=len(expected)-len(completed),
                 active_scan=str(NEW if (NEW/'protocol.json').exists() else CELLS),
                 combined_results=str(COMBINED/'completed_predictions.tsv'))
    if len(completed)==len(expected):
        state['stage']='complete'
    return merged,state


def aggregate():
    rows,state=collect()
    core.write_table(COMBINED/'predictions.tsv',list(rows[0]),rows)
    core.write_table(COMBINED/'completed_predictions.tsv',list(rows[0]),
                     [r for r in rows if r['Status']=='complete'])
    core.atomic_json(COMBINED/'status.json',state)
    return state


def install():
    pb_report=read(HERE/'pb_flags/validation/comparison.json')
    if set(pb_report)!={'wt','mutant'} or not all(v['bitwise_phi_equal'] and v['logged_results_equal']
                                                for v in pb_report.values()):
        raise ValueError('Frozen WT and mutant PB validation must pass')
    checks=['pb_flags/test_results.json','vdw/test_results.json','vdw/validation/small_report.json',
            'vdw/validation/large_audit_report.json','vdw/validation/large_performance_report.json',
            'transition_test_results.json']
    if not all(read(HERE/name).get('passed') for name in checks):
        raise ValueError('Native equivalence and continuation tests must pass')
    runtimes={}
    for mode,build,wrapper_name,install_name in [
            ('pb',HERE/'pb_flags',BACKEND,'charmm-pb-cuda-mstep'),
            ('vdw',HERE/'vdw/production',VDW_BACKEND,'charmm-vdw-cuda-lists')]:
        library=build/'lib/libchmm.so'
        if digest(library)!=read(build/'build_manifest.json')['library_sha256']:
            raise ValueError('Candidate library changed after build')
        cases=[v['backends']['mstep'] for v in pb_report.values()] if mode=='pb' else [
            v['backends']['gpu'] for v in read(HERE/'vdw/validation/large_performance_report.json')['cases']]
        if len(cases)!=2 or not all(c['runtime']['files'].get(str(library))==digest(library) for c in cases):
            raise ValueError('Validation used different libraries')
        target=INSTALL/'toolchain/external'/install_name/'lib/libchmm.so'
        target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists():
            if digest(target)!=digest(library):raise ValueError('Refusing to overwrite a deployed backend')
        else:shutil.copy2(library,target)
        ldpath=str(target.parent)+':/usr/local/cuda/lib64:/home/zek/miniforge3/lib'
        flags='export MUTABIND2_PB_BACKEND=cuda\n' if mode=='pb' else 'export MUTABIND2_VDW_COMPACT=1\nexport MUTABIND2_VDW_GPU=1\n'
        wrapper=('#!/bin/sh\nset -eu\nexport LD_LIBRARY_PATH='+shlex.quote(ldpath)
            +'\nexport OMP_NUM_THREADS=1\n'+flags+'exec '
            +str(INSTALL/'toolchain/external/charmm-51/bin/charmm')+' "$@"\n')
        if wrapper_name.exists() and wrapper_name.read_text()!=wrapper:
            raise ValueError('Refusing to overwrite a different wrapper')
        wrapper_name.write_text(wrapper);wrapper_name.chmod(0o755)
        runtime=pipeline.energy._runtime_identity(wrapper_name)
        if runtime['files'].get(str(target))!=digest(library):raise ValueError('Deployed loader resolved a different library')
        runtimes[mode]=runtime
    core.atomic_json(HERE/'deployment.json',dict(deployed_utc=core.now(),runtimes=runtimes,
        evidence_sha256={name:digest(HERE/name) for name in checks+['pb_flags/validation/comparison.json']}))


def prepare():
    if read(CELLS/'batch_status.json').get('stage')!='paused':
        raise ValueError('The preceding scan must finish its current batch first')
    old=pipeline.plan.load(CELLS)
    rows,_=collect()
    pending=[dict(id=r['Result_Id'],chain=r['Chain'],mutation=r['Mutation'])
             for r in rows if r['Status']!='complete' or r.get('Source_scan')==str(NEW)]
    if not pending:
        raise ValueError('No pending mutations')
    input_path=HERE/'remaining_mutations.tsv'
    core.write_table(input_path,['id','chain','mutation'],pending)
    if not (NEW/'scan.json').exists():
        command=[str(INSTALL/'mutabind2'),'scan','plan','--pdb',old['source'],
            '--partner1',','.join(old['partner1']),'--partner2',','.join(old['partner2']),
            '--mutations',str(input_path),'--output',str(NEW),'--pb-executable',str(BACKEND),
            '--vdw-executable',str(VDW_BACKEND),'--provean-cache',old['provean_cache'],
            '--config',str(CELLS/'config.snapshot.json')]
        for key in ('fragment_mode','gap_threshold','cap_mode','surface_backend','provean_mode',
                    'max_per_batch','workers','foldx_workers','backend'):
            command+=['--'+key.replace('_','-'),str(old[key])]
        subprocess.run(command,cwd=INSTALL,check=True)
    manifest=pipeline.plan.load(NEW)
    if manifest['mutation_identity']!=[[r['id'],r['chain'],r['mutation']] for r in pending]:
        raise ValueError('Continuation plan does not exactly cover the remaining requests')
    if manifest['plan']['prepared_pdb_sha256']!=old['plan']['prepared_pdb_sha256']:
        raise ValueError('Continuation changed the prepared structure')
    shutil.copy2(CELLS/'config.snapshot.json',NEW/'config.snapshot.json')
    config=pipeline.load_config(NEW/'config.snapshot.json')
    cache=prepare_scan_cache(NEW,manifest,config,8)
    identity,protocol=core.protocol(manifest,config,cache)
    prior=read(CELLS/'protocol.json')
    if {k:v for k,v in identity.items() if k not in ('pb','vdw','provean_scores')}!={
        k:v for k,v in prior['identity'].items() if k not in ('pb','vdw','provean_scores')}:
        raise ValueError('Unvalidated change outside validated PB/VDW backends and conservation subset')
    core.atomic_json(NEW/'protocol.json',dict(sha256=protocol,identity=identity))
    baseline=NEW/'baseline'
    original=CELLS/'baseline'
    pipeline.task.validate_baseline(original,original/'baseline.input',manifest['cap_mode'],prior['sha256'])
    if not (baseline/'status.json').exists():
        shutil.copytree(original,baseline)
        requests=core.read_table(NEW/'plan/mutations.input')
        core.write_table(baseline/'baseline.input',list(requests[0]),[requests[0]])
        # Preserve original calculation provenance and all cache hashes. This is
        # an explicitly validated donation; no new minimization is claimed.
        saved=read(baseline/'status.json')
        saved.update(protocol_sha256=protocol,donated_from=str(original),
            original_protocol_sha256=prior['sha256'],
            pb_equivalence_evidence=str(HERE/'pb_flags/validation/comparison.json'),
            vdw_equivalence_evidence=str(HERE/'vdw/validation/large_audit_report.json'))
        core.atomic_json(baseline/'status.json',saved)
    pipeline.task.validate_baseline(baseline,baseline/'baseline.input',manifest['cap_mode'],protocol)
    return manifest,protocol


def donate_preparation(record,manifest,protocol):
    donor=record.get('preparation')
    if not donor:
        return
    workspace=Path(donor['workspace'])
    state=read(workspace/'status.json')
    gate=workspace/'pipeline_energy_ready.json'
    if state.get('stage')!='waiting_pipeline_energy' or state.get('pb_runs') or gate.exists():
        raise ValueError('Donor is not safely waiting before its first energy operation')
    destination=NEW/'batches/batch_0001'
    if (workspace/'job/job.input').read_bytes()!=(destination/'batch.input').read_bytes():
        raise ValueError('Prepared batch does not match continuation batch 1')
    if core.engine_alive(donor['pid'],workspace):
        os.kill(donor['pid'],signal.SIGTERM)
        deadline=time.monotonic()+20
        while core.engine_alive(donor['pid'],workspace) and time.monotonic()<deadline:
            time.sleep(.2)
        if core.engine_alive(donor['pid'],workspace):
            raise ValueError('Donor worker did not stop')
    root=core.workspace_root(NEW)
    target=root/'batch_0001_attempt_001'
    if target.exists():
        raise ValueError('Migrated workspace already exists; inspect before retrying')
    shutil.copytree(workspace,target,symlinks=True)
    scientific={str(p.relative_to(target)):digest(p) for folder in ('job','jobout','sinputfiles')
        for p in (target/folder).rglob('*') if p.is_file()}
    core.atomic_json(target/'gpu_migration.json',dict(donated_from=str(workspace),
        old_state=state,new_protocol=protocol,scientific_files=scientific,
        evidence=str(HERE/'pb_flags/validation/comparison.json')))
    state.update(stage='prepared_for_gpu',protocol_sha256=protocol,baseline=str(NEW/'baseline'),
                 workspace=str(target),donated_preparation=str(workspace))
    core.atomic_json(target/'status.json',state)
    params=core.task_arguments(manifest,argparse.Namespace(threads=8,no_contributions=False),NEW/'provean',protocol)
    core.atomic_json(destination/'worker_parameters.json',params)
    # Validate the complete frozen checkpoint before launching its worker.
    from resume_prepared import validate
    validate(target,destination/'batch.input',NEW/'baseline',argparse.Namespace(**params),state)
    command=[sys.executable,'-u',str(HERE/'resume_prepared.py'),'mutabind2_tool.scan.task',
        '--source',str(destination/'batch.input'),'--workspace',str(target),
        '--parameters',str(destination/'worker_parameters.json'),'--config',str(NEW/'config.snapshot.json'),
        '--baseline',str(NEW/'baseline')]
    with (destination/'run_001.log').open('w') as log:
        proc=subprocess.Popen(command,cwd=INSTALL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    core.atomic_json(destination/'active_process.json',dict(workspace=str(target),pid=proc.pid,
        command=command,started_utc=core.now()))


def run():
    install()
    transition=HERE/'transition.json'
    record=read(transition)
    marker=CELLS/'STOP_AFTER_BATCH'
    if read(CELLS/'batch_status.json').get('stage')=='running':
        if marker.exists() and not record:
            raise ValueError('A separate pause request is already pending')
        if not record:
            record=dict(requested_utc=core.now(),preparation=read(CELLS/'pipeline_status.json').get('preparation'))
            core.atomic_json(transition,record)
            marker.write_text('Validated GPU lists and MSTEP transition\n')
        print(json.dumps(dict(stage='waiting_for_batch_boundary',preparation=record['preparation'])),flush=True)
        while read(CELLS/'batch_status.json').get('stage')=='running':
            time.sleep(2)
    manifest,protocol=prepare()
    files={str(p):digest(p) for p in (Path(__file__).resolve(),HERE/'resume_prepared.py')}
    extension=NEW/'gpu_continuation_deployment.json'
    if extension.exists() and read(extension)['files']!=files:
        raise ValueError('Continuation scripts changed; restore their recorded versions')
    core.atomic_json(extension,dict(files=files,protocol_sha256=protocol))
    if marker.exists() and marker.read_text()=='Validated GPU lists and MSTEP transition\n':
        marker.write_text('Superseded by validated GPU continuation\n')
    if not (NEW/'pipeline_deployment.json').exists():
        donate_preparation(record,manifest,protocol)
    args=argparse.Namespace(output=NEW,combined=None,handoff=False,threads=8,no_contributions=False,
                            keep_workspace=False,poll=2)
    coordinator=pipeline.Pipeline(args)
    coordinator.merge=sys.modules[__name__]
    record.update(activated_utc=core.now(),new_scan=str(NEW),protocol_sha256=protocol,
                  frozen_wt_preserved=True)
    core.atomic_json(transition,record)
    core.atomic_json(COMBINED/'continuation_gpu.json',record)
    coordinator.run()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('run','status'))
    if parser.parse_args().action=='status':
        print(json.dumps(collect()[1],ensure_ascii=False,indent=2))
    else:
        with (HERE/'coordinator.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            run()
