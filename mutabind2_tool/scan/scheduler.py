"""Serial GPU energy scheduling, immutable plans, verified checkpoints and recovery."""
from datetime import datetime, timezone
import argparse
import csv
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import tarfile
import time

from ..core.runtime import ALLOWED_KEYS, PACKAGE_ROOT, TOOL_ROOT, load_config
from . import plan as plan_module
from . import task as task_module
from .provean import prepare_scan_cache, validated_scores
from .optimize import energy


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    return energy._hash(path)


def atomic_json(path, data):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')
    temporary.replace(path)


def read_table(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write_table(path, fields, rows):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def parse_batches(value):
    if value is None:
        return None
    if not re.fullmatch(r'[1-9][0-9]*(?:,[1-9][0-9]*)*', value):
        raise ValueError('--only-batches must be comma-separated positive integers')
    return {int(x) for x in value.split(',')}


def update(output, **changes):
    path = output / 'batch_status.json'
    data = json.loads(path.read_text()) if path.exists() else {}
    data.update(changes, updated_utc=now())
    atomic_json(path, data)
    return data


def workspace_root(output, explicit=None):
    base = Path(explicit or os.environ.get('MUTABIND2_WORKSPACE_ROOT', TOOL_ROOT / 'workspaces')).expanduser().resolve()
    if not re.fullmatch(r'/[a-z0-9_./-]+', str(base)):
        raise ValueError('Native scratch path must use lowercase ASCII without spaces; use --workspace-root')
    scan_id = hashlib.sha256(str(Path(output).resolve()).encode()).hexdigest()[:16]
    root = base / scan_id
    root.mkdir(parents=True, exist_ok=True)
    owner = {'scan': str(Path(output).resolve()), 'id': scan_id}
    marker = root / 'owner.json'
    if marker.exists() and json.loads(marker.read_text()) != owner:
        raise ValueError('Scratch directory belongs to another scan')
    atomic_json(marker, owner)
    return root


def owned_workspace(workspace, root):
    if (workspace.is_symlink() or workspace.parent != root or
            not re.fullmatch(r'(?:baseline|batch_\d{4})_attempt_\d+', workspace.name)):
        raise ValueError(f'Refusing to remove an unowned workspace: {workspace}')


def task_arguments(manifest, args, cache, protocol_hash):
    identities = {rid: [chain, mutation] for rid, chain, mutation in manifest['mutation_identity']}
    return {'cap_mode': manifest['cap_mode'], 'mutated_partner': manifest['mutated_partner'],
            'surface_backend': manifest['surface_backend'], 'provean_mode': 'cached',
            'provean_cache': str(cache), 'pb_executable': manifest['pb_executable'],
            'vdw_executable': manifest['vdw_executable'], 'pb_threads': 1,
            'foldx_workers': manifest['foldx_workers'], 'threads': args.threads,
            'no_contributions': args.no_contributions, 'identities': identities,
            'protocol_sha256': protocol_hash}


def protocol(manifest, config, cache):
    required = ['vmd', 'charmm', 'namd', 'dssp', 'foldx', 'rscript']
    for key in required:
        value = config.get(key)
        if not value or not Path(value).is_file() or not os.access(value, os.X_OK):
            raise ValueError(f'Required tool is not executable: {key}={value}')
    pb = energy._runtime_identity(manifest['pb_executable'])
    if (manifest.get('backend', 'gpu') == 'gpu') != (pb['pb_backend'] == 'cuda'):
        raise ValueError('Selected PB executable does not match --backend gpu/cpu')
    vdw = energy._runtime_identity(manifest['vdw_executable']) if manifest.get('vdw_executable') else None
    data = {key: manifest[key] for key in ('cap_mode', 'mutated_partner', 'surface_backend', 'provean_mode')}
    data.update(pb=pb, vdw=vdw, charmm=energy._runtime_identity(config['charmm']),
                tools={key: {'path': config[key], 'sha256': sha256(config[key])} for key in required},
                resources={str(p.relative_to(PACKAGE_ROOT)):sha256(p)
                           for p in (PACKAGE_ROOT/'engines/resources/single').rglob('*') if p.is_file()},
                code={str(p.relative_to(PACKAGE_ROOT)):sha256(p) for p in PACKAGE_ROOT.rglob('*.py')},
                score_script=sha256(PACKAGE_ROOT/'core/score.R'),
                provean_scores=sha256(cache/'scores.tsv'),
                plan=manifest['plan']['prepared_pdb_sha256'])
    encoded = json.dumps(data, sort_keys=True).encode()
    return data, hashlib.sha256(encoded).hexdigest()


def engine_alive(pid, workspace):
    if not isinstance(pid, int) or pid <= 0:
        return False
    try:
        process = Path('/proc')/str(pid)
        argv = (process/'cmdline').read_bytes().decode().split('\0')
        state = (process/'stat').read_text().rsplit(')', 1)[1].split()[0]
        return (state not in {'Z', 'X'} and 'mutabind2_tool.scan.task' in argv and
                '--workspace' in argv and Path(argv[argv.index('--workspace')+1]) == workspace)
    except (OSError, ValueError, IndexError):
        return False


def execute_task(output, destination, label, source, root, parameters, baseline=None):
    launch_path = destination/'active_process.json'
    launch = json.loads(launch_path.read_text()) if launch_path.exists() else {}
    workspace = Path(launch['workspace']) if launch.get('workspace') else None
    child = None
    if workspace:
        owned_workspace(workspace, root)
        # Recover the narrow crash window between spawning and recording PID.
        matches = [int(p.name) for p in Path('/proc').iterdir()
                   if p.name.isdigit() and engine_alive(int(p.name), workspace)]
        if len(matches) > 1:
            raise RuntimeError('Multiple engines own the same scratch directory')
        if matches:
            launch['pid'] = matches[0]
            atomic_json(launch_path, launch)
        saved = json.loads((workspace/'status.json').read_text()) if (workspace/'status.json').exists() else {}
        if not engine_alive(launch.get('pid'), workspace) and saved.get('stage') != 'complete':
            workspace = None
    if workspace is None:
        attempts = list(root.glob(label+'_attempt_*'))
        number = max([int(p.name.rsplit('_',1)[1]) for p in attempts]+[0])+1
        workspace = root/f'{label}_attempt_{number:03d}'
        owned_workspace(workspace, root)
        args_file = destination/'worker_parameters.json'
        atomic_json(args_file, parameters)
        argv = [sys.executable, '-u', '-m', 'mutabind2_tool.scan.task',
                '--source', str(source), '--workspace', str(workspace),
                '--parameters', str(args_file), '--config', str(output/'config.snapshot.json')]
        argv += ['--baseline', str(baseline)] if baseline else ['--wt-only']
        launch = {'workspace':str(workspace), 'pid':None, 'command':argv, 'started_utc':now()}
        atomic_json(launch_path, launch)
        with (destination/f'run_{number:03d}.log').open('w') as log:
            env = dict(os.environ, PYTHONPATH=str(TOOL_ROOT)+os.pathsep+os.environ.get('PYTHONPATH',''))
            child = subprocess.Popen(argv, stdout=log, stderr=subprocess.STDOUT,
                                     start_new_session=True, env=env)
        launch['pid'] = child.pid
        atomic_json(launch_path, launch)
    update(output, stage='running', workspace=str(workspace), active_pid=launch['pid'], native_stage='starting')
    try:
        while (child.poll() is None if child else engine_alive(launch.get('pid'), workspace)):
            saved_path = workspace/'status.json'
            if saved_path.exists():
                saved = json.loads(saved_path.read_text())
                update(output, native_stage=saved.get('stage'), current_mutation=saved.get('pb_current_mutation'))
            time.sleep(2)
    except BaseException:
        if engine_alive(launch.get('pid'), workspace):
            os.kill(launch['pid'], signal.SIGTERM)
        raise
    saved = json.loads((workspace/'status.json').read_text()) if (workspace/'status.json').exists() else {}
    if saved.get('stage') != 'complete' or (child and child.returncode):
        raise RuntimeError(f"{label} failed: {saved.get('error', 'inspect the batch run log')}")
    if saved.get('protocol_sha256') != parameters['protocol_sha256']:
        raise ValueError('Recovered calculation protocol differs from this scan')
    return workspace, saved


def ensure_baseline(output, manifest, parameters, root):
    baseline = output/'baseline'
    baseline.mkdir(exist_ok=True)
    rows = read_table(output/'plan/mutations.input')
    source = baseline/'baseline.input'
    write_table(source, list(rows[0]), [rows[0]])
    shutil.copyfile(output/'plan/fragments.pdb', baseline/'fragments.pdb')
    if (baseline/'status.json').exists():
        task_module.validate_baseline(baseline, source, manifest['cap_mode'], parameters['protocol_sha256'])
        return baseline
    workspace, saved = execute_task(output, baseline, 'baseline', source, root, parameters)
    for name in ('job', 'jobout', 'sinputfiles', 'wt_cache'):
        shutil.copytree(workspace/name, baseline/name, dirs_exist_ok=True)
    if (workspace/'commands.log').exists():
        shutil.copyfile(workspace/'commands.log', baseline/'commands.log')
    atomic_json(baseline/'status.json', saved)
    task_module.validate_baseline(baseline, source, manifest['cap_mode'], parameters['protocol_sha256'])
    owned_workspace(workspace, root)
    shutil.rmtree(workspace)
    return baseline


def _validate_batch(destination, expected_rows, expected_protocol=None):
    state = json.loads((destination/'status.json').read_text())
    if state.get('stage') != 'complete':
        raise ValueError(f'Incomplete checkpoint: {destination}')
    if expected_protocol is not None and state.get('protocol_sha256') != expected_protocol:
        raise ValueError('Checkpoint calculation protocol differs from this scan')
    rows = read_table(destination/'results.tsv')
    wanted = {r['Result_Id']:r for r in expected_rows}
    ids = [r['Result_Id'] for r in rows]
    if len(wanted)!=len(expected_rows) or len(ids)!=len(set(ids)) or set(ids)!=set(wanted):
        raise ValueError('Checkpoint result IDs differ from plan')
    if set(state.get('result_ids', [])) != set(wanted):
        raise ValueError('Checkpoint status IDs differ from plan')
    identity = ('PDBfile','Partner1','Partner2','MutChain','Mutation_PDB','isPI')
    for r in rows:
        if any(r.get(k)!=wanted[r['Result_Id']].get(k) for k in identity):
            raise ValueError('Checkpoint mutation identity differs from plan')
        if not math.isfinite(float(r['SunDDG'])) or not r.get('Confidence'):
            raise ValueError('Invalid prediction or missing model labels')
    for name, expected in state.get('saved_sha256', {}).items():
        if sha256(destination/name) != expected:
            raise ValueError(f'Checkpoint checksum mismatch: {name}')
    return rows


def _archive_batch(workspace, destination, root, expected, protocol_hash, keep_workspace=False):
    _validate_batch(workspace, expected, protocol_hash)
    features = workspace/'job/job.input.cleaned.outdata'
    if not features.is_file():
        raise ValueError('Completed batch has no features')
    archive = destination/'scientific_intermediates.tar.gz'
    temporary = archive.with_suffix('.gz.tmp')
    with tarfile.open(temporary, 'w:gz', compresslevel=3) as tar:
        tar.add(workspace/'job', arcname='job')
        for path in (workspace/'jobout').iterdir():
            if path.is_file() and (path.suffix in {'.out','.fxout','.json','.log','.inp','.conf','.dssp'} or
                                  path.name.endswith(('_min.pdb','_min.crd','_charmm.psf','_charmm.crd'))):
                tar.add(path, arcname='jobout/'+path.name)
        for name in ('jobout/optimized_energy','jobout/optimized_surface','foldx_runs'):
            if (workspace/name).is_dir(): tar.add(workspace/name, arcname=name)
        for name in ('status.json','commands.log','results.tsv'):
            if (workspace/name).is_file():tar.add(workspace/name, arcname=name)
    with tarfile.open(temporary,'r:gz') as tar:
        for name in ('results.tsv','job/job.input.cleaned.outdata'):
            if hashlib.sha256(tar.extractfile(name).read()).hexdigest()!=sha256(workspace/name):
                raise ValueError('Archive verification failed')
    temporary.replace(archive)
    shutil.copyfile(workspace/'results.tsv', destination/'results.tsv')
    shutil.copyfile(features, destination/'features.tsv')
    saved = json.loads((workspace/'status.json').read_text())
    saved.update(workspace=str(workspace), workspace_removed=False,
                 saved_sha256={n:sha256(destination/n) for n in ('results.tsv','features.tsv',archive.name)})
    atomic_json(destination/'status.json', saved)
    _validate_batch(destination, expected, protocol_hash)
    if not keep_workspace:
        owned_workspace(workspace, root)
        if engine_alive(saved.get('pid'), workspace):raise ValueError('Refusing to remove an active workspace')
        shutil.rmtree(workspace)
        saved['workspace_removed']=True
        atomic_json(destination/'status.json',saved)


def aggregate(output, manifest):
    output = Path(output)
    requests = read_table(output/'plan/mutations.input')
    details = {r['result_id']:r for r in read_table(output/'plan/mutations.tsv')}
    annotations = {(r['original_chain'],r['original_number']):r for r in read_table(output/'plan/site_annotations.tsv')} if (output/'plan/site_annotations.tsv').exists() else {}
    batches = json.loads((output/'batches.json').read_text())['batches']
    protocol_hash = json.loads((output/'protocol.json').read_text())['sha256'] if (output/'protocol.json').exists() else None
    completed, batch_count = {}, 0
    for batch in batches:
        destination = output/'batches'/f"batch_{batch['index']:04d}"
        if not (destination/'status.json').exists():continue
        if json.loads((destination/'status.json').read_text()).get('stage')!='complete':continue
        expected = [r for r in requests if r['Result_Id'] in batch['result_ids']]
        for r in _validate_batch(destination,expected,protocol_hash):
            if r['Result_Id'] in completed:raise ValueError('Duplicate completed result')
            completed[r['Result_Id']]=r
        batch_count += 1
    scores = validated_scores(output/'provean') if (output/'provean/status.json').exists() else {}
    rows=[]
    for request in requests:
        rid=request['Result_Id'];detail=details[rid];result=completed.get(rid)
        chain,mutation=detail['chain'],detail['mutation']
        annotation=annotations.get((chain,mutation[1:-1]),{})
        value=float(result['SunDDG']) if result else None
        rows.append({'Result_Id':rid,'Chain':chain,'Mutation':mutation,'Status':'complete' if result else 'pending',
            'DeltaDeltaG_kcal_mol':result['SunDDG'] if result else '',
            'Predicted_binding_effect':('weaker' if value>0 else 'stronger' if value<0 else 'zero_rounded') if result else '',
            'Model_Interface':result.get('Interface','') if result else '',
            'Model_Confidence':result.get('Confidence','') if result else '',
            'PROVEAN':scores.get((chain,mutation),''),'PROVEAN_mode':manifest['provean_mode'],
            'Truncated_terminus':annotation.get('is_truncated_terminus',''),
            'Near_truncated_terminus_10A':annotation.get('near_truncated_terminus_10A',''),
            'Partner_min_heavy_distance_A':annotation.get('partner_min_heavy_distance_A',''),
            'Structure_note':manifest['fragment_mode']+'_fragments_'+manifest['cap_mode']+'_termini'})
    write_table(output/'predictions.tsv',list(rows[0]),rows)
    write_table(output/'completed_predictions.tsv',list(rows[0]),[r for r in rows if r['Status']=='complete'])
    update(output,total=len(rows),completed=len(completed),pending=len(rows)-len(completed),
           batches=len(batches),batches_complete=batch_count)
    return rows


def print_status(output):
    output=Path(output).resolve()
    manifest=plan_module.load(output)
    path=output/'batch_status.json'
    data=json.loads(path.read_text()) if path.exists() else {'stage':'planned','total':manifest['plan']['mutations']}
    workspace=Path(data['workspace']) if data.get('workspace') else None
    if data.get('stage')=='running' and workspace and (workspace/'status.json').exists():
        native=json.loads((workspace/'status.json').read_text())
        data.update(native_stage=native.get('stage'),current_mutation=native.get('pb_current_mutation'),
                    current_energies_complete=sum(k.startswith('mutant_') and v.get('outcome')=='complete' for k,v in native.get('pb_runs',{}).items()))
    print(json.dumps(data,ensure_ascii=False,indent=2))
    return 0


def _run_locked(args,output):
    if (output/'scan.json').exists():
        if not args.resume:raise ValueError('Existing scan: use --resume --output PATH')
        if any(getattr(args,k,None) for k in ('pdb','partner1','partner2','sites','mutations','fragments','provean_cache','pb_executable','vdw_executable')):
            raise ValueError('Resume uses the saved immutable plan; omit new structure/protocol arguments')
        manifest=plan_module.load(output)
    else:
        if args.resume:raise ValueError('Cannot resume without scan.json')
        manifest,_=plan_module.build(args)
    if args.plan_only:
        print(json.dumps(plan_module.summary(manifest),ensure_ascii=False,indent=2));return 0
    if args.threads<1 or (args.max_batches is not None and args.max_batches<1):
        raise ValueError('Thread and batch limits must be positive')
    batches=json.loads((output/'batches.json').read_text())['batches']
    selected=parse_batches(args.only_batches)
    if selected is not None and not selected <= {b['index'] for b in batches}:
        raise ValueError('Requested batch index is outside the plan')
    update(output,stage='starting',pid=os.getpid(),error='',started_utc=now())
    aggregate(output,manifest)
    if (output/'STOP_AFTER_BATCH').exists():
        update(output,stage='paused',active_pid=None);return 0
    config=load_config(args.config or manifest.get('config'))
    atomic_json(output/'config.snapshot.json',{k:config.get(k) for k in ALLOWED_KEYS})
    config=load_config(output/'config.snapshot.json')
    root=workspace_root(output,args.workspace_root)
    update(output,stage='conservation')
    cache=prepare_scan_cache(output,manifest,config,args.threads)
    identity,protocol_hash=protocol(manifest,config,cache)
    protocol_path=output/'protocol.json'
    if protocol_path.exists() and json.loads(protocol_path.read_text())['sha256']!=protocol_hash:
        raise ValueError('Tools, model, code or conservation data changed; use a new scan directory')
    atomic_json(protocol_path,{'sha256':protocol_hash,'identity':identity})
    parameters=task_arguments(manifest,args,cache,protocol_hash)
    baseline=ensure_baseline(output,manifest,parameters,root)
    if args.baseline_only:
        update(output,stage='baseline_complete',active_pid=None);return 0
    started=0
    requests=read_table(output/'plan/mutations.input')
    for batch in batches:
        index=batch['index']
        if selected is not None and index not in selected:continue
        destination=output/'batches'/f'batch_{index:04d}'
        expected=[r for r in requests if r['Result_Id'] in batch['result_ids']]
        if read_table(destination/'batch.input')!=expected:raise ValueError('Batch input differs from immutable plan')
        if (destination/'status.json').exists() and json.loads((destination/'status.json').read_text()).get('stage')=='complete':
            _validate_batch(destination,expected,protocol_hash);continue
        if (output/'STOP_AFTER_BATCH').exists():break
        update(output,batch=index,batch_mutations=len(expected))
        workspace,_=execute_task(output,destination,f'batch_{index:04d}',destination/'batch.input',root,parameters,baseline)
        _archive_batch(workspace,destination,root,expected,protocol_hash,args.keep_workspace)
        aggregate(output,manifest)
        started+=1
        if args.stop_after_batch or (args.max_batches is not None and started>=args.max_batches):break
    rows=aggregate(output,manifest)
    done=sum(r['Status']=='complete' for r in rows)
    stage='complete' if done==len(rows) else 'paused'
    update(output,stage=stage,active_pid=None,native_stage='',completed_utc=now() if stage=='complete' else None)
    print(json.dumps({'stage':stage,'completed':done,'total':len(rows),'results':str(output/'completed_predictions.tsv')},indent=2))
    return 0


def run_scan(args):
    output=Path(args.output).expanduser().resolve()
    output.mkdir(parents=True,exist_ok=True)
    with (output/'scan.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError('Another scheduler already owns this scan') from None
        previous=signal.getsignal(signal.SIGTERM)
        def interrupted(signum,frame):raise InterruptedError(f'Scheduler interrupted by signal {signum}')
        signal.signal(signal.SIGTERM,interrupted)
        try:return _run_locked(args,output)
        except BaseException as exc:
            status_path=output/'batch_status.json'
            state=json.loads(status_path.read_text()) if status_path.exists() else {}
            if state.get('pid') == os.getpid():
                update(output,stage='failed',error=str(exc))
            raise
        finally:signal.signal(signal.SIGTERM,previous)
