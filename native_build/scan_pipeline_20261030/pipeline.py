"""Prepare one batch ahead while retaining the original scientific protocol."""
import argparse
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
INSTALL = HERE.parents[1]
sys.path.insert(0, str(INSTALL))
from mutabind2_tool.core.runtime import load_config
from mutabind2_tool.scan import plan, scheduler as core, task
from mutabind2_tool.scan.optimize import energy


def read(path):
    return json.loads(path.read_text()) if path.exists() else {}


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


class Pipeline:
    def __init__(self, args):
        self.args = args
        self.output = args.output.resolve()
        self.manifest = plan.load(self.output)
        self.root = core.workspace_root(self.output)
        self.baseline = self.output / 'baseline'
        self.batches = read(self.output / 'batches.json')['batches']
        self.requests = core.read_table(self.output / 'plan/mutations.input')
        self.protocol = read(self.output / 'protocol.json')['sha256']
        config = load_config(self.output / 'config.snapshot.json')
        _, actual = core.protocol(self.manifest, config, self.output / 'provean')
        if actual != self.protocol:
            raise ValueError('Scientific code, tools or input changed; cannot continue this scan')
        self.parameters = core.task_arguments(self.manifest, args, self.output / 'provean', self.protocol)
        task.validate_baseline(self.baseline, self.baseline / 'baseline.input',
                               self.manifest['cap_mode'], self.protocol)
        self.worker = HERE / 'pipeline_worker.py'
        self.marker = self.output / 'STOP_AFTER_BATCH'
        self.marker_text = 'Pipeline handoff ' + str(os.getpid()) + '\n'
        self.started = time.monotonic()
        self.children = {}
        self.handles = {}
        self.owned = False
        self.merge = None
        self.last_merge = 0
        self.minimum_available = None
        self.merged_lock = None
        if args.combined:
            spec = importlib.util.spec_from_file_location('scan_merge', INSTALL / 'native_build/pb_reen_warp_20261030/continue_scan.py')
            self.merge = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(self.merge)
            if self.merge.NEW != self.output or self.merge.COMBINED != args.combined.resolve():
                raise ValueError('Combined result coordinator does not match this scan')
        record = self.output / 'pipeline_deployment.json'
        hashes = {str(p): digest(p) for p in [Path(__file__).resolve(), self.worker]}
        previous = read(record)
        if previous and previous['files'] != hashes:
            raise ValueError('Pipeline extension changed; restore its recorded version before resuming')
        core.atomic_json(record, {'files': hashes, 'protocol_sha256': self.protocol,
                                 'scientific_package_unchanged': True, 'depth': 2,
                                 'energy_concurrency': 1, 'preparation_admission_gib': 30})

    def event(self, name, **details):
        with (self.output / 'pipeline_events.jsonl').open('a') as stream:
            stream.write(json.dumps(dict(event=name, utc=core.now(), **details)) + '\n')

    def expected(self, batch):
        rows = [r for r in self.requests if r['Result_Id'] in batch['result_ids']]
        if core.read_table(self.destination(batch) / 'batch.input') != rows:
            raise ValueError('Batch input differs from the immutable plan')
        return rows

    def destination(self, batch):
        return self.output / 'batches' / f"batch_{batch['index']:04d}"

    def pending(self):
        result = []
        for batch in self.batches:
            destination = self.destination(batch)
            if read(destination / 'status.json').get('stage') == 'complete':
                core._validate_batch(destination, self.expected(batch), self.protocol)
            else:
                result.append(batch)
        return result

    def start(self, batch):
        destination = self.destination(batch)
        self.expected(batch)
        launch_path = destination / 'active_process.json'
        launch = read(launch_path)
        workspace = Path(launch['workspace']) if launch.get('workspace') else None
        if workspace:
            core.owned_workspace(workspace, self.root)
            matches = [int(p.name) for p in Path('/proc').iterdir()
                       if p.name.isdigit() and core.engine_alive(int(p.name), workspace)]
            if len(matches) > 1:
                raise ValueError('Multiple workers own the same batch')
            if matches:
                launch['pid'] = matches[0]
                core.atomic_json(launch_path, launch)
            saved = read(workspace / 'status.json')
            if matches or saved.get('stage') == 'complete':
                if read(destination / 'worker_parameters.json') != self.parameters:
                    raise ValueError('Recovered worker parameters differ from this scan')
                if saved and saved.get('protocol_sha256') != self.protocol:
                    raise ValueError('Recovered worker protocol differs from this scan')
                handle = dict(batch=batch, workspace=workspace, pid=launch['pid'])
                self.handles[launch['pid']] = handle
                return handle
        label = f"batch_{batch['index']:04d}"
        attempts = [int(p.name.rsplit('_', 1)[1]) for p in self.root.glob(label + '_attempt_*')]
        workspace = self.root / f'{label}_attempt_{max(attempts, default=0)+1:03d}'
        core.owned_workspace(workspace, self.root)
        parameters = destination / 'worker_parameters.json'
        core.atomic_json(parameters, self.parameters)
        command = [sys.executable, '-u', str(self.worker), 'mutabind2_tool.scan.task',
                   '--source', str(destination / 'batch.input'), '--workspace', str(workspace),
                   '--parameters', str(parameters), '--config', str(self.output / 'config.snapshot.json'),
                   '--baseline', str(self.baseline)]
        launch = dict(workspace=str(workspace), pid=None, command=command, started_utc=core.now())
        core.atomic_json(launch_path, launch)
        with (destination / f'run_{max(attempts, default=0)+1:03d}.log').open('w') as log:
            proc = subprocess.Popen(command, cwd=INSTALL, stdin=subprocess.DEVNULL,
                                    stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                                    env=dict(os.environ, PYTHONPATH=str(INSTALL) + os.pathsep + os.environ.get('PYTHONPATH', '')))
        self.children[proc.pid] = proc
        launch['pid'] = proc.pid
        core.atomic_json(launch_path, launch)
        self.event('preparation_started', batch=batch['index'], pid=proc.pid, workspace=str(workspace))
        handle = dict(batch=batch, workspace=workspace, pid=proc.pid)
        self.handles[proc.pid] = handle
        return handle

    def alive(self, handle):
        child = self.children.get(handle['pid'])
        if child is not None and child.poll() is not None:
            return False
        return core.engine_alive(handle['pid'], handle['workspace'])

    def release(self, handle):
        while not (handle['workspace'] / 'status.json').exists():
            if not self.alive(handle):
                raise RuntimeError('Batch failed before initialization')
            time.sleep(0.1)
        gate = handle['workspace'] / 'pipeline_energy_ready.json'
        if not gate.exists():
            core.atomic_json(gate, {'protocol_sha256': self.protocol, 'released_utc': core.now()})
            self.event('energy_released', batch=handle['batch']['index'])

    def status(self, current=None, future=None, phase='running'):
        available = energy._available()
        self.minimum_available = min(self.minimum_available or available, available)
        def detail(handle):
            if handle is None:
                return None
            state = read(handle['workspace'] / 'status.json')
            return dict(batch=handle['batch']['index'], pid=handle['pid'],
                        workspace=str(handle['workspace']), stage=state.get('stage', 'starting'),
                        mutation=state.get('pb_current_mutation'),
                        energy_complete=sum(k.startswith('mutant_') and v.get('outcome') == 'complete'
                                            for k, v in state.get('pb_runs', {}).items()))
        core.atomic_json(self.output / 'pipeline_status.json', dict(
            stage=phase, pid=os.getpid(), updated_utc=core.now(), energy=detail(current),
            preparation=detail(future), available_gib=round(available / 1024**3, 3),
            lowest_available_gib=round(self.minimum_available / 1024**3, 3),
            protocol_sha256=self.protocol, depth=2, energy_concurrency=1))
        if self.owned and current:
            native = read(current['workspace'] / 'status.json')
            core.update(self.output, stage='running', pid=os.getpid(), error='',
                        batch=current['batch']['index'], batch_mutations=len(current['batch']['result_ids']),
                        workspace=str(current['workspace']), active_pid=current['pid'],
                        native_stage=native.get('stage'), current_mutation=native.get('pb_current_mutation'))
        if self.owned and self.merge and time.monotonic() - self.last_merge >= 10:
            self.merge.aggregate()
            self.last_merge = time.monotonic()

    def maybe_prepare(self, current, batch):
        state = read(current['workspace'] / 'status.json')
        if (state.get('pb_runs') or state.get('stage') == 'waiting_pipeline_energy') and energy._available() >= 30 * 1024**3:
            return self.start(batch)
        return None

    def handoff(self, lock):
        state = read(self.output / 'batch_status.json')
        if not self.args.handoff:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            return
        handoff_record = self.output / 'pipeline_handoff.json'
        if self.marker.exists():
            previous = read(handoff_record)
            if not previous or self.marker.read_text() != previous['marker_text']:
                raise ValueError('A separate stop request already exists; preserve the requested pause')
            self.marker_text = previous['marker_text']
        else:
            if state.get('stage') != 'running':
                raise ValueError('Handoff requires a running original scheduler')
            with self.marker.open('x') as stream:
                stream.write(self.marker_text)
            core.atomic_json(handoff_record, dict(marker_text=self.marker_text, requested_utc=core.now()))
        current_batch = next(b for b in self.batches if b['index'] == state['batch'])
        current = self.start(current_batch) if state.get('stage') != 'paused' else None
        following = [b for b in self.pending() if b['index'] > current_batch['index']]
        future = None
        self.event('handoff_requested', current_batch=current_batch['index'])
        while True:
            state = read(self.output / 'batch_status.json')
            if state.get('stage') == 'paused':
                break
            if state.get('stage') == 'failed' or time.monotonic() - self.started > 7200:
                raise RuntimeError('Original scheduler did not complete its batch cleanly')
            if future is None and following and current:
                future = self.maybe_prepare(current, following[0])
            if future and not self.alive(future) and read(future['workspace'] / 'status.json').get('stage') != 'complete':
                raise RuntimeError('Prefetched batch preparation failed')
            self.status(current, future, 'handoff')
            time.sleep(self.args.poll)
        fcntl.flock(lock, fcntl.LOCK_EX)
        if self.merge:
            self.merged_lock = (self.merge.COMBINED / 'continuation.lock').open('a')
            fcntl.flock(self.merged_lock, fcntl.LOCK_EX)
        if self.marker.read_text() != self.marker_text:
            raise ValueError('Stop marker changed during handoff; preserve the requested pause')
        self.marker.unlink()
        self.event('handoff_complete', original_completed=state.get('completed'))

    def run_owned(self):
        self.owned = True
        core.update(self.output, stage='running', pid=os.getpid(), error='', scheduler='pipeline',
                    pipeline_status=str(self.output / 'pipeline_status.json'))
        core.aggregate(self.output, self.manifest)
        pending = self.pending()
        current = None
        for position, batch in enumerate(pending):
            if self.marker.exists():
                break
            current = self.start(batch)
            self.release(current)
            future = None
            while self.alive(current):
                if future is None and position + 1 < len(pending) and not self.marker.exists():
                    future = self.maybe_prepare(current, pending[position + 1])
                if future and not self.alive(future) and read(future['workspace'] / 'status.json').get('stage') != 'complete':
                    raise RuntimeError('Prefetched batch preparation failed')
                self.status(current, future)
                time.sleep(self.args.poll)
            saved = read(current['workspace'] / 'status.json')
            if saved.get('stage') != 'complete' or saved.get('protocol_sha256') != self.protocol:
                raise RuntimeError(f"Batch {batch['index']} failed: {saved.get('error', 'inspect its run log')}")
            core._archive_batch(current['workspace'], self.destination(batch), self.root,
                                self.expected(batch), self.protocol, self.args.keep_workspace)
            self.event('batch_archived', batch=batch['index'])
            core.aggregate(self.output, self.manifest)
            if self.merge:
                self.merge.aggregate()
        remaining = self.pending()
        phase = 'paused' if remaining else 'complete'
        core.update(self.output, stage=phase, active_pid=None, native_stage='',
                    completed_utc=core.now() if not remaining else None)
        self.status(phase=phase)
        if self.merge:
            self.merge.aggregate()

    def stop_children(self):
        targets = dict(self.children)
        if self.owned:
            targets.update({pid: self.children.get(pid) for pid in self.handles if pid is not None})
        def running(pid, child):
            if child is not None:
                return child.poll() is None
            return core.engine_alive(pid, self.handles[pid]['workspace'])
        for pid, child in targets.items():
            if running(pid, child):
                try:
                    os.killpg(pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
        deadline = time.monotonic() + 10
        while any(running(pid, child) for pid, child in targets.items()) and time.monotonic() < deadline:
            time.sleep(0.1)
        for pid, child in targets.items():
            if running(pid, child):
                try:
                    os.killpg(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                if child is not None:
                    child.wait()

    def run(self):
        with (self.output / 'pipeline.lock').open('a') as own, (self.output / 'scan.lock').open('a') as lock:
            fcntl.flock(own, fcntl.LOCK_EX | fcntl.LOCK_NB)
            try:
                self.handoff(lock)
                if self.merge and self.merged_lock is None:
                    self.merged_lock = (self.merge.COMBINED / 'continuation.lock').open('a')
                    fcntl.flock(self.merged_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                self.run_owned()
            except BaseException as exc:
                self.stop_children()
                if self.marker.exists() and self.marker.read_text() == self.marker_text:
                    self.marker.unlink()
                self.status(phase='failed')
                if self.owned:
                    core.update(self.output, stage='failed', error=str(exc))
                    if self.merge:
                        self.merge.aggregate()
                self.event('failed', error=str(exc))
                raise
            finally:
                if self.merged_lock:
                    self.merged_lock.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--combined', type=Path)
    parser.add_argument('--handoff', action='store_true')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--no-contributions', action='store_true')
    parser.add_argument('--keep-workspace', action='store_true')
    parser.add_argument('--poll', type=float, default=2)
    args = parser.parse_args()
    if args.poll <= 0 or args.threads < 1:
        parser.error('poll and threads must be positive')
    def interrupted(signum, frame):
        raise InterruptedError(f'Pipeline interrupted by signal {signum}')
    signal.signal(signal.SIGTERM, interrupted)
    Pipeline(args).run()


if __name__ == '__main__':
    main()
