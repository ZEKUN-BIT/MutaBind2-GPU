"""Resume the GPU scan with waiting memory admission and frozen recovery."""
import argparse
import fcntl
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
PIPE = INSTALL / 'native_build/scan_pipeline_20261030'
sys.path[:0] = [str(HERE), str(PIPE), str(INSTALL)]
import pipeline
from recover import freeze

spec = importlib.util.spec_from_file_location('gpu_continuation',
    INSTALL / 'native_build/gpu_pipeline_20261030/continue_scan.py')
legacy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(legacy)
core, read = pipeline.core, pipeline.read


class MemoryPipeline(pipeline.Pipeline):
    def start(self, batch):
        destination = self.destination(batch)
        launch = read(destination / 'active_process.json')
        workspace = Path(launch['workspace']) if launch.get('workspace') else None
        if workspace and workspace.exists():
            core.owned_workspace(workspace, self.root)
            saved = read(workspace / 'status.json')
            alive = any(p.name.isdigit() and core.engine_alive(int(p.name), workspace)
                        for p in Path('/proc').iterdir())
            if not alive and saved.get('stage') == 'failed' and 'energy_decks' in saved.get('timings_seconds', {}):
                self.expected(batch)
                parameters = destination / 'worker_parameters.json'
                if read(parameters) != self.parameters:
                    raise ValueError('Recovered worker parameters differ from this scan')
                freeze(workspace, destination / 'batch.input', self.baseline,
                       argparse.Namespace(**self.parameters), saved)
                attempts = list(destination.glob('memory_recovery_*.log'))
                command = [sys.executable, '-u', str(HERE / 'recover.py'), 'mutabind2_tool.scan.task',
                    '--source', str(destination / 'batch.input'), '--workspace', str(workspace),
                    '--parameters', str(parameters), '--config', str(self.output / 'config.snapshot.json'),
                    '--baseline', str(self.baseline)]
                with (destination / f'memory_recovery_{len(attempts)+1:03d}.log').open('x') as log:
                    proc = subprocess.Popen(command, cwd=INSTALL, stdin=subprocess.DEVNULL,
                        stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                        env=dict(os.environ, PYTHONPATH=str(INSTALL) + os.pathsep + os.environ.get('PYTHONPATH', '')))
                self.children[proc.pid] = proc
                launch.update(pid=proc.pid, command=command, recovery_started_utc=core.now())
                core.atomic_json(destination / 'active_process.json', launch)
                handle = dict(batch=batch, workspace=workspace, pid=proc.pid)
                self.handles[proc.pid] = handle
                self.event('frozen_recovery_started', batch=batch['index'], pid=proc.pid,
                           workspace=str(workspace))
                return handle
        return super().start(batch)

    def run_owned(self):
        self.owned = True
        core.update(self.output, stage='running', pid=os.getpid(), error='', scheduler='pipeline',
                    pipeline_status=str(self.output / 'pipeline_status.json'))
        core.aggregate(self.output, self.manifest)
        pending = self.pending()
        for position, batch in enumerate(pending):
            if self.marker.exists():
                break
            current = self.start(batch)
            self.release(current)
            future = None
            preparation_error = None
            while self.alive(current):
                if future is None and position + 1 < len(pending) and not self.marker.exists():
                    future = self.maybe_prepare(current, pending[position + 1])
                if future and not preparation_error and not self.alive(future):
                    saved = read(future['workspace'] / 'status.json')
                    if saved.get('stage') != 'complete':
                        preparation_error = (f"Batch {future['batch']['index']} preparation failed: "
                                             + saved.get('error', 'inspect its run log'))
                        self.event('preparation_failed_finishing_current', error=preparation_error,
                                   current_batch=batch['index'])
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
            if preparation_error:
                raise RuntimeError(preparation_error + f"; batch {batch['index']} completed and saved")
        remaining = self.pending()
        phase = 'paused' if remaining else 'complete'
        core.update(self.output, stage=phase, active_pid=None, native_stage='',
                    completed_utc=core.now() if not remaining else None)
        self.status(phase=phase)
        if self.merge:
            self.merge.aggregate()


def run():
    args = argparse.Namespace(output=legacy.NEW, combined=None, handoff=False, threads=8,
        no_contributions=False, keep_workspace=False, poll=2)
    coordinator = MemoryPipeline(args)
    coordinator.worker = HERE / 'worker.py'
    coordinator.merge = legacy
    files = {str(HERE / name): pipeline.digest(HERE / name)
             for name in ('continue_scan.py', 'worker.py', 'recover.py', 'memory_admission.py')}
    path = legacy.NEW / 'memory_admission_deployment.json'
    previous = read(path)
    if previous and previous['files'] != files:
        raise ValueError('Memory admission extension changed after deployment')
    core.atomic_json(path, dict(files=files, protocol_sha256=coordinator.protocol,
        original_pipeline=read(legacy.NEW / 'pipeline_deployment.json'), deployed_utc=core.now(),
        scientific_package_unchanged=True, minimum_available_gib=30,
        finish_current_batch_on_prefetch_failure=True))
    coordinator.run()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('run', 'status'))
    if parser.parse_args().action == 'status':
        print(json.dumps(legacy.collect()[1], ensure_ascii=False, indent=2))
    else:
        def interrupted(signum, frame):
            raise InterruptedError(f'Pipeline interrupted by signal {signum}')
        signal.signal(signal.SIGTERM, interrupted)
        with (legacy.HERE / 'coordinator.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            run()
