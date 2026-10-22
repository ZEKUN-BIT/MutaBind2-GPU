"""Guard independent FoldX calls; retain the engine's per-mutation scheduler."""
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import tempfile
import threading
import time

GIB = 1024**3
RESERVE = 10 * GIB
ADDRESS_LIMIT = 8 * GIB
# Two concurrent FoldX workers may each use up to ADDRESS_LIMIT, plus the
# runtime reserve.  The previous 30 GiB budget was copied from the PB stage and
# is higher than the steady-state level while a worker runs, so the next
# admission is refused as soon as the first worker starts (self-defeating).
MIN_AVAILABLE = RESERVE + ADDRESS_LIMIT * 2
TIMEOUT = 1800
ADMISSION_TIMEOUT = 1800
ADMISSION_POLL_SECONDS = 15
POLL_SECONDS = 0.2
# Shared by installations in this interpreter, including test engines.
_SLOTS = threading.BoundedSemaphore(2)
_LOG_LOCK = threading.Lock()


def available_memory():
    return int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                    if line.startswith('MemAvailable:'))) * 1024


def wait_for_admission(cancelled):
    """Wait until at least MIN_AVAILABLE bytes are available, up to
    ADMISSION_TIMEOUT.  Returns (admission_bytes, wait_seconds)."""
    started = time.monotonic()
    admission = available_memory()
    while admission < MIN_AVAILABLE:
        if cancelled.is_set():
            raise InterruptedError('FoldX engine has been cancelled')
        if time.monotonic() - started > ADMISSION_TIMEOUT:
            raise RuntimeError(
                f'FoldX admission failed: {admission/GIB:.1f} GiB available, '
                f'requires {MIN_AVAILABLE/GIB:.0f} GiB after waiting {ADMISSION_TIMEOUT} s')
        time.sleep(ADMISSION_POLL_SECONDS)
        admission = available_memory()
    return admission, time.monotonic() - started


def _stop_group(process):
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=1)
    except subprocess.TimeoutExpired:
        pass
    # Also kill descendants if their group leader already exited after TERM.
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait()


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def install(engine, workspace):
    """Wrap independent BuildModel calls; return a main-thread restore callback."""
    if getattr(engine, '_guarded_foldx_installed', False):
        raise ValueError('FoldX guard already installed on this engine')
    workspace = Path(workspace).resolve()
    original = engine.run_command
    cancelled = threading.Event()
    def interrupted(signum, frame):
        cancelled.set()
        raise InterruptedError(f'FoldX engine interrupted by signal {signum}')
    previous_handlers = {}
    prlimit = shutil.which('prlimit')
    if not prlimit:
        raise RuntimeError('prlimit is required for the FoldX address-space guard')

    def guarded(command):
        if not command.startswith('./foldx '):
            return original(command)
        argv = shlex.split(command)
        if (len(argv) != 4 or argv[:2] != ['./foldx', '--command=BuildModel']
                or not argv[2].startswith('--pdb=') or not argv[3].startswith('--mutant-file=')):
            raise ValueError('Unsupported internal FoldX command')
        pdb_name = argv[2].split('=', 1)[1]
        mutant_name = argv[3].split('=', 1)[1]
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*\.pdb', pdb_name):
            raise ValueError('FoldX PDB must be a local basename')
        stem = Path(pdb_name).stem
        if mutant_name != f'individual_list_{stem}.txt':
            raise ValueError('FoldX mutation list does not match PDB identity')
        for name in [pdb_name, mutant_name]:
            if not (workspace/name).is_file():
                raise ValueError(f'Missing FoldX input: {name}')
        with _SLOTS:
            if cancelled.is_set():
                raise InterruptedError('FoldX engine has been cancelled')
            admission, admission_wait_seconds = wait_for_admission(cancelled)
            records = workspace/'foldx_runs'
            records.mkdir(exist_ok=True)
            run = Path(tempfile.mkdtemp(prefix=stem+'-', dir=records))
            command_argv = [prlimit, f'--as={ADDRESS_LIMIT}', '--core=0', '--',
                            str(workspace/'foldx'), *argv[1:]]
            record = {'input_pdb': pdb_name, 'mutation_list': mutant_name,
                      'input_sha256': _sha(workspace/pdb_name),
                      'mutations_sha256': _sha(workspace/mutant_name),
                      'command': command_argv, 'address_limit_bytes': ADDRESS_LIMIT,
                      'admission_available_bytes': admission, 'minimum_available_bytes': MIN_AVAILABLE,
                      'admission_wait_seconds': round(admission_wait_seconds, 3),
                      'reserve_bytes': RESERVE, 'timeout_seconds': TIMEOUT,
                      'started_unix': time.time()}
            start = time.monotonic()
            process = None
            peak = 0
            lowest = admission
            error = None
            try:
                with (run/'native.log').open('w') as output:
                    process = subprocess.Popen(command_argv, cwd=workspace, stdout=output,
                                               stderr=subprocess.STDOUT, start_new_session=True)
                    record['pid'] = process.pid
                    while process.poll() is None:
                        if cancelled.is_set():
                            raise InterruptedError('FoldX engine has been cancelled')
                        try:
                            status = Path(f'/proc/{process.pid}/status').read_text()
                            hwm = next((line.split()[1] for line in status.splitlines()
                                        if line.startswith('VmHWM:')), '0')
                            peak = max(peak, int(hwm)*1024)
                        except FileNotFoundError:
                            pass
                        lowest = min(lowest, available_memory())
                        if lowest < RESERVE:
                            raise RuntimeError('FoldX stopped to preserve 10 GiB available system memory')
                        if time.monotonic()-start > TIMEOUT:
                            raise TimeoutError('FoldX exceeded its 1800-second timeout')
                        time.sleep(POLL_SECONDS)
                    if process.returncode:
                        raise RuntimeError(f'FoldX failed with exit status {process.returncode}; see {run}/native.log')
            except BaseException as exc:
                error = str(exc)
                if process is not None:
                    _stop_group(process)
                raise
            finally:
                record.update(elapsed_seconds=round(time.monotonic()-start, 3), peak_rss_bytes=peak,
                              lowest_available_bytes=lowest, error=error,
                              returncode=process.returncode if process else None)
                (run/'result.json').write_text(json.dumps(record, indent=2)+'\n')
                payload = shlex.join(argv)+'\n'
                if (run/'native.log').exists():
                    payload += (run/'native.log').read_text(errors='replace')
                with _LOG_LOCK, (workspace/'commands.log').open('a') as log:
                    log.write(payload)
            return 0

    if threading.current_thread() is threading.main_thread():
        for signum in [signal.SIGINT, signal.SIGTERM]:
            previous_handlers[signum] = signal.signal(signum, interrupted)
    engine.run_command = guarded
    engine._guarded_foldx_installed = True

    def restore():
        # Call after parallel_tasks has joined all workers, including on errors.
        for signum, handler in previous_handlers.items():
            signal.signal(signum, handler)
        engine.run_command = original
        engine._guarded_foldx_installed = False

    return restore
