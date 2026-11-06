"""Operational resource waits around the unchanged native energy runner.

The original decks, executables, convergence checks and cache validation remain
in charge of every accepted result. Only a recognized CUDA memory admission
rejection or a pre-launch host-memory rejection can be retried.
"""
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import time

from mutabind2_tool.scan import task
from mutabind2_tool.scan.optimize import energy

GIB = 1024 ** 3
MIN_HOST = energy.ADDRESS_LIMIT + energy.MEMORY_RESERVE
MIN_GPU = 5 * GIB
MAX_GPU_RETRIES = 3
POLL = 2
HOST_REJECTION = 'Requires 30 GiB available memory before each native process'
GPU_REJECTION = re.compile(
    r'Insufficient GPU memory: need ([0-9.]+) GiB plus ([0-9.]+) GiB reserve, free ([0-9.]+) GiB')


def gpu_available():
    # This scan uses the sole GPU (device 0). nvidia-smi reports MiB.
    result = subprocess.run(['nvidia-smi', '--id=0', '--query-gpu=memory.free',
                             '--format=csv,noheader,nounits'],
                            capture_output=True, text=True, check=True, timeout=10)
    lines = result.stdout.strip().splitlines()
    if len(lines) != 1 or not lines[0].strip().isdigit():
        raise RuntimeError('Cannot determine free GPU memory for device 0')
    return int(lines[0]) * 1024 ** 2


def memory_rejection(exc, case, mode):
    if mode != 'pb' or str(exc) != 'pb: CHARMM exited with code 84':
        return None
    log = Path(case) / 'pb.log'
    if not log.exists():
        return None
    text = log.read_text(errors='replace')
    match = GPU_REJECTION.search(text)
    if not match or 'ERROR STOP 84' not in text:
        return None
    required, reserve, free = map(float, match.groups())
    if not all(math.isfinite(v) and v >= 0 for v in (required, reserve, free)):
        return None
    # Keep the native 2 GiB reserve and add admission headroom for rounding.
    return math.ceil((required + reserve + 0.25) * GIB)


def wait_resources(workspace, stem, mode, required_gpu, retries):
    workspace = Path(workspace)
    started = time.monotonic()
    state_path = workspace / 'status.json'
    previous = json.loads(state_path.read_text()) if state_path.exists() else {}
    previous_stage = previous.get('stage')
    waited = False
    while True:
        host = energy._available()
        gpu = gpu_available() if host >= MIN_HOST else None
        admitted = host >= MIN_HOST and gpu is not None and gpu >= required_gpu
        info = dict(stage='admitted' if admitted else 'waiting_native_resources',
                    updated_utc=task.now(), stem=stem, mode=mode,
                    host_minimum_bytes=MIN_HOST, host_available_bytes=host,
                    gpu_minimum_bytes=required_gpu, gpu_available_bytes=gpu,
                    gpu_retry_count=retries, wait_seconds=round(time.monotonic() - started, 3))
        task.atomic_json(workspace / 'native_resource_admission.json', info)
        if not admitted or waited:
            if state_path.exists():
                state = json.loads(state_path.read_text())
                state['stage'] = previous_stage if admitted else 'waiting_native_resources'
                state['native_resource_admission'] = info
                task.atomic_json(state_path, state)
        if admitted:
            return
        waited = True
        time.sleep(POLL)


def install():
    original = energy._run_process
    if getattr(original, '_resource_admission', False):
        return

    def admitted(executable, workspace, stem, case, mode, threads):
        required_gpu = MIN_GPU
        retries = 0
        while True:
            wait_resources(workspace, stem, mode, required_gpu, retries)
            try:
                return original(executable, workspace, stem, case, mode, threads)
            except RuntimeError as exc:
                # This exact host rejection occurs before Popen, so no native
                # output has to be discarded or interpreted.
                if str(exc) == HOST_REJECTION:
                    continue
                needed = memory_rejection(exc, case, mode)
                if needed is None or retries >= MAX_GPU_RETRIES:
                    raise
                retries += 1
                required_gpu = max(required_gpu, needed)
                saved_log = Path(case) / f'{mode}.gpu_oom.{time.time_ns()}.log'
                shutil.copy2(Path(case) / f'{mode}.log', saved_log)
                with (Path(workspace) / 'native_resource_retries.jsonl').open('a') as stream:
                    stream.write(json.dumps(dict(utc=task.now(), mode=mode, stem=stem,
                        retry=retries, error=str(exc), log=str(saved_log),
                        gpu_minimum_bytes=required_gpu)) + '\n')
                # Re-enter the same native runner on the identical deck. Its
                # normal result audit is still mandatory after it succeeds.
                time.sleep(POLL)

    admitted._resource_admission = True
    energy._run_process = admitted
