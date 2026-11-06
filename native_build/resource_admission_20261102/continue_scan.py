"""Resume 256-374 with resource waits and validated frozen energy caches."""
import argparse
import importlib.util
from pathlib import Path
import signal
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MEMORY = ROOT / 'native_build/pipeline_memory_20261030'
OUTPUT = ROOT / 'projects/7mkj_256_374_20261102'
sys.path[:0] = [str(MEMORY), str(ROOT / 'native_build/scan_pipeline_20261030'), str(ROOT)]
spec = importlib.util.spec_from_file_location('memory_scheduler', MEMORY / 'continue_scan.py')
memory = importlib.util.module_from_spec(spec)
spec.loader.exec_module(memory)


def run():
    # Fail if a pinned deployment changed; do not rewrite historical evidence.
    previous = memory.read(OUTPUT / 'memory_admission_deployment.json')
    for name, expected in previous['files'].items():
        if memory.pipeline.digest(Path(name)) != expected:
            raise ValueError('Previous deployment changed: ' + name)
    args = argparse.Namespace(output=OUTPUT, combined=None, handoff=False,
        threads=8, no_contributions=False, keep_workspace=False, poll=2)
    coordinator = memory.MemoryPipeline(args)
    coordinator.worker = HERE / 'worker.py'
    # The existing start() uses HERE only to select its recovery executable.
    # Redirect this process's operational path, keeping its freeze validator.
    memory.HERE = HERE
    files = {str(HERE / name): memory.pipeline.digest(HERE / name)
             for name in ('continue_scan.py', 'worker.py', 'recover.py', 'resource_admission.py')}
    deployment = OUTPUT / 'resource_admission_deployment.json'
    recorded = memory.read(deployment)
    if recorded and recorded['files'] != files:
        raise ValueError('Resource admission extension changed after deployment')
    memory.core.atomic_json(deployment, dict(files=files,
        protocol_sha256=coordinator.protocol, deployed_utc=memory.core.now(),
        native_gpu_reserve_gib=2, native_host_admission_gib=30,
        gpu_prelaunch_admission_gib=5, gpu_memory_retries=3,
        scientific_package_unchanged=True))
    coordinator.run()


if __name__ == '__main__':
    def interrupted(signum, frame):
        raise InterruptedError(f'Pipeline interrupted by signal {signum}')
    signal.signal(signal.SIGTERM, interrupted)
    run()
