#!/usr/bin/env python3
"""Run a prepared scan with batch prefetch, memory waits and frozen recovery."""
import argparse
import importlib.util
from pathlib import Path
import signal
import sys

ROOT = Path(__file__).resolve().parents[1]
MEMORY = ROOT / 'native_build/pipeline_memory_20261030'
RESOURCE = ROOT / 'native_build/resource_admission_20261102'
sys.path[:0] = [str(MEMORY), str(ROOT / 'native_build/scan_pipeline_20261030'), str(ROOT)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='scan with a completed WT baseline')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--poll', type=float, default=2)
    parser.add_argument('--keep-workspace', action='store_true')
    parser.add_argument('--no-contributions', action='store_true')
    args = parser.parse_args()
    if args.threads < 1 or args.poll <= 0:
        parser.error('threads and poll must be positive')
    args.combined = None
    args.handoff = False
    spec = importlib.util.spec_from_file_location('memory_pipeline', MEMORY / 'continue_scan.py')
    controller = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(controller)
    controller.HERE = RESOURCE
    pipeline = controller.MemoryPipeline(args)
    pipeline.worker = RESOURCE / 'worker.py'
    files = [Path(__file__).resolve()]
    for directory, names in [(MEMORY, ('continue_scan.py', 'recover.py', 'memory_admission.py')),
                             (RESOURCE, ('worker.py', 'recover.py', 'resource_admission.py'))]:
        files.extend(directory / name for name in names)
    identity = {str(p.relative_to(ROOT)): controller.pipeline.digest(p) for p in files}
    record = pipeline.output / 'resource_pipeline_deployment.json'
    previous = controller.read(record)
    if previous and previous['files'] != identity:
        raise ValueError('Operational extension changed; use the recorded version or a new scan')
    controller.core.atomic_json(record, dict(files=identity, protocol_sha256=pipeline.protocol))
    def interrupted(signum, frame):
        raise InterruptedError(f'Pipeline interrupted by signal {signum}')
    signal.signal(signal.SIGTERM, interrupted)
    pipeline.run()


if __name__ == '__main__':
    main()
