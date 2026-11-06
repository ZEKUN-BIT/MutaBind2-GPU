"""Existing pipeline worker with host and GPU memory admission."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'native_build/pipeline_memory_20261030'),
               str(ROOT / 'native_build/scan_pipeline_20261030'), str(ROOT)]
import pipeline_worker
from memory_admission import install as install_surface_wait
from resource_admission import install

if __name__ == '__main__':
    install_surface_wait()
    install()
    pipeline_worker.main()
