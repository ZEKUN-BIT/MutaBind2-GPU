"""Original batch worker with resource admission waiting."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
INSTALL = HERE.parents[1]
sys.path[:0] = [str(HERE), str(INSTALL / 'native_build/scan_pipeline_20261030'), str(INSTALL)]
import pipeline_worker
from memory_admission import install

if __name__ == '__main__':
    install()
    pipeline_worker.main()
