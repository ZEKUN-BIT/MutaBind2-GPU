"""Recover frozen geometry using the deployed validator and resource waits."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'native_build/pipeline_memory_20261030'),
               str(ROOT / 'native_build/scan_pipeline_20261030'), str(ROOT)]
import recover as previous
from resource_admission import install

if __name__ == '__main__':
    previous.install()
    install()
    previous.recovery.validate_prepared = previous.validate_prepared
    previous.recovery.main()
