"""Job staging helpers shared by the CLI and the scan workflow."""

import shutil
from pathlib import Path

from .structure import inspect_job


def stage_job(source, workspace, mode):
    rows, _ = inspect_job(source, mode)
    directory = Path(workspace) / "job"
    directory.mkdir()
    staged = directory / "job.input"
    shutil.copyfile(source, staged)
    shutil.copyfile(Path(source).parent / rows[0]["PDBfile"],
                    directory / rows[0]["PDBfile"])
    return staged
