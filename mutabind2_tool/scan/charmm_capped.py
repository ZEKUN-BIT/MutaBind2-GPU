#!/usr/bin/env python3
"""Supply the matching neutral N-terminal patch to a CHARMM setup deck."""

import os
from pathlib import Path
import sys

TOOL_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(TOOL_ROOT))
from mutabind2_tool.core.runtime import load_config  # noqa: E402


def main():
    args = sys.argv[1:]
    parameters = dict(arg.split("=", 1) for arg in args if "=" in arg)
    if "I" in parameters:
        source = Path(parameters["path"]) / f"{parameters['protname']}_ch{parameters['I']}.pdb"
        first = next((line for line in source.read_text().splitlines()
                      if line.startswith("ATOM  ")), None)
        if first is None:
            raise ValueError(f"No ATOM record in {source}")
        args.append("npatch=" + ("ACP" if first[17:20] == "PRO" else "ACE"))
    config = load_config(os.environ.get("MUTABIND2_CONFIG"))
    charmm = config.get("charmm")
    if not charmm:
        raise ValueError("charmm is not configured")
    os.execv(charmm, [charmm, *args])


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as exc:
        print(f"CHARMM cap wrapper error: {exc}", file=sys.stderr)
        raise SystemExit(2)
