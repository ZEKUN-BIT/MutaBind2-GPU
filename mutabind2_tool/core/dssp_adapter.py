"""Feed classic PDB records to old/new mkdssp and request classic DSSP output."""

import re
import subprocess
import sys
import tempfile
from pathlib import Path


def normalize_pdb(source):
    result = ["HEADER    MUTABIND2 PREPARED PROTEIN COMPLEX\n"]
    for line in Path(source).read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atom = re.sub(r"[^A-Za-z]", "", line[12:16])
            if not atom:
                raise ValueError("Cannot infer the element for a PDB atom")
            # Inputs contain standard protein residues; CA is alpha carbon.
            element = atom[0].upper()
            result.append(f"{line[:76]:<76}{element:>2}\n")
        elif line.startswith("TER"):
            result.append("TER\n")
    result.append("END\n")
    return "".join(result)


def main(argv=None):
    executable, source, output = (argv or sys.argv[1:])
    help_text = subprocess.run([executable, "--help"], text=True, capture_output=True)
    flags = ["--output-format=dssp"] if "--output-format" in help_text.stdout + help_text.stderr else []
    with tempfile.TemporaryDirectory(prefix="mutabind2_dssp_") as directory:
        prepared = Path(directory) / "input.pdb"
        prepared.write_text(normalize_pdb(source))
        result = subprocess.run([executable, *flags, str(prepared), str(Path(output).resolve())])
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
