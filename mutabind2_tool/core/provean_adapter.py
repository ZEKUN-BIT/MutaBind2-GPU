"""Run PROVEAN 1.1.5 with configurable, modern BLAST/CD-HIT companions.

The native program builds some commands with a shell.  This adapter therefore
copies user inputs into a controlled temporary directory and refuses paths that
contain shell metacharacters.
"""

import argparse
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

TOOL_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(TOOL_ROOT))
from mutabind2_tool.core.runtime import load_config, configured_database  # noqa: E402

SAFE = re.compile(r"/[A-Za-z0-9_./+-]+")


def _safe(path, label):
    value = str(Path(path).resolve())
    if not SAFE.fullmatch(value):
        raise ValueError(f"{label} path must be absolute and contain no spaces or shell metacharacters: {value}")
    return value


def _tool(config, name, fallback):
    value = config.get(name)
    if value:
        return value
    return shutil.which(fallback)


def main(argv=None):
    parser = argparse.ArgumentParser(description="PROVEAN with accession and modern BLAST compatibility")
    parser.add_argument("-q", "--query", type=Path, required=True)
    parser.add_argument("-v", "--variation", type=Path, required=True)
    parser.add_argument("-d", "--database")
    parser.add_argument("--num_threads", type=int, default=4)
    parser.add_argument("--subject_sequences", type=Path)
    parser.add_argument("--supporting_set", type=Path)
    parser.add_argument("--save_supporting_set", type=Path)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)
    if args.num_threads < 1:
        parser.error("--num_threads must be positive")

    config = load_config(os.environ.get("MUTABIND2_CONFIG"))
    database = args.database or configured_database(config)
    if not database:
        parser.error("Set blast_db in the configuration or BLAST_DB_NR")
    if not Path(database).is_absolute():
        database = Path(database).resolve()
    database = _safe(database, "Database")

    provean = _tool(config, "provean_binary", "provean")
    if not provean:
        parser.error("provean_binary is not configured and no `provean` was found on PATH")
    provean = _safe(provean, "PROVEAN")
    psiblast = _tool(config, "psiblast", "psiblast")
    blastdbcmd = _tool(config, "blastdbcmd", "blastdbcmd")
    cdhit = _tool(config, "cdhit", "cd-hit")
    for label, value in (("psiblast", psiblast), ("blastdbcmd", blastdbcmd), ("cd-hit", cdhit)):
        if value:
            _safe(value, label)

    state = Path(database).parent / "status.json"
    if state.is_file() and json.loads(state.read_text()).get("stage") != "ready":
        parser.error("Database is incomplete; inspect " + str(state))

    with tempfile.TemporaryDirectory(prefix="mutabind2_provean_", dir="/tmp") as directory:
        work = Path(directory)
        command = [provean, "-d", database, "--num_threads", str(args.num_threads),
                   "--tmp_dir", str(work), "--quiet"]
        if psiblast:
            command += ["--psiblast", psiblast]
        if blastdbcmd:
            command += ["--blastdbcmd", blastdbcmd]
        if cdhit:
            command += ["--cdhit", cdhit]
        for source, option, name in [
            (args.query, "-q", "query.fasta"),
            (args.variation, "-v", "variants.txt"),
            (args.subject_sequences, "--subject_sequences", "subjects.fasta"),
            (args.supporting_set, "--supporting_set", "subjects.sss"),
        ]:
            if source:
                shutil.copyfile(source, work / name)
                command += [option, str(work / name)]
        if args.save_supporting_set:
            command += ["--save_supporting_set", str(work / "saved.sss")]

        environment = dict(os.environ)
        environment["BLAST_DB_NR"] = database
        result = subprocess.run(command, cwd=work, text=True, capture_output=True,
                                env=environment)
        sys.stderr.write(result.stderr)
        if result.returncode:
            sys.stderr.write(result.stdout)
            return result.returncode
        scores = result.stdout.partition("## PROVEAN scores ##")[2]
        rows = [line.split("\t") for line in scores.splitlines()
                if line and not line.startswith("#")]
        if not rows or any(len(row) != 2 or not math.isfinite(float(row[1])) for row in rows):
            raise ValueError("PROVEAN did not return finite scores")
        if args.save_supporting_set:
            shutil.copyfile(work / "saved.sss", args.save_supporting_set)
            fasta = work / "saved.sss.fasta"
            if fasta.exists():
                shutil.copyfile(fasta, str(args.save_supporting_set) + ".fasta")
        print(result.stdout, end="")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as exc:
        print("PROVEAN error: " + str(exc), file=sys.stderr)
        raise SystemExit(2)
