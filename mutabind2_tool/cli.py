"""Command-line interface for the general MutaBind2 tool."""

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from . import __version__
from .core.runtime import (PACKAGE_ROOT, RESOURCES, TOOLS, config_path, configured_database,
                      load_config, score_file, stage_tools)
from .core.structure import inspect_job, prepare_structure

ENGINES = PACKAGE_ROOT / "engines"


def input_path(job):
    job = Path(job).resolve()
    if job.is_file():
        return job
    expected = job / (job.name + ".input")
    if expected.is_file():
        return expected
    matches = list(job.glob("*.input"))
    if len(matches) != 1:
        raise ValueError("Specify a job directory with one .input table, or the .input file itself")
    return matches[0]


def _usable(path):
    return bool(path) and Path(path).is_file() and os.access(path, os.X_OK)


def dependency_status(config, scope="full", contributions=True):
    tools = load_config(config)
    required = ["rscript"] if scope == "score" else list(TOOLS)
    status = {}
    for name in required:
        path = tools.get(name)
        if name == "provean" and not path:
            path = tools.get("provean_binary")
        status[name] = {"ok": _usable(path), "path": path}
    if tools.get("rscript"):
        packages = ["randomForest"] + (["forestFloor"] if contributions else [])
        expression = "pkgs <- c(" + ",".join(json.dumps(x) for x in packages) + "); "
        expression += 'for(p in pkgs) cat(p, requireNamespace(p,quietly=TRUE), "\\n")'
        result = subprocess.run([tools["rscript"], "--vanilla", "-e", expression],
                                text=True, capture_output=True,
                                env=dict(os.environ, RGL_USE_NULL="TRUE"))
        for package in packages:
            status["R:" + package] = {
                "ok": result.returncode == 0 and f"{package} TRUE" in result.stdout}
    if scope != "score":
        database = configured_database(tools)
        status["PROVEAN database"] = {"ok": False, "path": database}
        blastdbcmd = tools.get("blastdbcmd") or shutil.which("blastdbcmd")
        state = Path(database).parent / "status.json" if database else None
        incomplete = state and state.is_file() and \
            json.loads(state.read_text()).get("stage") != "ready"
        if incomplete:
            status["PROVEAN database"]["detail"] = \
                "Managed database is incomplete; inspect " + str(state)
        elif database and blastdbcmd:
            result = subprocess.run([blastdbcmd, "-db", database, "-info"],
                                    text=True, capture_output=True)
            status["PROVEAN database"]["ok"] = result.returncode == 0
            if result.returncode:
                status["PROVEAN database"]["detail"] = result.stderr.strip()[-500:]
        for helper, key in (("psiblast", "psiblast"), ("blastdbcmd", "blastdbcmd"),
                            ("cd-hit", "cdhit")):
            executable = tools.get(key)
            status[helper] = {"ok": _usable(executable), "path": executable}
        status["provean_binary"] = {"ok": _usable(tools.get("provean_binary")),
                                    "path": tools.get("provean_binary")}
    if scope == "gpu":
        from .scan.optimize.energy import _runtime_identity
        for key in ("pb_executable", "vdw_executable"):
            try:
                identity = _runtime_identity(tools[key])
                status[key] = {"ok": key != "pb_executable" or identity["pb_backend"] == "cuda",
                               "path": tools[key]}
            except (OSError, ValueError, TypeError, subprocess.SubprocessError) as exc:
                status[key] = {"ok": False, "detail": str(exc)}
        try:
            gpu = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader"],
                                 text=True, capture_output=True, timeout=15)
            status["NVIDIA GPU"] = {"ok": gpu.returncode == 0, "detail": gpu.stdout.strip() or gpu.stderr.strip()}
        except (OSError, subprocess.SubprocessError) as exc:
            status["NVIDIA GPU"] = {"ok": False, "detail": str(exc)}
    return tools, status


def ensure_output(path):
    path = Path(path).resolve()
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise ValueError(f"Output directory must be new or empty: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def stage_job(source, workspace, mode):
    rows, _ = inspect_job(source, mode)
    directory = workspace / "job"
    directory.mkdir()
    staged = directory / "job.input"
    shutil.copyfile(source, staged)
    shutil.copyfile(source.parent / rows[0]["PDBfile"], directory / rows[0]["PDBfile"])
    return staged


def _engine_resources(mode):
    return RESOURCES / ("single" if mode == "single" else "multiple")


def prepare(args):
    source = input_path(args.job)
    inspect_job(source, args.mode)
    output = ensure_output(args.output)
    staged = stage_job(source, output, args.mode)
    rows = prepare_structure(staged, output / "jobout", args.mode)
    manifest = {"mode": args.mode, "source": str(source), "rows": len(rows),
                "input": str(staged), "stage": "structure_preparation_only"}
    (output / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    print(f"Prepared {len(rows)} substitutions: {output}")
    print("This stage validates and renumbers the structure; it does not calculate ΔΔG.")


def run(args):
    source = input_path(args.job)
    inspect_job(source, args.mode)
    tools, status = dependency_status(args.config, contributions=not args.no_contributions)
    missing = [name for name, item in status.items() if not item["ok"]]
    if missing:
        raise RuntimeError("Full structure prediction is not ready. Missing/unusable: "
                           + ", ".join(missing)
                           + ". Run ./mutabind2 doctor and configure config/tools.json.")
    output = ensure_output(args.output)
    workspace = Path(tempfile.mkdtemp(prefix="mutabind2_", dir="/tmp"))
    staged = stage_job(source, workspace, args.mode)
    resource = _engine_resources(args.mode)
    shutil.copytree(resource, workspace / ("sinputfiles" if args.mode == "single"
                                           else "minputfiles"))
    stage_tools(workspace, tools)
    script = ENGINES / ("single.py" if args.mode == "single" else "multiple.py")
    argv = [sys.executable, "-u", str(script), "-i", "job", "--workdir", str(workspace),
            "--workers", str(args.workers), "--threads", str(args.threads)]
    if args.no_contributions:
        argv.append("--no-contributions")
    manifest = {"version": __version__, "mode": args.mode, "source": str(source),
                "workspace": str(workspace), "tools": tools, "command": argv,
                "stage": "running", "blast_db": configured_database(tools)}
    (output / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    print(f"Running {args.mode} prediction; logs: {output / 'run.log'}", flush=True)
    with (output / "run.log").open("w") as log:
        env = dict(os.environ, BLAST_DB_NR=configured_database(tools) or "")
        result = subprocess.run(argv, cwd=workspace, stdout=log, stderr=subprocess.STDOUT,
                                env=env)
    shutil.copytree(staged.parent, output / "job", dirs_exist_ok=True)
    if (workspace / "jobout").exists():
        shutil.copytree(workspace / "jobout", output / "intermediates", dirs_exist_ok=True)
    if (workspace / "commands.log").exists():
        shutil.copyfile(workspace / "commands.log", output / "commands.log")
    result_file = output / "job/job.sunddg"
    manifest["stage"] = "complete" if result.returncode == 0 and result_file.is_file() else "failed"
    (output / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    if result.returncode:
        raise RuntimeError(f"Prediction failed; inspect {output / 'run.log'}. "
                           f"Workspace retained: {workspace}")
    if not result_file.is_file():
        raise RuntimeError(f"The engine exited without a prediction table; inspect {output / 'run.log'}")
    shutil.copyfile(result_file, output / "results.tsv")
    print(f"Prediction table: {output / 'results.tsv'}")


def add_scan(subparsers):
    from .scan.cli import add_parser
    add_parser(subparsers)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="General MutaBind2 tool for structure-based ΔΔG prediction")
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    doctor = subparsers.add_parser("doctor", help="check executables, R packages and database")
    doctor.add_argument("--scope", choices=["full", "score", "gpu"], default="gpu")
    doctor.add_argument("--config", type=Path, default=None)
    doctor.add_argument("--no-contributions", action="store_true")
    doctor.add_argument("--json", action="store_true")
    for command in ["validate", "prepare", "run"]:
        child = subparsers.add_parser(command)
        child.add_argument("--mode", required=True, choices=["single", "multiple"])
        child.add_argument("--job", required=True, type=Path, help="job folder or .input file")
        if command != "validate":
            child.add_argument("--output", required=True, type=Path)
        if command == "run":
            child.add_argument("--config", type=Path, default=None)
            child.add_argument("--workers", type=int, default=2)
            child.add_argument("--threads", type=int, default=4)
            child.add_argument("--no-contributions", action="store_true")
    score = subparsers.add_parser("score-features",
                                  help="score precomputed physical features with the real model")
    score.add_argument("--mode", required=True, choices=["single", "multiple"])
    score.add_argument("--features", required=True, type=Path)
    score.add_argument("--output", required=True, type=Path)
    score.add_argument("--config", type=Path, default=None)
    score.add_argument("--contributions", action="store_true")
    add_scan(subparsers)
    args = parser.parse_args(argv)
    try:
        if args.command == "doctor":
            _, status = dependency_status(args.config, args.scope, not args.no_contributions)
            if args.json:
                print(json.dumps(status, ensure_ascii=False, indent=2))
            else:
                for name, item in status.items():
                    print(f"{'OK     ' if item['ok'] else 'MISSING'} {name:20} "
                          f"{item.get('path') or ''}")
                    if item.get("detail"):
                        print(f"        {item['detail']}")
            return 0 if all(item["ok"] for item in status.values()) else 1
        if args.command == "validate":
            source = input_path(args.job)
            rows, _ = inspect_job(source, args.mode)
            print(f"Valid {args.mode} input: {len(rows)} substitutions; {source}")
        elif args.command == "prepare":
            prepare(args)
        elif args.command == "run":
            if args.workers < 1 or args.threads < 1:
                raise ValueError("workers and threads must be positive")
            run(args)
        elif args.command == "scan":
            from .scan.cli import dispatch
            return dispatch(args)
        else:
            tools = load_config(args.config)
            destination = args.output.resolve()
            if destination.exists():
                raise ValueError(f"Output already exists: {destination}")
            destination.parent.mkdir(parents=True, exist_ok=True)
            rows = score_file(args.mode, args.features, destination, tools["rscript"],
                              args.contributions, tools)
            print(f"Scored {len(rows)} rows: {destination}")
        return 0
    except (OSError, ValueError, RuntimeError, csv.Error) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
