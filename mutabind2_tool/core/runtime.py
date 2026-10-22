"""Configuration, checked subprocess execution and model scoring.

The code keeps the upstream behaviour; machine-specific locations belong in a
JSON configuration file rather than in the source tree.
"""

import argparse
import csv
import glob
import io
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
TOOL_ROOT = PACKAGE_ROOT.parent
CORE_DIR = Path(__file__).resolve().parent
RESOURCES = PACKAGE_ROOT / "engines" / "resources"

# Executables used by the calculation engines.  The value is the wrapper name
# created inside each staged workspace, not necessarily the upstream filename.
TOOLS = {"vmd": "vmd", "charmm": "charmm", "namd": "namd2", "dssp": "mkdssp",
         "provean": "provean.sh", "foldx": "foldx", "rscript": "Rscript"}
EXTRA_KEYS = {"blast_db", "cdhit", "psiblast", "blastdbcmd", "provean_binary",
              "pb_executable", "vdw_executable", "toolchain_root"}
ALLOWED_KEYS = set(TOOLS) | EXTRA_KEYS
ALLOWED_COMMANDS = {"cat", "grep", "cut", "sed", "cp", "mv", "rm", "mkdir", *TOOLS.values()}

CONFIG_ENV = "MUTABIND2_CONFIG"
CONFIG_CANDIDATES = (TOOL_ROOT / "config" / "tools.json", TOOL_ROOT / "tools.json")


def config_path(explicit=None):
    """Return the configuration file in use, or None when all defaults apply."""
    if explicit:
        return Path(explicit).expanduser().resolve()
    from_environment = os.environ.get(CONFIG_ENV)
    if from_environment:
        return Path(from_environment).expanduser().resolve()
    for candidate in CONFIG_CANDIDATES:
        if candidate.is_file():
            return candidate.resolve()
    return None


def _path_value(value, base, name):
    if value is None or value == "":
        return None
    if not isinstance(value, str):
        raise ValueError(f"Configuration value {name} must be a path string or null")
    expanded = Path(os.path.expandvars(value)).expanduser()
    if not expanded.is_absolute():
        expanded = base / expanded
    return str(expanded.resolve())


def _executable_value(value, base, name):
    if value is None or value == "":
        return None
    if not isinstance(value, str):
        raise ValueError(f"Configuration value {name} must be a path string or null")
    expanded = os.path.expandvars(value)
    # A bare command name is looked up on PATH; anything with a separator is a
    # path relative to the configuration file.
    if "/" not in expanded and "\\" not in expanded:
        return shutil.which(expanded)
    return _path_value(expanded, base, name)


def load_config(path=None):
    """Return a resolved configuration dictionary.

    Missing files are allowed so that ``doctor`` can explain what is absent.
    Relative paths are resolved against the configuration file directory.
    """
    resolved_path = config_path(path)
    raw = {}
    if resolved_path is not None:
        if not resolved_path.is_file():
            raise FileNotFoundError(f"Configuration file not found: {resolved_path}")
        raw = json.loads(resolved_path.read_text())
        if not isinstance(raw, dict):
            raise ValueError("Tool configuration must be a JSON object")
        unknown = set(raw) - ALLOWED_KEYS
        if unknown:
            raise ValueError("Unknown tool configuration keys: " + ", ".join(sorted(unknown)))
    base = resolved_path.parent if resolved_path is not None else TOOL_ROOT
    result = dict(TOOLS)
    result.update({key: None for key in EXTRA_KEYS})
    for name, value in raw.items():
        result[name] = value

    for name in TOOLS:
        environment = os.environ.get("MUTABIND2_" + name.upper())
        if environment:
            result[name] = environment
        result[name] = _executable_value(result[name], base, name)
    for name in ("cdhit", "psiblast", "blastdbcmd", "provean_binary",
                 "pb_executable", "vdw_executable"):
        environment = os.environ.get("MUTABIND2_" + name.upper())
        if environment:
            result[name] = environment
        result[name] = _executable_value(result[name], base, name)
    if os.environ.get("BLAST_DB_NR"):
        result["blast_db"] = os.environ["BLAST_DB_NR"]
    result["blast_db"] = _path_value(result["blast_db"], base, "blast_db")
    result["toolchain_root"] = _path_value(result["toolchain_root"], base, "toolchain_root")
    result["config_path"] = str(resolved_path) if resolved_path is not None else None
    result["config_dir"] = str(base)
    return result


def configured_database(config=None):
    """Return the BLAST database prefix, honouring ``BLAST_DB_NR`` first."""
    environment = os.environ.get("BLAST_DB_NR")
    if environment:
        return str(Path(os.path.expandvars(environment)).expanduser().resolve())
    if config is None:
        config = load_config()
    elif not isinstance(config, dict):
        config = load_config(config)
    value = config.get("blast_db")
    if not value:
        return None
    return str(Path(value).expanduser())


def checked(argv, **kwargs):
    result = subprocess.run([str(value) for value in argv], text=True, capture_output=True, **kwargs)
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}): {shlex.join(map(str, argv))}\n"
                           + result.stderr[-6000:] + result.stdout[-2000:])
    return result.stdout


def _natural_key(value):
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", value)]


def command_output(command):
    """Interpret only the pipelines/redirections used upstream, without a shell.

    Input names are validated before reaching this function. Shell expansion,
    command substitution, environment expansion and command separators are absent.
    """
    lexer = shlex.shlex(command, posix=True, punctuation_chars="|<>")
    lexer.whitespace_split = True
    tokens = list(lexer)
    stages, stage = [], []
    for token in tokens:
        if token == "|":
            stages.append(stage)
            stage = []
        else:
            stage.append(token)
    stages.append(stage)
    data = None
    for stage in stages:
        argv, input_name, output_name = [], None, None
        position = 0
        while position < len(stage):
            token = stage[position]
            if token in {"<", ">"}:
                if position + 1 >= len(stage):
                    raise ValueError("Missing redirection filename")
                if token == "<":
                    input_name = stage[position + 1]
                else:
                    output_name = stage[position + 1]
                position += 2
            elif token in {"||", "<<", ">>", "&", ";"}:
                raise ValueError("Unsupported command operator")
            else:
                argv.append(token)
                position += 1
        if not argv or Path(argv[0]).name not in ALLOWED_COMMANDS:
            raise ValueError("Unsupported command: " + command)
        executable = argv[0]
        expanded = [executable]
        for value in argv[1:]:
            matches = sorted(glob.glob(value), key=_natural_key) if glob.has_magic(value) else []
            expanded.extend(matches or [value])
        argv = expanded
        if Path(executable).name == "rm":
            for value in argv[1:]:
                Path(value).unlink(missing_ok=True)
            data = ""
            continue
        with Path("commands.log").open("a") as log:
            log.write(shlex.join(argv) + "\n")
        if input_name:
            with Path(input_name).open("rb") as input_stream:
                result = subprocess.run(argv, stdin=input_stream, text=True, capture_output=True)
        else:
            result = subprocess.run(argv, input=data, text=True, capture_output=True)
        if output_name:
            Path(output_name).write_text(result.stdout)
        if result.returncode and not (Path(executable).name == "grep" and result.returncode == 1):
            raise RuntimeError(f"Command failed ({result.returncode}): {shlex.join(argv)}\n"
                               + result.stderr[-6000:] + result.stdout[-2000:])
        if Path(executable).name == "vmd" and "psfgen) ERROR" in result.stdout:
            raise RuntimeError("VMD/psfgen reported a structure preparation error:\n"
                               + result.stdout[-6000:])
        data = result.stdout
        with Path("commands.log").open("a") as log:
            log.write(result.stderr)
            if not output_name and Path(executable).name not in {"grep", "cut", "cat", "sed"}:
                log.write(data)
    return data or ""


def run_command(command):
    command_output(command)
    return 0


def command_stream(command):
    return io.StringIO(command_output(command))


def parallel_tasks(function, values, workers):
    """Each mutant/combination owns one set of temporary files."""
    unique = list(dict.fromkeys(values))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(function, unique))


def stage_tools(workspace, config):
    """Create the executable wrappers expected by the legacy engines."""
    workspace = Path(workspace)
    bindir = workspace / "bin"
    bindir.mkdir(parents=True, exist_ok=True)
    runtime_config = workspace / "tools.runtime.json"
    runtime_config.write_text(json.dumps(
        {key: config.get(key) for key in sorted(ALLOWED_KEYS)}, indent=2) + "\n")
    python = sys.executable
    for name, filename in TOOLS.items():
        executable = config.get(name)
        if name == "provean":
            executable = config.get("provean_binary") or executable
        if not executable:
            raise RuntimeError(f"Required tool is not configured: {name}")
        if not Path(executable).is_file() or not os.access(executable, os.X_OK):
            raise RuntimeError(f"Configured tool is missing or not executable: {name}={executable}")
        if name == "dssp":
            argv = [python, str(CORE_DIR / "dssp_adapter.py"), executable]
        elif name == "provean":
            argv = [python, str(CORE_DIR / "provean_adapter.py")]
        else:
            argv = [executable]
        script = "#!/bin/sh\n"
        if name == "provean":
            script += "export MUTABIND2_CONFIG=" + shlex.quote(str(runtime_config)) + "\n"
        script += "exec " + shlex.join(argv) + ' "$@"\n'
        wrapper = bindir / filename
        wrapper.write_text(script)
        wrapper.chmod(0o755)
    (workspace / "foldx").unlink(missing_ok=True)
    (workspace / "foldx").symlink_to(bindir / "foldx")
    # FoldX distributes its rotamer database beside the executable.
    foldx = config.get("foldx")
    if foldx:
        support_dir = Path(foldx).resolve().parent
        for name in ("rotabase.txt", "molecules"):
            support = support_dir / name
            target = workspace / name
            if support.exists() and not target.exists():
                target.symlink_to(support)
    return runtime_config


def load_score_model(mode):
    if mode not in {"single", "multiple"}:
        raise ValueError("mode must be single or multiple")
    name = "mutabinds.RData" if mode == "single" else "mutabindm.RData"
    return RESOURCES / mode / name


def score_file(mode, features, output, rscript=None, contributions=False, config=None):
    config = config or load_config()
    sibling = Path(sys.executable).parent / "Rscript"
    rscript = rscript or config.get("rscript") or (str(sibling) if sibling.is_file() else None) \
        or shutil.which("Rscript")
    if not rscript:
        raise RuntimeError("Rscript is missing; configure it or run `mutabind2 doctor --scope score`")
    model = load_score_model(mode)
    checked([rscript, "--vanilla", str(CORE_DIR / "score.R"), mode, model,
             Path(features).resolve(), Path(output).resolve(), "yes" if contributions else "no"])
    with Path(output).open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def predict_values(mode, features, rscript=None, config=None):
    with tempfile.TemporaryDirectory(prefix="mutabind2_score_") as directory:
        rows = score_file(mode, features, Path(directory) / "prediction.tsv",
                          rscript, False, config)
    return [float(row["SunDDG"]) for row in rows]


def merge_contributions(mode, features, result_file, rscript=None, config=None):
    scored_path = Path(str(features) + ".contribution")
    score_file(mode, features, scored_path, rscript, True, config)
    with scored_path.open(newline="") as stream:
        scores = list(csv.DictReader(stream, delimiter="\t"))
    with Path(result_file).open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields, rows = reader.fieldnames, list(reader)
    columns = [name for name in scores[0] if name in
               {"DDE_vdw", "DDG_solv", "DDG_fold", "CS", "SA_com_wt", "SA_part_wt",
                "N_cont_wt", "dE_vdw_wt"}]
    if len(rows) != len(scores):
        raise ValueError("Contribution and prediction row counts differ")
    with Path(result_file).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields + columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row, score in zip(rows, scores):
            writer.writerow(dict(row, **{name: score[name] for name in columns}))


def configure_engine(mode):
    parser = argparse.ArgumentParser(description=f"MutaBind2 Python 3 {mode} mutation engine")
    parser.add_argument("-i", "--jobid", required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--no-contributions", action="store_true")
    args = parser.parse_args()
    if not re.fullmatch(r"[a-z0-9][a-z0-9_-]*", args.jobid):
        parser.error("jobid must contain lowercase letters, digits, underscores or hyphens")
    workdir = str(args.workdir.resolve()) + "/"
    if not re.fullmatch(r"/[a-z0-9_./-]+/", workdir):
        parser.error("Engine workdir must contain only lowercase ASCII letters, digits, _, . and -")
    if args.workers < 1 or args.threads < 1:
        parser.error("workers and threads must be positive")
    os.chdir(workdir)
    jobpath = workdir + args.jobid
    output = jobpath + "out"
    Path(output).mkdir(exist_ok=True)
    resource = "sinputfiles" if mode == "single" else "minputfiles"
    bindir = workdir + "bin/"
    return dict(workdir=workdir, jobid=args.jobid, jobpath=jobpath, pathoutput=output,
                pathinput=workdir + resource, pathpara=workdir + resource,
                in_file=jobpath + "/" + args.jobid + ".input",
                out_file=jobpath + "/" + args.jobid + ".sunddg",
                pathvmd=bindir, pathcharmm=bindir, pathnamd=bindir,
                pathdssp=bindir, pathprovean=bindir, pathrscript=bindir,
                command_workers=args.workers, provean_threads=args.threads,
                use_contributions=not args.no_contributions)


def check_error_files(output_file):
    for suffix in [".error1", ".error2", ".error3"]:
        path = Path(output_file + suffix)
        if path.exists():
            raise RuntimeError(path.read_text() or f"Engine reported {suffix}")


def stage_foldx_difference(stem, output):
    candidates = [p for p in [Path(f"Dif_{stem}.fxout"), Path(f"Dif_{stem}_BM.fxout")]
                  if p.is_file()]
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one FoldX difference file for {stem}; found {len(candidates)}")
    shutil.move(candidates[0], Path(output) / f"Dif_{stem}.fxout")


def read_foldx_difference(path):
    """Read the documented Total Energy column, independent of banner length."""
    column, values = None, []
    with Path(path).open(newline="") as stream:
        for row in csv.reader(stream, delimiter="\t"):
            fields = [field.strip().lower() for field in row]
            if "pdb" in fields and "total energy" in fields:
                column = fields.index("total energy")
            elif column is not None and row and any(field.strip() for field in row) \
                    and not row[0].startswith("#"):
                if len(row) <= column:
                    raise ValueError(f"Incomplete FoldX energy row: {path}")
                value = float(row[column])
                if not math.isfinite(value):
                    raise ValueError(f"Nonfinite FoldX energy: {path}")
                values.append(value)
    if len(values) != 1:
        raise ValueError(f"Expected one FoldX mutant energy under Total Energy header: {path}")
    return values[0]
