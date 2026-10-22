"""Resource-bounded, reproducible split CHARMM energy calculation.

The caller MUST hold task_engine.run_pb's existing global lock. This module
does not acquire a second lock. It uses cold SOR and the original physical
parameters; splitting releases VDW neighbor-list memory before PB starts.
"""
from datetime import datetime, timezone
from pathlib import Path
import csv
import hashlib
import json
import math
import os
import re
import resource
import shlex
import signal
import shutil
import subprocess
import threading
import time


GENERATION = "split-cold-sor-clcg4-v1"
# The preceding deployed implementation used the same decks/audits. Only this
# known version may donate a stage after the adapter itself has changed.
LEGACY_IMPLEMENTATION_SHA256 = "eff1fe6e858b794055e6c4ffc3a8c20d8557449417e502271ea09bc37361e382"
SEEDS = (314159265, 22222222, 33333333, 44444444)
GIB = 1024**3
ADDRESS_LIMIT = 20 * GIB
MEMORY_RESERVE = 10 * GIB
TIMEOUT = 1800
_HASH_CACHE = {}
_RUNTIME_CACHE = {}


def _hash(path):
    path = Path(path)
    stat = path.stat()
    identity = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
    cached = _HASH_CACHE.get(str(path))
    if cached and cached[0] == identity:
        return cached[1]
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    value = digest.hexdigest()
    _HASH_CACHE[str(path)] = (identity, value)
    return value


def _text_hash(text):
    return hashlib.sha256(text.encode()).hexdigest()


def _ascii_path(path):
    # abspath preserves the ASCII staging alias; resolve() could expose spaces
    # in the original installation, which CHARMM command parsing cannot handle.
    value = os.path.abspath(os.fspath(path))
    if not re.fullmatch(r"[A-Za-z0-9_./-]+", value):
        raise ValueError(f"CHARMM needs an ASCII path without spaces: {value}")
    return Path(value)


def _name(value):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", value):
        raise ValueError("Invalid structure/result identifier")
    return value


def _atomic_json(path, data):
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(json.dumps(data, indent=2) + "\n")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _runtime_identity(executable):
    """Fingerprint the wrapper, actual ELF and every resolved shared library.

    Accept direct ELF executables and the project's simple exec wrappers. Never
    source/eval a wrapper to discover its paths. Unknown wrappers fail closed.
    """
    executable = _ascii_path(executable)
    with executable.open("rb") as handle:
        prefix = handle.read(4)
    target = executable
    pb_backend = "cpu"
    library_path = os.environ.get("LD_LIBRARY_PATH", "")
    if prefix != b"\x7fELF":
        text = executable.read_text()
        if re.search(r'^export MUTABIND2_PB_BACKEND=cuda$', text, re.M):
            pb_backend = "cuda"
        matches = re.findall(r"^exec ([^\n]+)$", text, re.M)
        libraries = re.findall(r"^export LD_LIBRARY_PATH=([^\n]+)$", text, re.M)
        if len(matches) != 1 or len(libraries) != 1:
            raise ValueError("Unsupported CHARMM wrapper: expected explicit exec and LD_LIBRARY_PATH")
        command = shlex.split(matches[0])
        if len(command) != 2 or command[1] != "$@":
            raise ValueError("Unsupported CHARMM wrapper command")
        root = str(executable.parent.parent)
        target_text = command[0].replace("$mutabind2_root", root)
        library_path = libraries[0].strip('"').replace("$mutabind2_root", root)
        library_path = library_path.replace('${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}',
                                            ':' + os.environ['LD_LIBRARY_PATH'] if os.environ.get('LD_LIBRARY_PATH') else '')
        if "$" in target_text or "$" in library_path:
            raise ValueError("Unsupported variable in CHARMM wrapper")
        target = _ascii_path(target_text)
    with target.open("rb") as handle:
        if handle.read(4) != b"\x7fELF":
            raise ValueError("CHARMM wrapper must exec a native ELF binary")
    key = (str(executable), _hash(executable), str(target), _hash(target), library_path)
    cached = _RUNTIME_CACHE.get(key)
    if cached and all(_hash(path) == digest for path, digest in cached["files"].items()):
        return cached
    env = dict(os.environ, LD_LIBRARY_PATH=library_path)
    result = subprocess.run(["ldd", str(target)], env=env, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=20, check=True)
    if "not found" in result.stdout:
        raise ValueError("Unresolved CHARMM shared library")
    paths = {str(executable), str(target)}
    paths.update(re.findall(r"=>\s+(/\S+)\s+\(", result.stdout))
    paths.update(re.findall(r"^\s*(/\S+)\s+\(", result.stdout, re.M))
    if not any(Path(path).name.startswith("libchmm.so") for path in paths):
        raise ValueError("Could not identify the loaded libchmm.so")
    value = {"executable": str(executable), "target": str(target),
             "pb_backend": pb_backend,
             "ld_library_path": library_path,
             "files": {path: _hash(path) for path in sorted(paths)}}
    _RUNTIME_CACHE[key] = value
    return value


def _make_decks(workspace, label):
    resources = workspace / "sinputfiles"
    template = (workspace / "jobout" / f"energy_{label}.inp").read_text()
    stream = (resources / "vdw_pb.str").read_text()
    match = re.search(r"^stream @pathpara/vdw_pb\.str\s*$", template, re.M | re.I)
    if not match or template[match.end():].strip():
        raise ValueError("Unexpected energy template stream/tail")
    preamble, count = re.subn(r"^set pathpara .*?$", "set pathpara " + str(resources),
                             template[:match.start()], flags=re.M | re.I)
    surface_marker = "! molecular surface (MS) of the complex"
    if count != 1 or stream.count("! Inititalzation") != 1 or stream.count(surface_marker) != 1:
        raise ValueError("Unexpected original energy template markers")
    preamble += "\nrandom clcg iseed " + " ".join(map(str, SEEDS)) + "\n"
    vdw, pb = stream.split("! Inititalzation", 1)
    pb = "! Inititalzation" + pb.split(surface_marker, 1)[0]
    solves = re.findall(r"^SOLVE[^\n]*\n[^\n]+", pb, re.M)
    if (len(solves) != 3 or re.findall(r"^SET (Gab|Ga|Gb) = \?enpb", pb, re.M) != ["Gab", "Ga", "Gb"]
            or any(re.search(r"\b(?:OSOR|FMGR|KEEP\w*|MAXI\w*|DEPS)\b", solve, re.I) for solve in solves)):
        raise ValueError("Expected the original three cold SOR solves")
    for solve in solves:
        first, rest = solve.split("\n", 1)
        if not first.rstrip().endswith("-"):
            raise ValueError("Unexpected SOLVE continuation")
        pb = pb.replace(solve, first.rstrip()[:-1].rstrip() + " MAXI 2000 DEPS 0.000002 -\n" + rest, 1)
    if pb.count("prnlev 5") != 1:
        raise ValueError("Unexpected PB diagnostics initialization")
    pb = pb.replace("prnlev 5", "prnlev 11")
    return {"vdw": preamble + vdw + "\nstop\n", "pb": preamble + pb + "\nstop\n"}


def _context(workspace, stem, label, executable, threads, vdw_executable=None):
    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        raise ValueError("threads must be a positive integer")
    resources = workspace / "sinputfiles"
    decks = _make_decks(workspace, label)
    inputs = [workspace / "jobout" / f"{stem}_{suffix}" for suffix in ("charmm.psf", "min.crd")]
    inputs += [workspace / "jobout" / f"energy_{label}.inp"]
    inputs += sorted(path for path in resources.rglob("*") if path.is_file())
    for required in ("top_all36_prot.rtf", "par_all36_prot.prm", "radii_prot_na_.str"):
        if not (resources / required).is_file():
            raise ValueError(f"Missing PB resource: {required}")
    context = {"generation": GENERATION, "implementation_sha256": _hash(__file__),
               "workspace": str(workspace), "stem": stem, "label": label, "threads": threads,
               "clcg_seeds": list(SEEDS), "runtime": _runtime_identity(executable),
               "inputs": {str(path.relative_to(workspace)): _hash(path) for path in inputs},
               "decks": {name: _text_hash(text) for name, text in decks.items()}}
    if vdw_executable is not None:
        context["vdw_runtime"] = _runtime_identity(vdw_executable)
    return context, decks


def _available():
    return int(next(line.split()[1] for line in Path("/proc/meminfo").read_text().splitlines()
                    if line.startswith("MemAvailable:"))) * 1024


def _limits():
    resource.setrlimit(resource.RLIMIT_AS, (ADDRESS_LIMIT, ADDRESS_LIMIT))


def _terminate_group(proc):
    # Kill the group even if the launcher has already exited: its descendants
    # may still be alive. Reap the direct child and bound the grace period.
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    proc.wait(timeout=5)


def _run_process(executable, workspace, stem, case, mode, threads):
    if _available() < ADDRESS_LIMIT + MEMORY_RESERVE:
        raise RuntimeError("Requires 30 GiB available memory before each native process")
    stage_threads = 1 if mode == "vdw" else threads
    env = dict(os.environ, MUTABIND2_CHARMM_THREADS=str(stage_threads), OMP_NUM_THREADS=str(stage_threads),
               OMP_PROC_BIND="close", OMP_PLACES="cores", OMP_DYNAMIC="FALSE")
    argv = [str(executable), f"path={workspace / 'jobout'}", f"protname={stem}"]
    started = time.monotonic()
    previous_sigterm = None
    proc = None
    def interrupted(signum, frame):
        raise InterruptedError(f"Native calculation interrupted by signal {signum}")
    try:
        if threading.current_thread() is threading.main_thread():
            previous_sigterm = signal.signal(signal.SIGTERM, interrupted)
        with (case / f"{mode}.inp").open("rb") as stdin, (case / f"{mode}.log").open("wb") as stdout:
            proc = subprocess.Popen(argv, cwd=case, stdin=stdin, stdout=stdout, stderr=subprocess.STDOUT,
                                    env=env, preexec_fn=_limits, start_new_session=True)
            while proc.poll() is None:
                if time.monotonic() - started > TIMEOUT:
                    raise RuntimeError(f"{mode}: exceeded {TIMEOUT} second timeout")
                if _available() < MEMORY_RESERVE:
                    raise RuntimeError(f"{mode}: available memory fell below 10 GiB reserve")
                time.sleep(0.5)
            if proc.returncode:
                raise RuntimeError(f"{mode}: CHARMM exited with code {proc.returncode}")
    except BaseException:
        if proc is not None:
            _terminate_group(proc)
        raise
    finally:
        if previous_sigterm is not None:
            signal.signal(signal.SIGTERM, previous_sigterm)
    return {"returncode": proc.returncode, "wall_seconds": time.monotonic() - started, "threads": stage_threads,
            "address_limit_bytes": ADDRESS_LIMIT, "reserve_bytes": MEMORY_RESERVE,
            "timeout_seconds": TIMEOUT}


def _number(value):
    number = float(value.replace("D", "E").replace("d", "e"))
    if not math.isfinite(number):
        raise ValueError("Nonfinite CHARMM output")
    return number


def _values(text, names):
    result = {}
    for name in names:
        found = re.findall(r'Parameter:\s+' + name + r'\s+<-\s+"([^"\n]+)"', text)
        if len(found) != 1:
            raise ValueError(f"Expected exactly one {name} assignment")
        result[name] = _number(found[0])
    return result


def _audit_log(path, mode, workspace, stem, require_gpu=False):
    text = path.read_text()
    if (len(re.findall(r"^\s*NORMAL TERMINATION\b", text, re.M)) != 1
            or re.search(r"ABNORMAL TERMINATION|calculation not converged|Cannot allocate memory|Fortran runtime error|\b(?:NaN|Infinity)\b", text, re.I)):
        raise ValueError(f"Incomplete/failed native {mode} log")
    identities = re.findall(r'Parameter:\s+PROTNAME\s+<-\s+"([^"]+)"', text, re.I)
    if [value.lower() for value in identities] != [stem.lower()]:
        raise ValueError("Native log has the wrong structure identity")
    for suffix in ("charmm.psf", "min.crd"):
        if f"::{workspace}/jobout/{stem}_{suffix}::".lower() not in text.lower():
            raise ValueError("Native log did not read the expected PSF/CRD")
    if mode == "vdw":
        if "Number of iterations:" in text or re.search(r"Parameter:\s+DELTAGELEC\s+<-", text):
            raise ValueError("VDW log contains unexpected PB calculation")
        return _values(text, ("VDWENER",))
    if re.search(r"Parameter:\s+VDWENER\s+<-", text):
        raise ValueError("PB log contains duplicate VDW calculation")
    gpu_successes = re.findall(r'^\s*PBEQ_CUDA_SUCCESS\s+iterations=\s*(\d+)', text, re.M)
    if require_gpu and len(gpu_successes) != 3:
        raise ValueError("Expected three successful CUDA PB solves; CPU fallback is not accepted")
    result = _values(text, ("GAB", "GA", "GB", "DELTAGELEC"))
    iterations = list(map(int, re.findall(r"Number of iterations:\s+(\d+)", text)))
    maximum = list(map(int, re.findall(r"Maximum # iterations\s+\(MAXITS\)\s*=\s*(\d+)", text)))
    tolerance = [_number(x) for x in re.findall(r"Tolerance of convergence\s+\(DEPS\)\s*=\s*(\S+)", text)]
    energies = [_number(x) for x in re.findall(r"Electrostatic energy \[KCAL/MOL\]\s*=\s*(\S+)", text)]
    grids = [list(map(int, x)) for x in re.findall(r"NUMBER OF GRID POINTS:\s+(\d+)\s+(\d+)\s+(\d+)", text)]
    if (len(iterations) != 3 or maximum != [2000] * 3 or tolerance != [0.000002] * 3
            or len(energies) != 3 or len(grids) != 3 or len({tuple(x) for x in grids}) != 1
            or any(size < 3 or size % 2 != 1 for grid in grids for size in grid)
            or any(n < 2 or n > 2000 for n in iterations)):
        raise ValueError("Incomplete or altered three-component cold PB calculation")
    if require_gpu and list(map(int, gpu_successes)) != iterations:
        raise ValueError("CUDA iteration counters disagree with CHARMM convergence logs")
    sections = re.split(r"Maximum # iterations", text)[1:]
    final_residuals = []
    for section, last in zip(sections, iterations):
        trace = re.findall(r"N,ANORM,DEPS,OMEGA,PHI\s+(\d+)([^\n]*)", section)
        if [int(row[0]) for row in trace] != list(range(1, last + 1)):
            raise ValueError("Missing PB residual trace")
        values = []
        for _, row in trace:
            # Fortran 4E11.4 has no separator before a negative PHI value.
            numbers = re.findall(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)[EeDd][-+]?\d+", row)
            if len(numbers) != 4:
                raise ValueError("Malformed PB residual trace")
            values.append([_number(x) for x in numbers])
        residual, deps = values[-1][:2]
        # A true value just below DEPS may round to equality in E11.4 output.
        if residual < 0 or residual > 0.000002 or deps != 0.000002:
            raise ValueError("PB residual does not meet the original threshold")
        final_residuals.append(residual)
    return {**result, "iterations": iterations, "energies": energies, "grids": grids,
            "final_residuals": final_residuals}


def valid_cached(workspace, stem, label, executable, threads=1, vdw_executable=None):
    """Return False for legacy, incomplete, changed or unverified energy output."""
    try:
        workspace = _ascii_path(workspace)
        stem, label = _name(stem), _name(label)
        case = workspace / "jobout" / "optimized_energy" / label
        manifest = json.loads((case / "manifest.json").read_text())
        context, _ = _context(workspace, stem, label, executable, threads, vdw_executable)
        if (manifest.get("status") != "complete" or manifest.get("generation") != GENERATION
                or manifest.get("context") != context):
            return False
        expected = {"vdw.inp", "pb.inp", "vdw.log", "pb.log"}
        if set(manifest["artifacts"]) != expected:
            return False
        if any(_hash(case / name) != digest for name, digest in manifest["artifacts"].items()):
            return False
        for mode in ("vdw", "pb"):
            process = manifest["processes"][mode]
            if process["returncode"] != 0 or process["threads"] != (1 if mode == "vdw" else threads):
                return False
            if _hash(case / f"{mode}.inp") != context["decks"][mode]:
                return False
            if _audit_log(case / f"{mode}.log", mode, workspace, stem,
                          context["runtime"].get("pb_backend") == "cuda") != manifest["audits"][mode]:
                return False
        merged = workspace / "jobout" / f"energy_{label}.out"
        if _hash(merged) != manifest["merged_sha256"]:
            return False
        values = _values(merged.read_text(), ("VDWENER", "DELTAGELEC"))
        return values == {"VDWENER": manifest["audits"]["vdw"]["VDWENER"],
                          "DELTAGELEC": manifest["audits"]["pb"]["DELTAGELEC"]}
    except (OSError, ValueError, KeyError, TypeError, RuntimeError, subprocess.SubprocessError):
        return False


def _reusable_stages(case, context, workspace, stem):
    """Reaudit each completed native stage against its own exact provenance.

    A VDW executable change must not invalidate an unchanged PB calculation.
    Unknown adapter versions, changed input files and incomplete manifests do
    not donate results. The old directory is archived by the caller afterward.
    """
    try:
        manifest = json.loads((case / "manifest.json").read_text())
        previous = manifest["context"]
        if (manifest.get("status") != "complete" or manifest.get("generation") != GENERATION
                or previous.get("implementation_sha256") not in
                {context["implementation_sha256"], LEGACY_IMPLEMENTATION_SHA256}):
            return {}
        for key in ("generation", "workspace", "stem", "label", "clcg_seeds", "inputs"):
            if previous.get(key) != context[key]:
                return {}
    except (OSError, ValueError, KeyError, TypeError):
        return {}
    reused = {}
    for mode in ("vdw", "pb"):
        try:
            old_runtime = previous.get("vdw_runtime", previous["runtime"]) if mode == "vdw" else previous["runtime"]
            new_runtime = context.get("vdw_runtime", context["runtime"]) if mode == "vdw" else context["runtime"]
            expected_threads = 1 if mode == "vdw" else context["threads"]
            process = manifest["processes"][mode]
            if (old_runtime != new_runtime or process["returncode"] != 0
                    or process["threads"] != expected_threads
                    or (mode == "pb" and previous["threads"] != context["threads"])
                    or previous["decks"][mode] != context["decks"][mode]):
                continue
            for suffix in ("inp", "log"):
                name = f"{mode}.{suffix}"
                if _hash(case / name) != manifest["artifacts"][name]:
                    raise ValueError("Changed native artifact")
            if _hash(case / f"{mode}.inp") != context["decks"][mode]:
                continue
            audit = _audit_log(case / f"{mode}.log", mode, workspace, stem,
                               new_runtime.get("pb_backend") == "cuda")
            if audit != manifest["audits"][mode]:
                continue
            reused[mode] = {"process": process, "audit": audit,
                            "log_sha256": manifest["artifacts"][f"{mode}.log"],
                            "source_manifest_sha256": _hash(case / "manifest.json")}
        except (OSError, ValueError, KeyError, TypeError):
            continue
    return reused


def calculate_energy(workspace, stem, label, executable, threads=1, vdw_executable=None):
    """Calculate and publish one energy; caller holds the outer PB lock."""
    workspace, executable = _ascii_path(workspace), _ascii_path(executable)
    stem, label = _name(stem), _name(label)
    if vdw_executable is not None:
        vdw_executable = _ascii_path(vdw_executable)
    if valid_cached(workspace, stem, label, executable, threads, vdw_executable):
        return json.loads((workspace / "jobout/optimized_energy" / label / "manifest.json").read_text())
    context, decks = _context(workspace, stem, label, executable, threads, vdw_executable)
    case = workspace / "jobout/optimized_energy" / label
    reusable = _reusable_stages(case, context, workspace, stem)
    merged = workspace / "jobout" / f"energy_{label}.out"
    if case.exists():
        previous_case = case.with_name(case.name + f".previous.{time.time_ns()}")
        case.rename(previous_case)
    if merged.exists():
        merged.rename(merged.with_name(merged.name + f".interrupted.{time.time_ns()}"))
    case.mkdir(parents=True)
    for mode, text in decks.items():
        (case / f"{mode}.inp").write_text(text)
    manifest = {"generation": GENERATION, "status": "running", "context": context,
                "started_utc": datetime.now(timezone.utc).isoformat(), "processes": {}, "audits": {}}
    _atomic_json(case / "manifest.json", manifest)
    temporary = merged.with_name(merged.name + f".tmp.{os.getpid()}")
    try:
        for mode in ("vdw", "pb"):
            selected = vdw_executable if mode == "vdw" and vdw_executable is not None else executable
            if mode in reusable:
                if _hash(previous_case / "manifest.json") != reusable[mode]["source_manifest_sha256"]:
                    raise ValueError("Source stage manifest changed after validation")
                shutil.copyfile(previous_case / f"{mode}.log", case / f"{mode}.log")
                if _hash(case / f"{mode}.log") != reusable[mode]["log_sha256"]:
                    raise ValueError("Source stage log changed after validation")
                manifest["processes"][mode] = reusable[mode]["process"]
                manifest.setdefault("reused_stages", {})[mode] = {
                    "source": str(previous_case),
                    "source_manifest_sha256": reusable[mode]["source_manifest_sha256"]}
            else:
                manifest["processes"][mode] = _run_process(selected, workspace, stem, case, mode, threads)
            manifest["audits"][mode] = _audit_log(case / f"{mode}.log", mode, workspace, stem,
                                                 context["runtime"].get("pb_backend") == "cuda")
            if mode in reusable and manifest["audits"][mode] != reusable[mode]["audit"]:
                raise ValueError("Reused stage audit changed after validation")
            _atomic_json(case / "manifest.json", manifest)
        if _context(workspace, stem, label, executable, threads, vdw_executable)[0] != context:
            raise ValueError("Input/runtime contents changed during native calculation")
        if any(_hash(case / f"{mode}.inp") != context["decks"][mode] for mode in ("vdw", "pb")):
            raise ValueError("Native input deck changed during calculation")
        vdw_text = (case / "vdw.log").read_text()
        assignment = re.findall(r'^.*Parameter:\s+VDWENER\s+<-\s+"[^"\n]+"[^\n]*$', vdw_text, re.M)
        if len(assignment) != 1:
            raise ValueError("Cannot merge ambiguous VDW assignment")
        combined = assignment[0] + "\n" + (case / "pb.log").read_text()
        _values(combined, ("VDWENER", "DELTAGELEC"))
        temporary.write_text(combined)
        manifest["artifacts"] = {name: _hash(case / name) for name in ("vdw.inp", "pb.inp", "vdw.log", "pb.log")}
        manifest["merged_sha256"] = _hash(temporary)
        temporary.replace(merged)
        manifest.update(status="complete", finished_utc=datetime.now(timezone.utc).isoformat())
        _atomic_json(case / "manifest.json", manifest)
        return manifest
    except BaseException as exc:
        manifest.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        _atomic_json(case / "manifest.json", manifest)
        raise
    finally:
        temporary.unlink(missing_ok=True)


def install(engine, workspace, executable, threads=1, vdw_executable=None):
    """Install adapters; return audit_cached(stem, label) -> bool for resume."""
    workspace, executable = _ascii_path(workspace), _ascii_path(executable)
    def energy():
        seen, results = set(), []
        with Path(str(engine.in_file) + ".cleaned").open() as handle:
            rows = csv.reader(handle, delimiter="\t")
            next(rows)
            for row in rows:
                stem, mutation = row[7].lower(), row[10].strip().lower()
                partner = "p1" if mutation[1].upper() in row[8] else "p2"
                label = stem + "_" + partner
                if label not in seen:
                    results.append(calculate_energy(workspace, stem, label, executable, threads, vdw_executable))
                    seen.add(label)
        return results
    def energy_mut_threads(stem):
        return calculate_energy(workspace, stem, stem, executable, threads, vdw_executable)
    engine.energy = energy
    engine.energy_mut_threads = energy_mut_threads
    return lambda stem, label: valid_cached(workspace, stem, label, executable, threads, vdw_executable)
