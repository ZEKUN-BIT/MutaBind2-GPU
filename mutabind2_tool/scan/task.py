"""Per-batch calculation engine with an optional reused wild-type baseline."""

from datetime import datetime, timezone
import argparse
import csv
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import sys
import tempfile
import time
import signal

from ..core.jobs import stage_job
from ..core.runtime import (PACKAGE_ROOT, TOOL_ROOT, RESOURCES, check_error_files, configured_database,
                            load_config, read_foldx_difference, stage_tools)
from ..core.structure import AMINO_ACIDS, inspect_job, prepare_structure
from .provean import write_engine_outputs

ENGINE = PACKAGE_ROOT / "engines" / "single.py"
PROC = Path("/proc")
PB_WAIT_SECONDS = 15


def now():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, data):
    temporary = Path(str(path) + ".tmp")
    temporary.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")
    temporary.replace(path)


def file_hash(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _segments(expected):
    return len(expected)


def patch_caps(workspace):
    vmd = Path(workspace) / "sinputfiles" / "vmd.pgn"
    text = vmd.read_text()
    old = 'segment CH${j} {pdb $path/${prot}_CH${j}.pdb}'
    new = '''set fp [open $path/${prot}_CH${j}.pdb r]
set firstres ""
while {[gets $fp line] >= 0} {
    if {[string range $line 0 5] eq "ATOM  "} {
        set firstres [string range $line 17 19]
        break
    }
}
close $fp
set ncap ACE
if {$firstres eq "PRO"} {set ncap ACP}
segment CH${j} {
    first $ncap
    last CT2
    pdb $path/${prot}_CH${j}.pdb
}'''
    if text.count(old) != 1:
        raise ValueError("Unexpected VMD template; cannot apply neutral caps")
    vmd.write_text(text.replace(old, new))
    setup = Path(workspace) / "sinputfiles" / "setup.inp"
    text = setup.read_text()
    if text.count("generate setu ch@I") != 1:
        raise ValueError("Unexpected CHARMM setup template; cannot apply neutral caps")
    setup.write_text(text.replace("generate setu ch@I",
                                  "generate setu ch@I first @npatch last CT2"))


def stage(source, workspace, config, cap_mode):
    workspace = Path(workspace)
    workspace.mkdir(parents=True, exist_ok=False)
    stage_job(Path(source), workspace, "single")
    shutil.copytree(RESOURCES / "single", workspace / "sinputfiles")
    stage_tools(workspace, config)
    if cap_mode == "neutral":
        patch_caps(workspace)
        capped = PACKAGE_ROOT / "scan" / "charmm_capped.py"
        wrapper = workspace / "bin" / "charmm"
        # The native adapter recognises this exact transparent wrapper.
        wrapper.write_text('#!/bin/sh\nexec ' + shlex.quote(str(capped)) + ' "$@"\n')
        wrapper.chmod(0o755)
    return workspace


def load_engine(workspace, workers, threads, no_contributions=False):
    spec = importlib.util.spec_from_file_location("mutabind2_single_engine", ENGINE)
    engine = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(engine)
    previous = sys.argv
    sys.argv = ["task_engine", "-i", "job", "--workdir", str(workspace),
                "--workers", str(workers), "--threads", str(threads)]
    if no_contributions:
        sys.argv.append("--no-contributions")
    try:
        engine.__dict__.update(engine.configure_engine("single"))
    finally:
        sys.argv = previous
    return engine


def audit_topology(workspace, stem, expected_segments, write_report=True):
    out = Path(workspace) / "jobout"
    lines = (out / f"{stem}_charmm.psf").read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if "!NATOM" in line)
    count = int(lines[start].split()[0])
    atoms = [line.split() for line in lines[start + 1:start + 1 + count]]
    crd = (out / f"{stem}_charmm.crd").read_text().splitlines()
    coords = []
    for line in crd:
        fields = line.split()
        if len(fields) >= 10 and fields[0].isdigit():
            coords.append(tuple(map(float, fields[4:7])))
    if len(coords) != count or not all(math.isfinite(v) and abs(v) < 10000
                                       for xyz in coords for v in xyz):
        raise ValueError("Uninitialized or inconsistent CHARMM coordinates")
    start = next(i for i, line in enumerate(lines) if "!NBOND" in line)
    bond_count = int(lines[start].split()[0])
    values = []
    for line in lines[start + 1:]:
        if not line.strip():
            break
        values.extend(map(int, line.split()))
    if len(values) != 2 * bond_count:
        raise ValueError("Incomplete PSF bond table")
    bonds = list(zip(values[::2], values[1::2]))
    lengths = [math.dist(coords[a - 1], coords[b - 1]) for a, b in bonds]
    long_bonds = [(a, b, d) for (a, b), d in zip(bonds, lengths) if d > 2.5]
    cross_segment = [(a, b) for a, b in bonds if atoms[a - 1][1] != atoms[b - 1][1]]
    caps = {}
    for atom in atoms:
        caps.setdefault(atom[1], set()).add(atom[4])
    report = {
        "atoms": count, "bonds": bond_count, "maximum_bond_A": max(lengths),
        "bonds_over_2_5_A": long_bonds, "cross_segment_bonds": cross_segment,
        "total_charge": round(sum(float(a[6]) for a in atoms), 6),
        "neutral_caps_present": all({"CAY", "CY", "OY", "NT"} <= names
                                    for names in caps.values()),
        "segments": len(caps), "expected_segments": expected_segments,
    }
    if write_report:
        atomic_json(out / f"topology_{stem}.json", report)
    if long_bonds or cross_segment or report["segments"] != expected_segments:
        raise ValueError(f"Invalid fragment topology: {report}")
    return report


def minimize(engine, stem):
    out = Path(engine.pathoutput)
    template = Path(engine.pathpara)
    text = (template / "rest.pgn").read_text().replace("protname", stem).replace(
        "pathoutput", str(out))
    deck = out / f"rest_{stem}.pgn"
    deck.write_text(text)
    engine.command_output(f"{engine.pathvmd}vmd -dispdev text -e {deck}")
    text = (template / "min.conf").read_text().replace("protname", stem).replace(
        "pathoutput", str(out)).replace("pathinput", str(template)).replace("runstep", "100")
    deck = out / f"min_{stem}.conf"
    deck.write_text(text)
    log = engine.command_output(f"{engine.pathnamd}namd2 {deck}")
    (out / f"min_{stem}.log").write_text(log)
    energies = [line.split() for line in log.splitlines() if line.startswith("ENERGY:")]
    if not energies or energies[-1][1] != "100":
        raise ValueError("NAMD did not finish 100 minimization steps")
    final = energies[-1]
    if not all(math.isfinite(float(value)) for value in final[1:]):
        raise ValueError("Nonfinite minimization energy")
    if float(final[2]) > 10000:
        raise ValueError(f"Excessive final bond energy: {final[2]}")
    engine.run_command(
        f"{engine.pathcharmm}charmm <{template}/frame.inp path={out} protname={stem} "
        f"startframe=100 pathinput={template} >{out}/frame_{stem}.out")
    before = (out / f"{stem}_charmm.crd").read_text().splitlines()
    after = (out / f"{stem}_min.crd").read_text().splitlines()

    def coordinates(lines):
        return [tuple(map(float, line.split()[4:7])) for line in lines
                if len(line.split()) >= 10 and line.split()[0].isdigit()]

    initial, final_coords = coordinates(before), coordinates(after)
    if len(initial) != len(final_coords) or not all(
            math.isfinite(v) and abs(v) < 10000 for xyz in final_coords for v in xyz):
        raise ValueError("Invalid coordinates extracted from NAMD trajectory")
    return {"final_bond_kcal_mol": float(final[2]),
            "final_total_kcal_mol": float(final[11])}


def audit_pb(path):
    text = Path(path).read_text()
    iterations = list(map(int, re.findall(r"Number of iterations:\s+(\d+)", text)))
    maximum = list(map(int, re.findall(
        r"Maximum # iterations\s+\(MAXITS\)\s*=\s*(\d+)", text)))
    values = {name: re.findall(r'Parameter:\s+' + name + r'\s+<-\s+"([^"]+)"', text)
              for name in ["VDWENER", "DELTAGELEC"]}
    if (not re.search(r"^\s*NORMAL TERMINATION\b", text, re.M)
            or "ABNORMAL TERMINATION" in text or len(iterations) != 3
            or len(maximum) != 3):
        raise ValueError(f"Incomplete PB calculation: {path}")
    if "calculation not converged" in text.lower() or any(
            n > limit for n, limit in zip(iterations, maximum)):
        raise ValueError(f"PB solver did not converge: {path}; {iterations}")
    if any(len(v) != 1 or not math.isfinite(float(v[0])) for v in values.values()):
        raise ValueError(f"Missing or nonfinite PB/VDW energies: {path}")
    return {"iterations": iterations, **{k: float(v[0]) for k, v in values.items()}}


def audit_pb_identity(path, workspace, stem):
    report = audit_pb(path)
    text = Path(path).read_text().lower()
    if re.findall(r'parameter:\s+protname\s+<-\s+"([^"]+)"', text) != [stem.lower()]:
        raise ValueError(f"PB output has a different structure identity: {path}")
    for suffix in ["charmm.psf", "min.crd"]:
        if f"::{workspace}/jobout/{stem}_{suffix}::".lower() not in text:
            raise ValueError(f"PB output did not read the expected geometry: {path}")
    return report


def run_pb(workspace, state, operation, function, *parameters):
    lock_path = TOOL_ROOT / "runtime" / "energy.lock"
    lock_path.parent.mkdir(exist_ok=True)
    minimum = 30 * 1024 ** 3

    def update(phase, **details):
        if state.get("stage") != phase:
            state.update(stage=phase, stage_started_utc=now())
        state.update(details, pb_operation=operation, pb_workers=1,
                     pb_global_max_concurrency=1,
                     pb_min_available_gib=round(minimum / 1024 ** 3, 3),
                     pb_lock=str(lock_path))
        atomic_json(Path(workspace) / "status.json", state)

    lock_started = time.monotonic()
    update("waiting_pb_lock", pb_lock_held=False)
    with lock_path.open("a") as lock:
        # Queue in the kernel so a continuously running scan cannot starve a
        # second scan by reacquiring between periodic nonblocking polls.
        fcntl.flock(lock, fcntl.LOCK_EX)
        lock_wait = time.monotonic() - lock_started
        memory_started = time.monotonic()
        try:
            while True:
                available = int(next(
                    line.split()[1] for line in Path("/proc/meminfo").read_text().splitlines()
                    if line.startswith("MemAvailable:"))) * 1024
                update("waiting_pb_memory", pb_lock_held=True,
                       pb_lock_wait_seconds=round(lock_wait, 3),
                       pb_memory_wait_seconds=round(time.monotonic() - memory_started, 3),
                       pb_available_gib=round(available / 1024 ** 3, 3))
                if available >= minimum:
                    break
                time.sleep(PB_WAIT_SECONDS)
            memory_wait = time.monotonic() - memory_started
            started = time.monotonic()
            record = {"lock_wait_seconds": round(lock_wait, 3),
                      "memory_wait_seconds": round(memory_wait, 3), "started_utc": now()}
            state.setdefault("pb_runs", {})[operation] = record
            update(operation, pb_memory_wait_seconds=round(memory_wait, 3))
            try:
                result = function(*parameters)
            except BaseException:
                record["outcome"] = "failed"
                raise
            else:
                record["outcome"] = "complete"
                return result
            finally:
                record.update(run_seconds=round(time.monotonic() - started, 3),
                              finished_utc=now())
                state.setdefault("timings_seconds", {})[operation] = record["run_seconds"]
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)
            state["pb_lock_held"] = False
            atomic_json(Path(workspace) / "status.json", state)


def reuse_wt_contacts(engine, baseline, table):
    resource = "pdb_distance_withH.inp"
    template = (Path(engine.pathpara) / resource).read_text()
    baseline_template = (Path(baseline) / "sinputfiles" / resource).read_text()
    if template != baseline_template or any(name in template for name in ["chmut", "resnum"]):
        raise ValueError("Contact deck depends on the mutated residue; cannot reuse baseline")
    source = Path(baseline) / "jobout" / "distance_job_len_mp_contact_op_10.txt"
    rows = list(csv.DictReader(source.open(), delimiter="\t"))
    counts = {int(row["p1_p2_wt"]) for row in rows}
    if len(counts) != 1:
        raise ValueError("Baseline WT contact counts disagree")
    stem = table[0]["PDBid"].lower()
    current = Path(engine.pathoutput) / f"{stem}_min.crd"
    previous = Path(baseline) / "jobout" / current.name
    if file_hash(current) != file_hash(previous):
        raise ValueError("WT coordinates differ from the contact baseline")
    count = counts.pop()
    destination = Path(engine.pathoutput) / "distance_job_len_mp_contact_op_10.txt"
    with destination.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["PDB_ID", "mutation_cleaned", "p1_p2_wt"])
        writer.writerows((stem, row["Mutation_cleaned"].lower(), count) for row in table)
    atomic_json(Path(engine.pathoutput) / "contact_reuse.json", {
        "baseline": str(baseline), "count": count,
        "source_sha256": file_hash(source),
        "reason": "Identical WT coordinates, partners and residue-independent contact deck."})


def finish_pb(engine, workspace, state, table, step, mutated_partner, baseline=None,
              wt_only=False, backend=None, executable=None, threads=1, vdw_executable=None):
    stem = table[0]["PDBid"].lower()
    label = stem + ("_p1" if mutated_partner == 1 else "_p2")
    mutants = [row["PDBid"].lower() + "_" + row["Mutation_cleaned"].lower()
               for row in table]
    cache = Path(workspace) / "wt_cache"
    cached = Path(baseline) / "wt_cache" / f"energy_{label}.out" if baseline else None
    output = Path(engine.pathoutput) / f"energy_{label}.out"
    if cached is not None and cached.is_file():
        audit = audit_pb(cached)
        shutil.copyfile(cached, output)
        state.setdefault("reused_pb_outputs", {})[output.name] = {
            "validated_utc": now(), **audit}
        state["wt_pb_cache_source"] = str(cached)
    else:
        run_pb(workspace, state, "wt_pb_energy", engine.energy)
        if backend is not None and not backend.valid_cached(workspace, stem, label,
                executable, threads, vdw_executable):
            raise ValueError("Optimized WT energy provenance validation failed")
        audit = step("pb_audit_" + label, audit_pb_identity, output, workspace, stem)
        state["wt_pb_audit"] = audit
    state.setdefault("pb_audits", {})[output.name] = {
        "sha256": file_hash(output), **audit}
    if wt_only or baseline is None:
        cache.mkdir(exist_ok=True)
        shutil.copyfile(output, cache / output.name)
    if wt_only:
        return
    for row, mutant in zip(table, mutants):
        state.update(pb_current_mutation=row["Mutation_PDB"], pb_current_result_id=row["Result_Id"])
        path = Path(engine.pathoutput) / f"energy_{mutant}.out"
        run_pb(workspace, state, "mutant_pb_energy_" + mutant,
               engine.energy_mut_threads, mutant)
        if backend is not None and not backend.valid_cached(
                workspace, mutant, mutant, executable, threads, vdw_executable):
            raise ValueError(f"Optimized PB provenance validation failed: {mutant}")
        audit = step("pb_audit_" + mutant, audit_pb_identity, path, workspace, mutant)
        state.setdefault("pb_audits", {})[path.name] = {
            "sha256": file_hash(path), **audit}


def cache_wt(workspace, table):
    stem = table[0]["PDBid"].lower()
    mutants = [row["PDBid"].lower() + "_" + row["Mutation_cleaned"].lower()
               for row in table]
    cache = Path(workspace) / "wt_cache"
    cache.mkdir(exist_ok=True)
    for path in (Path(workspace) / "jobout").iterdir():
        if not path.is_file() or path.suffix in {".seq", ".var"}:
            continue
        if any(mutant.lower() in path.name.lower() for mutant in mutants):
            continue
        shutil.copyfile(path, cache / path.name)
    return {"files": len(list(cache.iterdir())), "stem": stem}


def validate_baseline(baseline, source, cap_mode, protocol_sha256=None):
    baseline = Path(baseline).resolve()
    status_path = baseline / "status.json"
    if not status_path.is_file():
        raise ValueError(f"Baseline status is missing: {status_path}")
    saved = json.loads(status_path.read_text())
    if saved.get("stage") != "complete" or not saved.get("wt_only"):
        raise ValueError("Baseline is incomplete or was not produced in WT-only mode")
    if saved.get("cap_mode") != cap_mode:
        raise ValueError("Baseline cap mode differs from this batch")
    if protocol_sha256 is not None and saved.get("protocol_sha256") != protocol_sha256:
        raise ValueError("Baseline calculation protocol differs from this batch")
    rows, _ = inspect_job(source, "single")
    pdb_path = Path(source).parent / rows[0]["PDBfile"]
    if saved.get("pdb_sha256") != file_hash(pdb_path):
        raise ValueError("Baseline structure differs from this batch")
    if saved.get("partners") != [rows[0]["Partner1"], rows[0]["Partner2"]]:
        raise ValueError("Baseline partners differ from this batch")
    if not (baseline / "wt_cache").is_dir() or not any((baseline / "wt_cache").iterdir()):
        raise ValueError("Baseline WT cache is missing")
    hashes = saved.get("wt_cache_sha256", {})
    if not hashes or any(not (baseline / "wt_cache" / name).is_file() or
                         file_hash(baseline / "wt_cache" / name) != expected
                         for name, expected in hashes.items()):
        raise ValueError("Baseline cache checksum mismatch")
    return saved


def _expected_residues(source, table):
    rows, chains = inspect_job(source, "single")
    selected = rows[0]["Partner1"].split(".") + rows[0]["Partner2"].split(".")
    expected = {(f"CH{index}", str(number)): AMINO_ACIDS[residue["name"]]
                for index, chain in enumerate(selected, 1)
                for number, residue in enumerate(chains[chain].values(), 1)}
    return expected


def audit_geometry(workspace, stem, expected, cap_mode):
    out = Path(workspace) / "jobout"
    expected_segments = len({key[0] for key in expected})
    topology = audit_topology(workspace, stem, expected_segments, write_report=False)
    lines = (out / f"{stem}_charmm.psf").read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if "!NATOM" in line)
    atoms = [line.split() for line in lines[start + 1:start + 1 + topology["atoms"]]]
    residues, names = {}, {}
    codes = dict(AMINO_ACIDS, HSD="H", HSE="H", HSP="H")
    for index, atom in enumerate(atoms, 1):
        if int(atom[0]) != index or not all(math.isfinite(float(v)) for v in atom[6:8]):
            raise ValueError(f"Invalid PSF atom/charge/mass: {stem}")
        key = (atom[1], atom[2])
        residue = codes[atom[3]]
        if key in residues and residues[key] != residue:
            raise ValueError(f"Inconsistent PSF residue: {stem}, {key}")
        residues[key] = residue
        names.setdefault(key, set()).add(atom[4])
    if residues != expected:
        raise ValueError(f"PSF residue identities do not match requested structure: {stem}")
    for segment in {key[0] for key in names}:
        keys = [key for key in names if key[0] == segment]
        neutral = {"CAY", "CY", "OY"} <= names[keys[0]] and "NT" in names[keys[-1]]
        if neutral != (cap_mode == "neutral"):
            raise ValueError(f"PSF termini disagree with cap mode: {stem}/{segment}")
    for suffix in ["charmm.crd", "min.crd"]:
        path = out / f"{stem}_{suffix}"
        records = [line.split() for line in path.read_text().splitlines()
                   if len(line.split()) >= 10 and line.split()[0].isdigit()]
        if len(records) != len(atoms):
            raise ValueError(f"Coordinate atom count differs from PSF: {path}")
        for atom, record in zip(atoms, records):
            if (record[0], record[7], record[8], record[2], record[3]) != tuple(atom[:5]):
                raise ValueError(f"Coordinate atom identity differs from PSF: {path}")
            if not all(math.isfinite(float(v)) and abs(float(v)) < 10000
                       for v in record[4:7]):
                raise ValueError(f"Invalid coordinates: {path}")
    return topology


def install_backends(engine, workspace, config, state, surface_backend,
                     pb_executable=None, vdw_executable=None, pb_threads=1):
    if surface_backend == "local":
        os.environ["MUTABIND2_CHARMM_NATIVE"] = config.get("charmm") or ""
        if config.get("toolchain_root"):
            os.environ["MUTABIND2_TOOLCHAIN_ROOT"] = config["toolchain_root"]
        from .optimize import surface as surface_module
        surface_module.install(engine, workspace)
        state.update(surface_backend="local", surface_generation=surface_module.GENERATION,
                     surface_implementation_sha256=file_hash(surface_module.__file__))
    else:
        state.update(surface_backend="original",
                     surface_generation="upstream-full-surface")
    backend = None
    if pb_executable:
        from .optimize import energy as energy_module
        backend = energy_module
        backend.install(engine, workspace, Path(pb_executable), pb_threads,
                        Path(vdw_executable) if vdw_executable else None)
        state.update(pb_backend="optimized_split", pb_generation=backend.GENERATION,
                     pb_executable=str(pb_executable),
                     pb_executable_sha256=file_hash(pb_executable), pb_threads=pb_threads,
                     vdw_executable=str(vdw_executable) if vdw_executable else None,
                     vdw_executable_sha256=file_hash(vdw_executable)
                     if vdw_executable else None)
    else:
        state.update(pb_backend="original", pb_generation="upstream-unseeded", pb_threads=1)
    return backend


def _run_foldx(engine, workspace, mutants, workers, state, step):
    from .optimize import foldx
    restore = foldx.install(engine, workspace)
    state.update(foldx_workers=workers, foldx_address_limit_gib=8,
                 foldx_implementation_sha256=file_hash(foldx.__file__))
    try:
        step("foldx", engine.parallel_tasks, engine.runfoldx_mut_threads, mutants, workers)
    finally:
        restore()


def _features_and_predict(engine, workspace, state, table, step, provean_mode,
                          provean_cache, no_contributions, identities):
    if provean_mode in {"full", "cached"}:
        if not provean_cache:
            raise ValueError("--provean-mode full requires --provean-cache")
        step("provean_cached", write_engine_outputs, provean_cache, table, engine.pathoutput, identities)
    else:
        step("provean", engine.RunProvean)
    step("features", engine.getenergy)
    step("predict", engine.Prediction)
    if not no_contributions:
        step("contributions", engine.rf_contribution)
    check_error_files(engine.out_file)
    results = list(csv.DictReader(Path(engine.out_file).open(), delimiter="\t"))
    if len(results) != len(table) or {row["Result_Id"] for row in results} != {
            row["Result_Id"] for row in table}:
        raise ValueError("Prediction IDs do not match the requested batch")
    if not all(math.isfinite(float(row["SunDDG"])) for row in results):
        raise ValueError("Nonfinite binding prediction")
    shutil.copyfile(engine.out_file, Path(workspace) / "results.tsv")
    return results


def _execute(args, source, workspace, config, baseline=None, wt_only=False):
    mutated_partner = int(getattr(args, "mutated_partner", None) or 2)
    state = {
        "stage": "starting", "started_utc": now(), "pid": os.getpid(),
        "mode": "baseline" if wt_only else "batch",
        "source": str(source), "workspace": str(workspace),
        "cap_mode": args.cap_mode, "mutated_partner": mutated_partner,
        "protocol_sha256": args.protocol_sha256,
        "baseline": str(baseline) if baseline else None,
        "protocol": f"fragmented coordinates; {args.cap_mode} termini",
    }
    atomic_json(Path(workspace) / "status.json", state)
    os.environ["BLAST_DB_NR"] = configured_database(config) or ""
    os.environ["RGL_USE_NULL"] = "TRUE"
    os.environ["PYTHONNOUSERSITE"] = "1"
    os.environ["MUTABIND2_CONFIG"] = str(Path(workspace) / "tools.runtime.json")
    os.environ["PATH"] = os.pathsep.join(
        [str(Path(workspace) / "bin"), os.environ.get("PATH", "")])
    engine = load_engine(workspace, args.foldx_workers, args.threads,
                         args.no_contributions)

    def step(name, function, *parameters):
        state.update(stage=name, stage_started_utc=now())
        atomic_json(Path(workspace) / "status.json", state)
        started = time.monotonic()
        value = function(*parameters)
        state.setdefault("timings_seconds", {})[name] = round(
            time.monotonic() - started, 3)
        atomic_json(Path(workspace) / "status.json", state)
        return value

    try:
        step("prepare_structure", engine.prepare_structure, engine.in_file,
             engine.pathoutput, "single")
        table = list(csv.DictReader(Path(str(engine.in_file) + ".cleaned").open(),
                                    delimiter="\t"))
        stem = table[0]["PDBid"].lower()
        mutants = [row["PDBid"] + "_" + row["Mutation_cleaned"]
                   for row in table]
        expected_segments = len(table[0]["NewPartner1"]) + len(table[0]["NewPartner2"])
        state.update(result_ids=[row["Result_Id"] for row in table],
                     pdb_sha256=file_hash(Path(workspace) / "job" / table[0]["PDBfile"]),
                     partners=[table[0]["Partner1"], table[0]["Partner2"]],
                     expected_segments=expected_segments)
        if baseline:
            saved = validate_baseline(baseline, source, args.cap_mode, args.protocol_sha256)
            state.update(baseline_sha256=saved.get("pdb_sha256"))
            for path in (Path(baseline) / "wt_cache").iterdir():
                if path.is_file():
                    shutil.copyfile(path, Path(engine.pathoutput) / path.name)
            state["reused_wt_baseline"] = str(baseline)

        backend = install_backends(
            engine, workspace, config, state, args.surface_backend,
            args.pb_executable, args.vdw_executable, args.pb_threads)

        if not baseline:
            step("wt_pdb", engine.wtpdb)
            step("wt_vmd", engine.vmd_wt)
            if engine.vmderror:
                raise ValueError("Too many atoms for the original workflow")
            step("wt_charmm", engine.charmmfile_wt)
            topology = step("wt_topology", audit_topology, workspace, stem,
                            expected_segments)
            if args.cap_mode == "neutral" and not topology["neutral_caps_present"]:
                raise ValueError("Missing neutral cap atoms")
            step("wt_interface", engine.interface)
            if engine.dsasa >= -100:
                raise ValueError("Partners do not have a sufficient contact interface")
            state["wt_minimization"] = step("wt_minimize", minimize, engine, stem)
            expected = _expected_residues(source, table)
            step("wt_geometry", audit_geometry, workspace, stem, expected, args.cap_mode)
            step("wt_dssp", engine.DSSP)
            if wt_only:
                step("energy_decks", engine.energyfile)
                step("contacts", engine.distance_contact_resid)
                step("extract_contacts", engine.distance_contact_extract)
                step("count_contacts", engine.distance_protein_contact_protein)
                finish_pb(engine, workspace, state, table, step, mutated_partner,
                          baseline=None, wt_only=True, backend=backend,
                          executable=args.pb_executable, threads=args.pb_threads,
                          vdw_executable=args.vdw_executable)
                cached = cache_wt(workspace, table)
                state.update(stage="complete", completed_utc=now(), wt_only=True,
                             predictions=0, wt_cache_files=cached["files"],
                             wt_cache_sha256={p.name:file_hash(p) for p in (Path(workspace)/"wt_cache").iterdir() if p.is_file()})
                atomic_json(Path(workspace) / "status.json", state)
                return state

        step("foldx_inputs", engine.inputfoldx)
        _run_foldx(engine, workspace, mutants, args.foldx_workers, state, step)
        step("split_mutants", engine.splitchain_mut)
        step("mutant_vmd", engine.vmd_mut)
        step("mutant_charmm", engine.charmmfile_mut)
        wt_residues = _expected_residues(source, table)
        for row, mutant in zip(table, mutants):
            topology = step("topology_" + mutant.lower(), audit_topology, workspace,
                            mutant.lower(), expected_segments)
            if args.cap_mode == "neutral" and not topology["neutral_caps_present"]:
                raise ValueError("Missing mutant cap atoms")
            step("minimize_" + mutant.lower(), minimize, engine, mutant.lower())
            expected = dict(wt_residues)
            cleaned = row["Mutation_cleaned"]
            segment = "CH" + str(ord(cleaned[1].upper()) - ord("A") + 1)
            expected[(segment, cleaned[2:-1])] = cleaned[-1].upper()
            step("geometry_" + mutant.lower(), audit_geometry, workspace, mutant.lower(), expected, args.cap_mode)

        step("energy_decks", engine.energyfile)
        step("surface_decks", engine.surfacefile)
        step("surface", engine.surface)
        if baseline:
            try:
                step("reuse_wt_contacts", reuse_wt_contacts, engine, baseline, table)
            except (ValueError, OSError):
                step("contacts", engine.distance_contact_resid)
                step("extract_contacts", engine.distance_contact_extract)
                step("count_contacts", engine.distance_protein_contact_protein)
        else:
            step("contacts", engine.distance_contact_resid)
            step("extract_contacts", engine.distance_contact_extract)
            step("count_contacts", engine.distance_protein_contact_protein)

        finish_pb(engine, workspace, state, table, step, mutated_partner,
                  baseline=baseline, backend=backend, executable=args.pb_executable,
                  threads=args.pb_threads, vdw_executable=args.vdw_executable)
        _features_and_predict(engine, workspace, state, table, step, args.provean_mode,
                              args.provean_cache, args.no_contributions, args.identities)
        state.update(stage="complete", completed_utc=now(), predictions=len(table))
        atomic_json(Path(workspace) / "status.json", state)
        return state
    except BaseException as exc:
        state.update(stage="failed", failed_stage=state.get("stage"), error=str(exc),
                     finished_utc=now())
        atomic_json(Path(workspace) / "status.json", state)
        raise


def run_baseline(source, workspace, config, args):
    stage(source, workspace, config, args.cap_mode)
    return _execute(args, source, workspace, config, baseline=None, wt_only=True)


def run_batch(source, workspace, config, args, baseline=None):
    stage(source, workspace, config, args.cap_mode)
    return _execute(args, source, workspace, config, baseline=baseline, wt_only=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--parameters", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--wt-only", action="store_true")
    parsed = parser.parse_args()
    args = argparse.Namespace(**json.loads(parsed.parameters.read_text()))
    config = load_config(parsed.config)
    def interrupted(signum, frame):
        raise InterruptedError(f"Calculation interrupted by signal {signum}")
    signal.signal(signal.SIGTERM, interrupted)
    if parsed.wt_only:
        run_baseline(parsed.source.resolve(), parsed.workspace.resolve(), config, args)
    else:
        run_batch(parsed.source.resolve(), parsed.workspace.resolve(), config, args,
                  parsed.baseline.resolve() if parsed.baseline else None)


if __name__ == "__main__":
    main()
