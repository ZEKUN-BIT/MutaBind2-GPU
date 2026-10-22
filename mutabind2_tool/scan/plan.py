"""Build a resumable scan plan from a PDB, partner chains and mutation input."""

from datetime import datetime, timezone
import csv
import hashlib
import json
import os
import math
from pathlib import Path
import re

from .. import __version__
from ..core.runtime import load_config
from ..core.structure import COLUMNS, inspect_job
from . import pdb as pdb_module


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_chains(value):
    chains = [part.strip() for part in value.split(",") if part.strip()]
    for chain in chains:
        if not re.fullmatch(r"[A-Za-z0-9]", chain):
            raise ValueError(f"Invalid chain identifier: {chain!r}")
    if len(set(chains)) != len(chains):
        raise ValueError("Partner chain lists contain duplicates")
    if not chains:
        raise ValueError("Partner chain list is empty")
    return chains


def assign_result_ids(rows):
    explicit = [row for row in rows if row.get("result_id")]
    used = set()
    for row in explicit:
        value = str(row["result_id"]).strip()
        if not re.fullmatch(r"[0-9]+", value):
            raise ValueError("Explicit result IDs must be nonnegative integers")
        row["result_id"] = value
        if value in used:
            raise ValueError(f"Duplicate result ID: {value}")
        used.add(value)
    next_id = max((int(value) for value in used), default=0) + 1
    for row in rows:
        if not row.get("result_id"):
            while str(next_id) in used:
                next_id += 1
            row["result_id"] = str(next_id)
            used.add(str(next_id))
            next_id += 1
    return rows


def _write_table(path, fieldnames, rows):
    temporary = Path(str(path) + ".tmp")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _existing_scan(output):
    manifest = Path(output) / "scan.json"
    if manifest.is_file():
        return json.loads(manifest.read_text())
    return None


def ensure_output(output, resume=False):
    output = Path(output).resolve()
    if output.exists():
        if not output.is_dir():
            raise ValueError(f"Output path is not a directory: {output}")
        contents = [p for p in output.iterdir() if p.name != "scan.lock"]
        if contents and not (resume and (output / "scan.json").is_file()):
            raise ValueError(f"Output directory must be new or empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    return output


def build(args):
    if args.resume and (Path(args.output) / "scan.json").is_file():
        return load(args.output), {}
    if not args.pdb:
        raise ValueError("--pdb is required for a new scan")
    if bool(args.mutations) == bool(args.sites):
        raise ValueError("Provide exactly one of --mutations or --sites")
    if not math.isfinite(args.gap_threshold) or args.gap_threshold <= 0:
        raise ValueError("--gap-threshold must be finite and positive")
    if args.max_per_batch < 1 or args.workers not in (1, 2):
        raise ValueError("--max-per-batch must be positive and --workers must be 1 or 2")
    config = load_config(args.config)
    pb = args.pb_executable or config.get("pb_executable")
    if getattr(args, "backend", "gpu") == "cpu":
        pb = args.pb_executable or config.get("charmm")
    vdw = args.vdw_executable or config.get("vdw_executable")
    if not pb:
        raise ValueError("Configure pb_executable for GPU, or explicitly select --backend cpu")
    if args.provean_mode == "full" and not args.provean_cache:
        raise ValueError("Full-sequence mode requires --provean-cache; create it with scan provean")
    output = ensure_output(args.output, args.resume)
    previous = _existing_scan(output) if args.resume else None
    source = Path(args.pdb).resolve()
    if not source.is_file():
        raise ValueError(f"PDB file not found: {source}")
    if args.partner1 is None or args.partner2 is None:
        raise ValueError("Both --partner1 and --partner2 are required")
    partner1 = parse_chains(args.partner1)
    partner2 = parse_chains(args.partner2)
    structure = pdb_module.parse_pdb(source)
    explicit = None
    if args.fragment_mode == "explicit":
        if not args.fragments:
            raise ValueError("--fragment-mode explicit requires --fragments JSON")
        explicit = json.loads(Path(args.fragments).read_text())
        if not isinstance(explicit, list):
            raise ValueError("Fragment specification must be a JSON list")
    fragments, mapping, partner1_string, partner2_string = pdb_module.build_fragments(
        structure, partner1, partner2, mode=args.fragment_mode,
        gap_threshold=args.gap_threshold, explicit=explicit)

    if args.mutations:
        rows = pdb_module.read_mutation_table(args.mutations)
    elif args.sites:
        rows = pdb_module.expand_sites(structure, pdb_module.parse_sites(args.sites))
    else:
        raise ValueError("Provide --mutations or --sites")
    rows = pdb_module.validate_mutations(structure, rows)
    rows = assign_result_ids(rows)

    mutation_chains = {row["chain"] for row in rows}
    unknown = mutation_chains - set(partner1) - set(partner2)
    if unknown:
        raise ValueError("Mutation chains are not in either partner: " + ", ".join(sorted(unknown)))
    if args.mutated_partner:
        mutated_partner = int(args.mutated_partner)
        expected = set(partner1 if mutated_partner == 1 else partner2)
    elif mutation_chains <= set(partner1):
        mutated_partner = 1
        expected = set(partner1)
    elif mutation_chains <= set(partner2):
        mutated_partner = 2
        expected = set(partner2)
    else:
        raise ValueError("Mutations span both binding partners; run separate scans. "
                         "The WT baseline and energy orientation differ per partner.")
    if not mutation_chains <= expected:
        raise ValueError("Mutation chains do not belong to the selected mutated partner")

    plan_directory = output / "plan"
    plan_directory.mkdir(exist_ok=True)
    only = {(row["chain"], pdb_module.MUTATION.fullmatch(row["mutation"]).group(2))
            for row in rows}
    annotations = {} if args.no_annotations else pdb_module.annotate_sites(
        structure, fragments, mutated_partner, only=only)
    plan = pdb_module.write_plan(plan_directory, structure, fragments, mapping,
                                 partner1_string, partner2_string, annotations)
    # write_plan stores the source PDB text; use the remapped file for jobs.
    plan["prepared_pdb"] = str(plan_directory / "fragments.pdb")
    plan["prepared_pdb_sha256"] = sha256(plan_directory / "fragments.pdb")

    input_rows = []
    detailed = []
    for row in rows:
        chain = row["chain"]
        fragment_chain = mapping[(chain, pdb_module.MUTATION.fullmatch(row["mutation"]).group(2))]
        input_rows.append({
            "PDBfile": "fragments.pdb",
            "Partner1": partner1_string,
            "Partner2": partner2_string,
            "MutChain": fragment_chain + "_1",
            "Mutation_PDB": row["mutation"],
            "Result_Id": row["result_id"],
            "isPI": "0",
        })
        detailed.append({
            "result_id": row["result_id"],
            "chain": chain,
            "mutation": row["mutation"],
            "fragment_chain": fragment_chain,
            "note": row.get("note", ""),
        })
    table_path = plan_directory / "mutations.input"
    _write_table(table_path, COLUMNS, input_rows)
    _write_table(plan_directory / "mutations.tsv",
                 ["result_id", "chain", "mutation", "fragment_chain", "note"], detailed)
    # The engine's own validator is the final authority on the staged table.
    inspect_job(table_path, "single")

    config = load_config(args.config)
    manifest = {
        "schema": 1,
        "version": __version__,
        "created_utc": previous.get("created_utc") if previous else now(),
        "updated_utc": now(),
        "stage": "planned",
        "source": str(source),
        "source_sha256": sha256(source),
        "partner1": partner1,
        "partner2": partner2,
        "mutated_partner": mutated_partner,
        "fragment_mode": args.fragment_mode,
        "gap_threshold": args.gap_threshold,
        "cap_mode": args.cap_mode,
        "surface_backend": args.surface_backend,
        "provean_mode": args.provean_mode,
        "provean_cache": str(Path(args.provean_cache).resolve()) if args.provean_cache else None,
        "max_per_batch": args.max_per_batch,
        "workers": args.workers,
        "foldx_workers": args.foldx_workers,
        "backend": getattr(args, "backend", "gpu"),
        "pb_executable": str(Path(pb).resolve()),
        "vdw_executable": str(Path(vdw).resolve()) if vdw else None,
        "config": config.get("config_path"),
        "plan": {
            "directory": str(plan_directory),
            "prepared_pdb": plan["prepared_pdb"],
            "prepared_pdb_sha256": plan["prepared_pdb_sha256"],
            "input": str(table_path),
            "input_sha256": sha256(table_path),
            "mutations": len(input_rows),
            "sites": len({item["chain"] + ":" + pdb_module.MUTATION
                          .fullmatch(item["mutation"]).group(2) for item in detailed}),
            "fragments": len(fragments),
        },
        "mutation_identity": [
            [item["result_id"], item["chain"], item["mutation"]] for item in detailed
        ],
    }
    if previous:
        for key in ("source_sha256", "partner1", "partner2", "mutated_partner",
                    "mutation_identity"):
            if previous.get(key) != manifest.get(key):
                raise ValueError(f"Existing scan plan differs in {key}; use a new output directory")
    (output / "scan.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")

    batches = plan_batches(output, manifest, annotations)
    manifest["batches"] = batches["count"]
    manifest["batch_order"] = batches["order"]
    manifest["plan_hashes"] = {str(path.relative_to(output)): sha256(path)
                              for path in [table_path, plan_directory / "fragments.pdb",
                                           plan_directory / "mutations.tsv",
                                           plan_directory / "fragments.json", output / "batches.json"]}
    if (plan_directory / "site_annotations.tsv").exists():
        manifest["plan_hashes"]["plan/site_annotations.tsv"] = sha256(plan_directory / "site_annotations.tsv")
    (output / "scan.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    return manifest, annotations


def _read_rows(output):
    with (Path(output) / "plan" / "mutations.input").open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def _detailed_rows(output):
    path = Path(output) / "plan" / "mutations.tsv"
    with path.open(newline="") as stream:
        return {row["result_id"]: row for row in csv.DictReader(stream, delimiter="\t")}


def plan_batches(output, manifest, annotations=None):
    """Group substitutions by site, interface-proximal first."""
    output = Path(output)
    rows = _read_rows(output)
    detailed = _detailed_rows(output)
    if not rows:
        raise ValueError("No mutations to batch")
    if annotations is None:
        annotations = {}
        annotation_path = output / "plan" / "site_annotations.tsv"
        if annotation_path.is_file():
            with annotation_path.open(newline="") as stream:
                for row in csv.DictReader(stream, delimiter="\t"):
                    annotations[(row["original_chain"], row["original_number"])] = row
    sites = {}
    for row in rows:
        match = pdb_module.MUTATION.fullmatch(row["Mutation_PDB"])
        key = (row["MutChain"], match.group(2))
        sites.setdefault(key, []).append(row)
    order = []
    for key, site_rows in sites.items():
        mutation = pdb_module.MUTATION.fullmatch(site_rows[0]["Mutation_PDB"])
        original = mutation.group(2)
        original_chain = detailed.get(site_rows[0]["Result_Id"], {}).get("chain")
        annotation = annotations.get((original_chain, original))
        distance = float(annotation["partner_min_heavy_distance_A"]) \
            if annotation and annotation.get("partner_min_heavy_distance_A") else float("inf")
        order.append((distance, key, site_rows))
    order.sort(key=lambda item: (item[0], item[1]))
    max_per_batch = max(1, int(manifest.get("max_per_batch") or 19))
    batches = []
    for _, key, site_rows in order:
        for start in range(0, len(site_rows), max_per_batch):
            chunk = site_rows[start:start + max_per_batch]
            batches.append({
                "index": len(batches) + 1,
                "mut_chain": key[0],
                "residue": key[1],
                "result_ids": [row["Result_Id"] for row in chunk],
                "mutations": [row["Mutation_PDB"] for row in chunk],
            })
    batch_directory = output / "batches"
    batch_directory.mkdir(parents=True, exist_ok=True)
    for batch in batches:
        directory = batch_directory / f"batch_{batch['index']:04d}"
        directory.mkdir(exist_ok=True)
        fields = COLUMNS
        rows_for_batch = [row for row in rows
                          if row["Result_Id"] in set(batch["result_ids"])]
        _write_table(directory / "batch.input", fields, rows_for_batch)
        pdb_target = directory / "fragments.pdb"
        if not pdb_target.exists():
            os.link(output / "plan" / "fragments.pdb", pdb_target)
    (output / "batches.json").write_text(json.dumps(
        {"count": len(batches), "order": [batch["index"] for batch in batches],
         "batches": batches}, ensure_ascii=False, indent=2) + "\n")
    return {"count": len(batches), "order": [batch["index"] for batch in batches]}


def load(output):
    output = Path(output).resolve()
    manifest_path = output / "scan.json"
    if not manifest_path.is_file():
        raise ValueError(f"Not a scan directory (scan.json missing): {output}")
    manifest = json.loads(manifest_path.read_text())
    for name, expected in manifest.get("plan_hashes", {}).items():
        path = output / name
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"Scan plan file changed or missing: {path}")
    for key, path in (("input_sha256", output / "plan" / "mutations.input"),
                      ("prepared_pdb_sha256", output / "plan" / "fragments.pdb")):
        expected = manifest["plan"][key]
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"Scan plan file changed or missing: {path}")
    return manifest


def summary(manifest):
    plan = manifest["plan"]
    return {
        "source": manifest["source"],
        "mutations": plan["mutations"],
        "sites": plan["sites"],
        "fragments": plan["fragments"],
        "batches": manifest.get("batches"),
        "cap_mode": manifest["cap_mode"],
        "mutated_partner": manifest["mutated_partner"],
        "partner1": ".".join(manifest["partner1"]),
        "partner2": ".".join(manifest["partner2"]),
    }
