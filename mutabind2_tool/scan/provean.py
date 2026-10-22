"""Full-sequence PROVEAN cache and engine-output adapter.

The MutaBind2 engine normally runs PROVEAN separately for every fragment.  When
a complete sequence is available, conservation can be computed once and reused
for every batch.  Fragment queries are not stitched together; the cache is
validated against the exact full sequence and the exact mutation set.
"""

from datetime import datetime, timezone
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys

from ..core.runtime import PACKAGE_ROOT, configured_database, load_config

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
VARIANT = re.compile(r"([ACDEFGHIKLMNPQRSTVWY])([1-9][0-9]*)([ACDEFGHIKLMNPQRSTVWY])\Z")


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_sequence(path):
    lines = Path(path).read_text().splitlines()
    headers = [index for index, line in enumerate(lines) if line.startswith(">")]
    if headers != [0]:
        raise ValueError("Expected exactly one FASTA record with a leading header")
    sequence = "".join(line.strip() for line in lines[1:])
    if not sequence or set(sequence) - set(AMINO_ACIDS):
        raise ValueError("Sequence must contain only the 20 standard amino acids")
    return sequence


def validate_variant(variant, sequence):
    match = VARIANT.fullmatch(variant)
    if match is None:
        raise ValueError(f"Invalid single substitution: {variant!r}")
    wild, position, mutant = match.groups()
    position = int(position)
    if position > len(sequence) or sequence[position - 1] != wild or wild == mutant:
        raise ValueError(f"Substitution does not match the full sequence: {variant}")
    return wild, position, mutant


def _rows(path):
    with Path(path).open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        lowered = {name.lower(): name for name in (reader.fieldnames or [])}
        if "chain" not in lowered or "mutation" not in lowered:
            raise ValueError("Mutation table needs columns chain<TAB>mutation")
        rows = []
        for number, row in enumerate(reader, 2):
            chain = row[lowered["chain"]].strip()
            mutation = row[lowered["mutation"]].strip().upper()
            if not chain or not mutation:
                raise ValueError(f"Row {number}: empty chain or mutation")
            rows.append((chain, mutation))
    if not rows:
        raise ValueError("Mutation table is empty")
    return rows


def validated_scores(cache):
    cache = Path(cache)
    status = json.loads((cache / "status.json").read_text())
    if status.get("stage") != "complete":
        raise ValueError(f"Full-sequence PROVEAN cache is incomplete: {cache}")
    if status.get("query_file_sha256"):
        if sha256(cache / "query.fasta") != status["query_file_sha256"]:
            raise ValueError("PROVEAN query checksum mismatch")
        sequence = read_sequence(cache / "query.fasta")
        if hashlib.sha256(sequence.encode()).hexdigest() != status["query_sha256"]:
            raise ValueError("PROVEAN query sequence mismatch")
    for name, expected in status.get("source_hashes", {}).items():
        source = cache / name
        if not source.is_file() or sha256(source) != expected:
            raise ValueError(f"PROVEAN provenance changed: {name}")
    scores_path = cache / "scores.tsv"
    if sha256(scores_path) != status["scores_sha256"]:
        raise ValueError("PROVEAN score checksum mismatch")
    scores = {}
    with scores_path.open(newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            key = (row["Chain"], row["Mutation"])
            if key in scores or not math.isfinite(float(row["PROVEAN"])):
                raise ValueError("Duplicate or nonfinite PROVEAN score")
            scores[key] = float(row["PROVEAN"])
    if len(scores) != status.get("score_count", len(scores)) or not scores:
        raise ValueError("PROVEAN score count mismatch")
    return scores


def build_cache(sequence_path, chain, mutations_path, output, threads=8, config=None,
                force=False):
    output = Path(output).resolve()
    sequence = read_sequence(sequence_path)
    if not re.fullmatch(r"[A-Za-z0-9]", chain) or threads < 1:
        raise ValueError("Use a single PDB chain identifier and positive thread count")
    rows = _rows(mutations_path)
    if any(item_chain != chain for item_chain, _ in rows):
        raise ValueError(f"All mutations must belong to chain {chain}")
    for _, mutation in rows:
        validate_variant(mutation, sequence)
    identifiers = [mutation for _, mutation in rows]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("Duplicate substitutions in the PROVEAN request")
    if output.exists() and any(output.iterdir()) and not force:
        status_path = output / "status.json"
        if status_path.is_file() and json.loads(status_path.read_text()).get("stage") == "complete":
            saved = json.loads(status_path.read_text())
            existing = validated_scores(output)
            expected = {(chain, mutation) for _, mutation in rows}
            if (set(existing) == expected and saved.get("query_sha256") ==
                    hashlib.sha256(sequence.encode()).hexdigest()):
                return {"cache": str(output), "scores": len(existing), "reused": True}
        raise ValueError(f"Output directory must be new or empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    config = load_config(config)
    database = configured_database(config)
    if not database:
        raise ValueError("Set blast_db in the configuration or BLAST_DB_NR")
    provean = config.get("provean_binary")
    if not provean:
        raise ValueError("provean_binary is not configured")
    query = output / "query.fasta"
    variants = output / "variants.txt"
    query.write_text(f">{chain}_full_sequence\n{sequence}\n")
    variants.write_text("\n".join(identifiers) + "\n")
    argv = [sys.executable, str(PACKAGE_ROOT / "core" / "provean_adapter.py"),
            "-q", str(query), "-v", str(variants), "-d", database,
            "--num_threads", str(threads), "--quiet"]
    state = {
        "stage": "running", "started_utc": now(), "pid": os.getpid(),
        "query_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "query_file_sha256": sha256(query), "variants_sha256": sha256(variants),
        "sequence_length": len(sequence), "variants": len(identifiers),
        "chain": chain, "database": database, "command": argv,
        "scope": "PROVEAN conservation only; not a binding ddG prediction",
        "sequence_policy": "Complete sequence in original residue numbering",
    }
    (output / "status.json").write_text(json.dumps(state, indent=2) + "\n")
    try:
        environment = dict(os.environ)
        if config.get("config_path"):
            environment["MUTABIND2_CONFIG"] = config["config_path"]
        with (output / "native.log").open("w") as log:
            result = subprocess.run(argv, text=True, stdout=log,
                                    stderr=subprocess.STDOUT, cwd=output,
                                    env=environment)
        if result.returncode:
            raise RuntimeError(f"PROVEAN exited with status {result.returncode}; see native.log")
        native_text = (output / "native.log").read_text()
        if native_text.count("## PROVEAN scores ##") != 1:
            raise ValueError("Missing or repeated native PROVEAN score header")
        scores = {}
        for line in native_text.partition("## PROVEAN scores ##")[2].splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 2 or not math.isfinite(float(fields[1])):
                raise ValueError("Invalid PROVEAN score line")
            if fields[0] in scores:
                raise ValueError("Duplicate PROVEAN score")
            scores[fields[0]] = float(fields[1])
        if set(scores) != set(identifiers):
            raise ValueError(f"Expected {len(identifiers)} scores, obtained {len(scores)}")
        output_scores = output / "scores.tsv"
        with output_scores.open("w", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(["Chain", "Mutation", "PROVEAN"])
            for mutation in identifiers:
                writer.writerow([chain, mutation, scores[mutation]])
        state.update(stage="complete", completed_utc=now(), score_count=len(scores),
                     scores_sha256=sha256(output_scores),
                     native_log_sha256=sha256(output / "native.log"))
        (output / "status.json").write_text(json.dumps(state, indent=2) + "\n")
        validated_scores(output)
        return {"cache": str(output), "scores": len(scores), "reused": False}
    except BaseException as exc:
        state.update(stage="failed", finished_utc=now(), error=str(exc))
        (output / "status.json").write_text(json.dumps(state, indent=2) + "\n")
        raise


def write_engine_outputs(cache, cleaned_table, jobout, identities):
    """Write the per-chain files that the legacy engine reads in getenergy()."""
    scores = validated_scores(cache)
    grouped = {}
    for row in cleaned_table:
        key = tuple(identities[row["Result_Id"]])
        if key not in scores:
            raise ValueError(f"No cached PROVEAN score for {key[0]}:{key[1]}")
        cleaned = row["Mutation_cleaned"].upper()
        local_variant = cleaned[0] + cleaned[2:]
        grouped.setdefault(row["MutChain"], []).append((local_variant, scores[key]))
    for mut_chain, items in grouped.items():
        path = Path(jobout) / f"provean_{cleaned_table[0]['PDBid'].upper()}_{mut_chain}.out"
        lines = [
            "## PROVEAN v1.1 output ##",
            "# Scores reused from a validated full-sequence cache",
            "## PROVEAN scores ##",
            "# VARIATION\tSCORE",
        ]
        for mutation, score in items:
            lines.append(f"{mutation}\t{score}")
        path.write_text("\n".join(lines) + "\n")
    return grouped


def prepare_scan_cache(output, manifest, config, threads):
    """Compute each continuous sequence once for all requested substitutions."""
    from .pdb import MUTATION
    output = Path(output)
    with (output / "plan/mutations.tsv").open() as stream:
        requests = list(csv.DictReader(stream, delimiter="\t"))
    expected = {(r["chain"], r["mutation"]) for r in requests}
    fragments = json.loads((output / "plan/fragments.json").read_text())
    cache = output / "provean"
    identity = {"mode": manifest["provean_mode"], "requests": sorted(expected),
                "structure": manifest["plan"]["prepared_pdb_sha256"]}
    identity_hash = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    status_path = cache / "status.json"
    if status_path.is_file():
        saved = json.loads(status_path.read_text())
        if saved.get("stage") == "complete":
            scores = validated_scores(cache)
            if saved.get("identity_sha256") != identity_hash or set(scores) != expected:
                raise ValueError("PROVEAN cache belongs to a different scan")
            return cache
    cache.mkdir(exist_ok=True)
    scores, provenance = {}, {}
    if manifest["provean_mode"] == "full":
        source = Path(manifest["provean_cache"])
        available = validated_scores(source)
        sequence = read_sequence(source / "query.fasta")
        source_state = json.loads((source / "status.json").read_text())
        chains = {r["chain"] for r in requests}
        if chains != {source_state.get("chain")}:
            raise ValueError("Full-sequence cache must match the single mutated original chain")
        # Require numbering and WT identity to agree over all observed residues.
        for fragment in fragments["fragments"]:
            if fragment["original_chain"] not in chains:
                continue
            mapping = [m for m in fragments["mapping"] if m["fragment_chain"] == fragment["fragment_chain"]]
            for item, wild in zip(mapping, fragment["sequence"], strict=True):
                position = item["original_number"]
                if not position.isdigit() or not 1 <= int(position) <= len(sequence) or sequence[int(position)-1] != wild:
                    raise ValueError("Full-sequence cache numbering/WT differs from the PDB; provide a matching sequence or use fragment mode")
        if not expected <= set(available):
            raise ValueError("Full-sequence cache misses requested substitutions")
        scores = {key: available[key] for key in expected}
        import shutil
        source_dir = cache / "full_source"
        source_dir.mkdir(exist_ok=True)
        for name in ("query.fasta", "scores.tsv", "status.json"):
            shutil.copyfile(source / name, source_dir / name)
            provenance[f"full_source/{name}"] = sha256(source_dir / name)
    else:
        for fragment in fragments["fragments"]:
            chain = fragment["fragment_chain"]
            chosen = [r for r in requests if r["fragment_chain"] == chain]
            if not chosen:
                continue
            mapping = [m for m in fragments["mapping"] if m["fragment_chain"] == chain]
            positions = {m["original_number"]: i for i, m in enumerate(mapping, 1)}
            directory = cache / "queries" / chain
            directory.mkdir(parents=True, exist_ok=True)
            query = directory / "sequence.fasta"
            query.write_text(f">fragment_{chain}\n{fragment['sequence']}\n")
            variants = directory / "mutations.tsv"
            translated = []
            for row in chosen:
                wild, number, mutant = MUTATION.fullmatch(row["mutation"]).groups()
                local = f"{wild}{positions[number]}{mutant}"
                validate_variant(local, fragment["sequence"])
                translated.append((row, local))
            with variants.open("w") as stream:
                writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                writer.writerow(["chain", "mutation"])
                writer.writerows((chain, local) for _, local in translated)
            attempts = sorted(directory.glob("attempt_*"))
            complete = [a for a in attempts if (a / "status.json").exists() and
                        json.loads((a / "status.json").read_text()).get("stage") == "complete"]
            attempt = complete[-1] if complete else directory / f"attempt_{len(attempts)+1:03d}"
            build_cache(query, chain, variants, attempt, threads, config.get("config_path"))
            local_scores = validated_scores(attempt)
            scores.update({(r["chain"], r["mutation"]): local_scores[(chain, local)] for r, local in translated})
            for name in ("query.fasta", "scores.tsv", "status.json", "native.log"):
                provenance[str((attempt / name).relative_to(cache))] = sha256(attempt / name)
    if set(scores) != expected:
        raise ValueError("Incomplete scan conservation cache")
    with (cache / "scores.tsv").open("w") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["Chain", "Mutation", "PROVEAN"])
        writer.writerows((c, m, scores[(c, m)]) for c, m in sorted(scores))
    status_path.write_text(json.dumps({"stage": "complete", "mode": manifest["provean_mode"],
        "completed_utc": now(), "identity_sha256": identity_hash,
        "scores_sha256": sha256(cache / "scores.tsv"), "score_count": len(scores),
        "source_hashes": provenance}, indent=2) + "\n")
    validated_scores(cache)
    return cache
