"""Validate original job tables and prepare PDB/PROVEAN inputs without external tools."""

import csv
import re
from collections import OrderedDict
from pathlib import Path
from string import ascii_uppercase

COLUMNS = ["PDBfile", "Partner1", "Partner2", "MutChain", "Mutation_PDB", "Result_Id", "isPI"]
AMINO_ACIDS = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCQEGHILKMFPSTWYV",
))
MUTATION = re.compile(r"([ACDEFGHIKLMNPQRSTVWY])(-?\d+[A-Za-z]?)([ACDEFGHIKLMNPQRSTVWY])\Z")
CHAIN = re.compile(r"[A-Za-z0-9]_\d+\Z")
PDB_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]*\.pdb(?:\d+)?\Z")


def read_job(input_file, mode):
    input_file = Path(input_file)
    with input_file.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != COLUMNS:
            raise ValueError("Input must have exactly these tab-separated columns: " + "\t".join(COLUMNS))
        rows = list(reader)
    if not rows:
        raise ValueError("The mutation table is empty")
    seen_ids = set()
    group_members = {}
    for number, row in enumerate(rows, 2):
        if None in row or any(value is None or value != value.strip() for value in row.values()):
            raise ValueError(f"Row {number}: missing field, extra tab, or leading/trailing whitespace")
        if not PDB_NAME.fullmatch(row["PDBfile"]):
            raise ValueError(f"Row {number}: use a simple PDB filename, without a directory or spaces")
        if row["PDBfile"] != rows[0]["PDBfile"]:
            raise ValueError("Use one protein complex per job")
        if (row["Partner1"], row["Partner2"]) != (rows[0]["Partner1"], rows[0]["Partner2"]):
            raise ValueError("All rows in one job must use the same partners")
        partners = [row[name].split(".") for name in ["Partner1", "Partner2"]]
        chains = partners[0] + partners[1]
        if any(not CHAIN.fullmatch(chain) for chain in chains + [row["MutChain"]]):
            raise ValueError(f"Row {number}: chain labels must look like E_1 or I_1")
        if len(set(chains)) != len(chains) or row["MutChain"] not in chains:
            raise ValueError(f"Row {number}: partners must be disjoint and include the mutated chain")
        # Downstream legacy templates fold chain names to lower case; avoid collisions.
        if len(chains) > len(ascii_uppercase):
            raise ValueError("At most 26 selected chains are supported by this local port")
        match = MUTATION.fullmatch(row["Mutation_PDB"])
        if not match or match[1] == match[3]:
            raise ValueError(f"Row {number}: invalid amino-acid substitution {row['Mutation_PDB']!r}")
        if not re.fullmatch(r"[0-9]+", row["Result_Id"]) or row["Result_Id"] in seen_ids:
            raise ValueError(f"Row {number}: Result_Id must be a unique nonnegative integer")
        seen_ids.add(row["Result_Id"])
        if mode == "single":
            if row["isPI"] != "0":
                raise ValueError("Single mutation jobs must use isPI=0")
        elif not re.fullmatch(r"[A-Za-z][A-Za-z0-9]*", row["isPI"]):
            raise ValueError("Multiple mutation labels must start with a letter and contain letters/digits only")
        if mode == "multiple":
            label = row["isPI"].lower()
            members = group_members.setdefault(label, [])
            member = (row["MutChain"], match[2])
            if member in members:
                raise ValueError(f"Combination {row['isPI']}: the same residue is mutated more than once")
            members.append(member)
            if len(members) > 10:
                raise ValueError("A multiple mutation combination can contain at most 10 substitutions")
    if mode == "multiple":
        labels = {row["isPI"] for row in rows}
        if len(labels) != len({label.lower() for label in labels}):
            raise ValueError("Combination labels cannot differ only in letter case")
    pdb_file = input_file.parent / rows[0]["PDBfile"]
    if not pdb_file.is_file() or not pdb_file.resolve().is_relative_to(input_file.parent.resolve()):
        raise ValueError(f"PDB file is missing or points outside the job directory: {pdb_file.name}")
    return rows


def read_residues(pdb_file):
    """Return complete standard residues indexed by chain/model and PDB numbering."""
    pdb_file = Path(pdb_file)
    all_models = pdb_file.suffix.lower() != ".pdb"
    chains = OrderedDict()
    model = "1"
    first_model = None
    for line in pdb_file.read_text().splitlines():
        if line.startswith("MODEL"):
            model = line.split()[1]
            if first_model is None:
                first_model = model
            if not all_models and model != first_model:
                break
        if not line.startswith("ATOM  ") or len(line) < 54 or line[17:20] not in AMINO_ACIDS:
            continue
        chain = line[21] + "_" + (model if all_models else "1")
        number = line[22:27].strip()
        residues = chains.setdefault(chain, OrderedDict())
        residue = residues.setdefault(number, {"name": line[17:20], "atoms": OrderedDict()})
        # Keep the first alternate coordinate for each atom, as documented upstream.
        residue["atoms"].setdefault(line[12:16].strip(), line)
    for chain in list(chains):
        chains[chain] = OrderedDict(
            (key, residue) for key, residue in chains[chain].items()
            if {"N", "CA", "C", "O"}.issubset(residue["atoms"])
        )
    return chains


def inspect_job(input_file, mode):
    rows = read_job(input_file, mode)
    chains = read_residues(Path(input_file).parent / rows[0]["PDBfile"])
    selected = rows[0]["Partner1"].split(".") + rows[0]["Partner2"].split(".")
    for chain in selected:
        if not chains.get(chain):
            raise ValueError(f"Selected chain {chain} has no complete standard residues")
    for row in rows:
        wt, number, _ = MUTATION.fullmatch(row["Mutation_PDB"]).groups()
        residue = chains[row["MutChain"]].get(number)
        if residue is None:
            raise ValueError(f"Result_Id={row['Result_Id']}: residue {row['MutChain']}:{number} is missing/nonstandard/incomplete")
        actual = AMINO_ACIDS[residue["name"]]
        if actual != wt:
            raise ValueError(f"Result_Id={row['Result_Id']}: expected wild type {wt}, PDB contains {actual}")
    return rows, chains


def prepare_structure(input_file, output_dir, mode):
    rows, chains = inspect_job(input_file, mode)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pdb = rows[0]["PDBfile"].split(".")[0].upper()
    partner1 = rows[0]["Partner1"].split(".")
    selected = partner1 + rows[0]["Partner2"].split(".")
    mapping = dict(zip(selected, ascii_uppercase))
    numbering = {}
    processed = []
    for chain in selected:
        for number, residue in chains[chain].items():
            for line in residue["atoms"].values():
                processed.append(f"{line[:54]:<72}{chain}\n")
    (output_dir / f"{pdb.lower()}_p.pdb").write_text("".join(processed))
    mutated_chains = {row["MutChain"] for row in rows}
    for index, chain in enumerate(selected, 1):
        original, remapped, sequence = [], [], []
        for count, (number, residue) in enumerate(chains[chain].items(), 1):
            if count > 9999:
                raise ValueError("A selected chain exceeds the 9999-residue PDB field limit")
            numbering[(chain, number)] = count
            sequence.append(AMINO_ACIDS[residue["name"]])
            for line in residue["atoms"].values():
                original.append(f"{line[:54]:<72}{chain}\n")
                # Clear alternate-location and insertion fields after renumbering.
                remapped.append(line[:16] + " " + line[17:21] + mapping[chain]
                                + f"{count:4d} " + line[27:72].ljust(45) + chain + "\n")
        (output_dir / f"{pdb}_{chain}.pdb").write_text("".join(original))
        (output_dir / f"{pdb}_CH{index}.pdb").write_text("".join(remapped))
        if chain in mutated_chains:
            (output_dir / f"{pdb}_{chain}.seq").write_text(f"> {pdb}{chain}\n{''.join(sequence)}\n")
            variants = []
            for row in rows:
                if row["MutChain"] == chain:
                    wt, number, mutant = MUTATION.fullmatch(row["Mutation_PDB"]).groups()
                    variants.append(f"{wt}{numbering[(chain, number)]}{mutant}\n")
            (output_dir / f"{pdb}_{chain}.var").write_text("".join(variants))
    new_columns = ["PDBid", "NewPartner1", "NewPartner2", "Mutation_cleaned"]
    with Path(str(input_file) + ".cleaned").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, COLUMNS + new_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            wt, number, mutant = MUTATION.fullmatch(row["Mutation_PDB"]).groups()
            writer.writerow(dict(row, PDBid=pdb,
                NewPartner1="".join(mapping[chain] for chain in partner1),
                NewPartner2="".join(mapping[chain] for chain in selected[len(partner1):]),
                Mutation_cleaned=f"{wt}{mapping[row['MutChain']]}{numbering[(row['MutChain'], number)]}{mutant}"))
    return rows
