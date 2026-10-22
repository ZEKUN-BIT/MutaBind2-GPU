"""PDB parsing, automatic fragment splitting and mutation-site expansion.

A MutaBind2 job cannot represent a chain with missing loops: every selected
chain must be continuous.  This module detects physical peptide breaks and
splits each original chain into separate fragments, preserving original residue
numbers.  All downstream numbering and mutation strings stay in the original
PDB frame.
"""

from collections import OrderedDict
from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
import re
from string import ascii_uppercase
import numpy as np

STANDARD = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
MUTATION = re.compile(r"([ACDEFGHIKLMNPQRSTVWY])(-?\d+[A-Za-z]?)([ACDEFGHIKLMNPQRSTVWY])\Z")
CHAIN_ID = re.compile(r"[A-Za-z0-9]\Z")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _element(line, atom):
    value = line[76:78].strip() if len(line) >= 78 else ""
    if value:
        return value.upper()
    return atom.lstrip("0123456789")[0].upper()


@dataclass
class Residue:
    chain: str
    number: str
    name3: str
    atoms: "OrderedDict[str, str]" = field(default_factory=OrderedDict)
    heavy: list = field(default_factory=list)

    @property
    def name1(self):
        return STANDARD[self.name3]

    @property
    def number_int(self):
        match = re.match(r"(-?\d+)", self.number)
        return int(match.group(1)) if match else 0

    @property
    def insertion(self):
        match = re.match(r"-?\d+(.*)", self.number)
        return match.group(1) if match else ""

    @property
    def complete(self):
        return {"N", "CA", "C", "O"}.issubset(self.atoms)

    def atom_xyz(self, name):
        line = self.atoms.get(name)
        if line is None:
            return None
        return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


@dataclass
class Fragment:
    chain: str
    original_chain: str
    residues: list
    partner: int
    index: int

    @property
    def start(self):
        return self.residues[0].number

    @property
    def end(self):
        return self.residues[-1].number

    @property
    def sequence(self):
        return "".join(residue.name1 for residue in self.residues)


@dataclass
class Structure:
    path: Path
    chains: "OrderedDict[str, list]"
    seqres: dict
    atom_records: int
    hetatm_records: int
    models: int
    headers: list

    def residue(self, chain, number):
        for residue in self.chains.get(chain, []):
            if residue.number == number:
                return residue
        return None

    def residues(self, chain):
        return list(self.chains.get(chain, []))


def parse_pdb(path):
    path = Path(path)
    chains = OrderedDict()
    seqres = OrderedDict()
    headers = []
    atom_records = 0
    hetatm_records = 0
    models = 0
    for line in path.read_text().splitlines():
        if line.startswith("MODEL"):
            models += 1
            if models > 1:
                raise ValueError("Multi-model PDB files are not supported by the scan workflow; "
                                 "extract one model first")
            continue
        if line.startswith("SEQRES"):
            chain = line[11]
            seqres.setdefault(chain, []).extend(line[19:].split())
            continue
        if line.startswith("HETATM"):
            hetatm_records += 1
            continue
        if not line.startswith("ATOM  ") or len(line) < 54:
            continue
        atom_records += 1
        residue_name = line[17:20]
        if residue_name not in STANDARD:
            raise ValueError(f"Nonstandard ATOM residue {residue_name!r}; prepare a standard protein PDB first")
        if line[16] not in {" ", "A"}:
            continue
        line = line[:16] + " " + line[17:]
        chain = line[21]
        if not CHAIN_ID.fullmatch(chain):
            raise ValueError(f"PDB chain identifier is missing or unsupported at: {line[:30]}")
        number = line[22:27].strip()
        xyz = tuple(float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54)))
        if not all(math.isfinite(v) for v in xyz):
            raise ValueError(f"Nonfinite coordinates at {chain}:{number}")
        residues = chains.setdefault(chain, [])
        if not residues or residues[-1].number != number:
            if any(r.number == number for r in residues):
                raise ValueError(f"Repeated noncontiguous residue {chain}:{number}")
            residue = Residue(chain=chain, number=number, name3=residue_name)
            residues.append(residue)
            if len(headers) < 3:
                headers.append(line[:80])
        else:
            residue = residues[-1]
            if residue.name3 != residue_name:
                raise ValueError(f"Inconsistent residue name at {chain}:{number}")
        atom = line[12:16].strip()
        if atom in residue.atoms:
            continue
        residue.atoms.setdefault(atom, line)
        if _element(line, atom) not in {"H", "D"}:
            residue.heavy.append((float(line[30:38]), float(line[38:46]),
                                  float(line[46:54])))
    if not chains:
        raise ValueError(f"No standard protein ATOM records found: {path}")
    for chain in chains:
        for residue in chains[chain]:
            if not residue.atoms:
                raise ValueError(f"Empty residue at {chain}:{residue.number}")
    return Structure(path=path, chains=chains, seqres=seqres, atom_records=atom_records,
                     hetatm_records=hetatm_records, models=max(models, 1), headers=headers)


def _distance(a, b):
    return math.dist(a, b)


def _break(previous, current, threshold):
    if current.number_int < previous.number_int:
        return True
    # Author numbering may skip an integer in a continuous peptide. Physical
    # connectivity, rather than that numbering gap alone, defines a break.
    c_terminal = previous.atom_xyz("C")
    n_terminal = current.atom_xyz("N")
    if c_terminal is None or n_terminal is None:
        return True
    return _distance(c_terminal, n_terminal) > threshold


def _segments(residues, threshold):
    result = []
    current = []
    for residue in residues:
        if current and _break(current[-1], residue, threshold):
            result.append(current)
            current = []
        current.append(residue)
    if current:
        result.append(current)
    return result


def _explicit_fragments(structure, spec, partner_of):
    result = []
    for entry in spec:
        chain = entry["chain"]
        start = str(entry["start"])
        end = str(entry["end"])
        residues = structure.chains.get(chain, [])
        first = int(re.match(r"-?\d+", start).group())
        last = int(re.match(r"-?\d+", end).group())
        selected = [r for r in residues if first <= r.number_int <= last]
        if not selected:
            raise ValueError(f"Explicit fragment {chain}:{start}-{end} selects no residues")
        if selected[0].number != start or selected[-1].number != end:
            raise ValueError(f"Explicit fragment {chain}:{start}-{end} boundaries are not observed residues")
        result.append((chain, selected))
    return result


def build_fragments(structure, partner1, partner2, mode="auto", gap_threshold=2.0,
                    explicit=None):
    if not partner1 or not partner2:
        raise ValueError("Both partner chain lists must be non-empty")
    selected = list(partner1) + list(partner2)
    if len(set(selected)) != len(selected):
        raise ValueError("Partner chains must be disjoint")
    missing = [chain for chain in selected if chain not in structure.chains]
    if missing:
        raise ValueError("Selected chains are absent from the PDB: " + ", ".join(missing))
    partner_of = {chain: 1 for chain in partner1}
    partner_of.update({chain: 2 for chain in partner2})
    if mode == "explicit":
        groups = _explicit_fragments(structure, explicit or [], partner_of)
        if any(chain not in partner_of for chain, _ in groups):
            raise ValueError("Explicit fragments must belong to the selected partners")
        observed = {(c, r.number) for c in selected for r in structure.chains[c]}
        covered = {(c, r.number) for c, residues in groups for r in residues}
        if observed != covered:
            raise ValueError("Explicit fragments must cover every observed residue of both partners")
        if any(not r.complete for _, residues in groups for r in residues):
            raise ValueError("Explicit fragments contain incomplete backbone atoms")
    else:
        groups = []
        for chain in selected:
            residues = structure.chains[chain]
            bad = [residue.number for residue in residues if not residue.complete]
            if bad:
                raise ValueError(f"Chain {chain} has incomplete residues (missing N/CA/C/O): "
                                 + ", ".join(bad[:5]))
            if mode == "none":
                if len(_segments(residues, gap_threshold)) != 1:
                    raise ValueError(f"Chain {chain} has breaks; use --fragment-mode auto or explicit")
                groups.append((chain, residues))
            elif mode == "auto":
                for segment in _segments(residues, gap_threshold):
                    groups.append((chain, segment))
            else:
                raise ValueError("fragment mode must be auto, none or explicit")
    if len(groups) > len(ascii_uppercase):
        raise ValueError("At most 26 fragments are supported by the legacy engine")
    fragments = []
    mapping = {}
    for index, (original_chain, residues) in enumerate(groups):
        chain = ascii_uppercase[index]
        fragment = Fragment(chain=chain, original_chain=original_chain,
                            residues=list(residues), partner=partner_of[original_chain],
                            index=index + 1)
        fragments.append(fragment)
        for residue in residues:
            key = (original_chain, residue.number)
            if key in mapping:
                raise ValueError(f"Residue selected twice: {original_chain}:{residue.number}")
            mapping[key] = chain
    partner1_chains = [f.chain for f in fragments if f.partner == 1]
    partner2_chains = [f.chain for f in fragments if f.partner == 2]
    if not partner1_chains or not partner2_chains:
        raise ValueError("Fragmenting removed one entire binding partner")
    return fragments, mapping, ".".join(c + "_1" for c in partner1_chains), \
        ".".join(c + "_1" for c in partner2_chains)


def fragment_pdb_text(structure, fragments):
    lines = [
        "REMARK 900 MUTABIND2 GENERAL SCAN: DISCONTINUOUS CHAINS SPLIT",
        "REMARK 900 ORIGINAL RESIDUE NUMBERS AND INSERTION CODES RETAINED",
        "REMARK 900 TERMINAL CAPS ARE ADDED DOWNSTREAM, NOT IN THIS FILE",
    ]
    for fragment in fragments:
        for residue in fragment.residues:
            for line in residue.atoms.values():
                lines.append(line[:21] + fragment.chain + line[22:])
        lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n"


def nearest_distance(points, targets):
    targets = np.asarray(targets, dtype=float)
    best = float("inf")
    for point in points:
        delta = targets - point
        best = min(best, float(np.einsum("ij,ij->i", delta, delta).min()))
    return math.sqrt(best)


def annotate_sites(structure, fragments, mutated_partner, only=None):
    """Describe each residue in the mutated partner without claiming confidence."""
    partner_fragments = [f for f in fragments if f.partner != mutated_partner]
    partner_atoms = np.asarray([xyz for fragment in partner_fragments
                     for residue in fragment.residues for xyz in residue.heavy], dtype=float)
    terminal_residues = []
    for fragment in fragments:
        full = structure.seqres.get(fragment.original_chain, [])
        observed = structure.chains[fragment.original_chain]
        direct_mapping = full and all(not r.insertion and 1 <= r.number_int <= len(full)
                                      and full[r.number_int - 1] == r.name3 for r in observed)
        if direct_mapping:
            start_truncated = fragment.residues[0].number_int > 1
            end_truncated = fragment.residues[-1].number_int < len(full)
        else:
            start_truncated = fragment.residues[0] is not observed[0]
            end_truncated = fragment.residues[-1] is not observed[-1]
        if start_truncated:
            terminal_residues.append(fragment.residues[0])
        if end_truncated:
            terminal_residues.append(fragment.residues[-1])
    terminal_atoms = np.asarray([xyz for residue in terminal_residues for xyz in residue.heavy], dtype=float)
    annotations = {}
    for fragment in fragments:
        if fragment.partner != mutated_partner:
            continue
        for residue in fragment.residues:
            if only is not None and (residue.chain, residue.number) not in only:
                continue
            partner_distance = nearest_distance(residue.heavy, partner_atoms) \
                if len(partner_atoms) else float("inf")
            terminal_distance = nearest_distance(residue.heavy, terminal_atoms) \
                if len(terminal_atoms) else float("inf")
            annotations[(residue.chain, residue.number)] = {
                "original_chain": residue.chain,
                "original_number": residue.number,
                "wild_type": residue.name1,
                "fragment_chain": fragment.chain,
                "fragment_number": residue.number,
                "partner_min_heavy_distance_A": round(partner_distance, 4),
                "truncated_terminus_min_heavy_distance_A": round(terminal_distance, 4),
                "is_truncated_terminus": int(any(residue is r for r in terminal_residues)),
                "near_truncated_terminus_10A": int(terminal_distance <= 10.0),
            }
    return annotations


def read_mutation_table(path):
    """Read ``chain<TAB>mutation`` rows; optional id/note columns are accepted."""
    import csv
    path = Path(path)
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Mutation table has no header")
        lowered = {name.lower(): name for name in reader.fieldnames}
        for required in ("chain", "mutation"):
            if required not in lowered:
                raise ValueError("Mutation table needs columns: chain<TAB>mutation")
        rows = []
        for number, row in enumerate(reader, 2):
            chain = row[lowered["chain"]].strip()
            mutation = row[lowered["mutation"]].strip()
            if not chain or not mutation:
                raise ValueError(f"Row {number}: empty chain or mutation")
            result_id = row[lowered.get("id", "")].strip() if "id" in lowered else ""
            note = row[lowered.get("note", "")].strip() if "note" in lowered else ""
            rows.append({"chain": chain, "mutation": mutation,
                         "result_id": result_id, "note": note})
    if not rows:
        raise ValueError("Mutation table is empty")
    return rows


def parse_sites(value):
    """Parse ``L:90-166,200`` or ``A:10,12`` into (chain, first, last) tuples."""
    result = []
    for part in value.split(","):
        part = part.strip()
        if not part:
            continue
        match = re.fullmatch(r"([A-Za-z0-9]):(-?\d+)(?:-(-?\d+))?", part)
        if not match:
            raise ValueError(f"Invalid site selector: {part!r}; use CHAIN:START-END")
        chain = match.group(1)
        first = int(match.group(2))
        last = int(match.group(3) or first)
        if first > last:
            raise ValueError(f"Site range start exceeds end: {part}")
        result.append((chain, first, last))
    if not result:
        raise ValueError("No site selectors were provided")
    return result


def expand_sites(structure, selectors):
    rows = []
    seen = set()
    for chain, first, last in selectors:
        if chain not in structure.chains:
            raise ValueError(f"Site chain is absent from the PDB: {chain}")
        selected = [r for r in structure.chains[chain]
                    if first <= r.number_int <= last]
        if not selected:
            raise ValueError(f"No observed residues in {chain}:{first}-{last}")
        for residue in selected:
            if not residue.complete:
                raise ValueError(f"Site {chain}:{residue.number} is incomplete")
            for mutant in AMINO_ACIDS:
                if mutant == residue.name1:
                    continue
                mutation = f"{residue.name1}{residue.number}{mutant}"
                if (chain, residue.number, mutant) in seen:
                    continue
                seen.add((chain, residue.number, mutant))
                rows.append({"chain": chain, "mutation": mutation,
                             "result_id": "", "note": ""})
    if not rows:
        raise ValueError("Site expansion produced no substitutions")
    return rows


def validate_mutations(structure, rows):
    result = []
    seen = set()
    for number, row in enumerate(rows, 2):
        chain = row["chain"]
        mutation = row["mutation"].upper()
        match = MUTATION.fullmatch(mutation)
        if not match:
            raise ValueError(f"Row {number}: invalid substitution {row['mutation']!r}")
        wild, residue_number, mutant = match.groups()
        if wild == mutant:
            raise ValueError(f"Row {number}: wild type and mutant are identical")
        residue = structure.residue(chain, residue_number)
        if residue is None:
            raise ValueError(f"Row {number}: {chain}:{residue_number} is absent from the PDB")
        if not residue.complete:
            raise ValueError(f"Row {number}: {chain}:{residue_number} is incomplete")
        if residue.name1 != wild:
            raise ValueError(f"Row {number}: expected {wild} at {chain}:{residue_number}, "
                             f"PDB contains {residue.name1}")
        if (chain, mutation) in seen:
            raise ValueError(f"Duplicate substitution: {chain}:{mutation}")
        seen.add((chain, mutation))
        result.append(dict(row, chain=chain, mutation=mutation))
    return result


def write_plan(output, structure, fragments, mapping, partner1, partner2, annotations):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    pdb_path = output / "fragments.pdb"
    pdb_path.write_text(fragment_pdb_text(structure, fragments))
    plan = {
        "source": str(structure.path),
        "source_sha256": sha256(structure.path),
        "prepared_pdb": str(pdb_path),
        "prepared_pdb_sha256": sha256(pdb_path),
        "partner1": partner1,
        "partner2": partner2,
        "fragments": [
            {
                "fragment_chain": f.chain,
                "original_chain": f.original_chain,
                "start": f.start,
                "end": f.end,
                "partner": f.partner,
                "residue_count": len(f.residues),
                "sequence": f.sequence,
            }
            for f in fragments
        ],
        "mapping": [
            {"original_chain": key[0], "original_number": key[1], "fragment_chain": value}
            for key, value in mapping.items()
        ],
        "validation": {
            "atom_records": structure.atom_records,
            "hetatm_records": structure.hetatm_records,
            "fragments": len(fragments),
            "mutated_partner_residues": len(annotations),
        },
    }
    (output / "fragments.json").write_text(
        json.dumps(plan, ensure_ascii=False, indent=2) + "\n")
    import csv
    if annotations:
        fields = list(next(iter(annotations.values())))
        with (output / "site_annotations.tsv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(annotations.values())
    return plan
