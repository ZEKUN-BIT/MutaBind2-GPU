"""Validate one isolated REEN candidate against retained CPU PHI hashes.

Uses the existing sequential benchmark's memory guard and global lock. Large
PHI files are retained until the parent task records and performs cleanup.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re

from benchmark_pb import ROOT, run_case


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=('wt', 'e90a'), required=True)
    parser.add_argument('--label', required=True)
    parser.add_argument('--executable', type=Path, required=True)
    args = parser.parse_args()
    if not re.fullmatch(r'[a-z0-9_]+', args.label):
        parser.error('label must contain only lowercase letters, digits and underscores')
    reference_name, stem = {
        'wt': ('large_wt_cudab_save', '7mkj2610_frag'),
        'e90a': ('large_mut_e90a_cudab_save', '7mkj2610_frag_ee1a'),
    }[args.case]
    reference = ROOT / reference_name
    cache = ROOT / (args.label + '_phi')
    if (ROOT / args.label).exists():
        parser.error('benchmark label already exists')
    cache.mkdir(exist_ok=False)
    deck = (reference / 'input.inp').read_text()
    deck, replaced = re.subn(
        r'(?m)^(open write file unit 40 name )\S+/([^/\s]+\.phi)$',
        lambda match: match[1] + str(cache / match[2]), deck)
    if replaced != 3:
        raise ValueError('Expected three PHI output paths')
    result = run_case(args.label, deck, ROOT.parent / 'baseline_01/jobout',
                      stem, args.executable.resolve(),
                      reference=reference / 'output.log', require_gpu=True)
    checks = []
    expected = json.loads((reference / 'phi_validation.json').read_text())
    for item in expected['phi_files']:
        path = cache / item['file']
        with path.open('rb') as stream:
            checksum = hashlib.file_digest(stream, 'sha256').hexdigest()
        checks.append(dict(file=path.name, bytes=path.stat().st_size,
                           cpu_sha256=item['cpu_sha256'],
                           candidate_sha256=checksum,
                           bitwise_equal=checksum == item['cpu_sha256']))
    log = (ROOT / args.label / 'output.log').read_text()
    reen_solves = len(re.findall(r'^\s*PBEQ_REEN_CUDA_SUCCESS\b', log, re.M))
    comparison = result['comparison']
    identical = (all(value == 0 for value in comparison['values_difference'].values())
                 and all(value == 0 for value in comparison['energies_difference'])
                 and all(comparison[key] for key in
                         ('grids_equal', 'boxes_equal', 'iterations_equal')))
    validation = dict(all_phi_files_equal=all(item['bitwise_equal'] for item in checks),
                      phi_files=checks, reen_gpu_solves=reen_solves,
                      logged_results_identical=identical)
    (ROOT / args.label / 'phi_validation.json').write_text(
        json.dumps(validation, indent=2) + '\n')
    print(json.dumps(validation), flush=True)
    if not validation['all_phi_files_equal'] or not identical or reen_solves != 3:
        raise RuntimeError('Candidate did not satisfy the full native REEN validation')


if __name__ == '__main__':
    main()
