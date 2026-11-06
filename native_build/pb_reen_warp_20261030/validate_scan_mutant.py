"""Compare both backends on an archived mutant from the current scan."""
from pathlib import Path
import hashlib
import json
import re
import shutil
import sys
import tarfile

HERE = Path(__file__).resolve().parent
INSTALL = Path('/home/zek/mutabind2-tool')
sys.path.insert(0, str(INSTALL))
from mutabind2_tool.scan.optimize import energy
from benchmark_pb import ROOT, run_case


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    scan = INSTALL / 'projects/7mkj_128_255_20261030'
    workspace = HERE / 'validation/scan_mutant'
    workspace.mkdir(exist_ok=False)
    (workspace / 'jobout').mkdir()
    shutil.copytree(scan / 'baseline/sinputfiles', workspace / 'sinputfiles')
    archive = scan / 'batches/batch_0001/scientific_intermediates.tar.gz'
    with tarfile.open(archive) as tar:
        name = next(n for n in tar.getnames()
                    if n.startswith('jobout/optimized_energy/') and n.endswith('/manifest.json'))
        manifest = json.loads(tar.extractfile(name).read())
        context = manifest['context']
        stem, label = context['stem'], context['label']
        for name in [f'jobout/{stem}_charmm.psf', f'jobout/{stem}_min.crd',
                     f'jobout/energy_{label}.inp']:
            path = workspace / name
            path.write_bytes(tar.extractfile(name).read())
            assert digest(path) == context['inputs'][name], name
        reference = workspace / 'reference_pb.log'
        reference.write_bytes(tar.extractfile(f'jobout/optimized_energy/{label}/pb.log').read())
    deck = energy._make_decks(workspace, label)['pb']
    reports, phi_hashes = {}, {}
    for backend, executable in [('baseline', INSTALL / 'bin/charmm-pb-gpu'),
                                ('candidate', HERE / 'bin/charmm')]:
        case_label = f'warp_20261030_scan_mutant_{backend}'
        phi = workspace / (backend + '_phi')
        phi.mkdir()
        def save_phi(match):
            component = match[1].lower()
            return (match[0] + f'\nopen write file unit 40 name {phi / (component + ".phi")}\n'
                    'write phi unit 40\nclose unit 40')
        saved_deck, count = re.subn(r'^SET (Gab|Ga|Gb) = \?enpb[^\n]*$', save_phi, deck, flags=re.M)
        assert count == 3
        result = run_case(case_label, saved_deck, workspace / 'jobout', stem,
                          executable, reference=reference, require_gpu=True)
        cmp = result['comparison']
        assert all(v == 0 for v in cmp['values_difference'].values())
        assert all(v == 0 for v in cmp['energies_difference'])
        assert all(cmp[k] for k in ['grids_equal', 'boxes_equal', 'iterations_equal'])
        log = (ROOT / case_label / 'output.log').read_text()
        for marker in ['PBEQ_REEN_CUDA_SUCCESS', 'PBEQ_BOUNDARY_CUDA_SUCCESS', 'PBEQ_CUDA_SUCCESS']:
            assert len(re.findall(r'^\s*' + marker + r'\b', log, re.M)) == 3
        phi_hashes[backend] = {p.name: digest(p) for p in sorted(phi.glob('*.phi'))}
        assert len(phi_hashes[backend]) == 3
        reports[backend] = result
    assert phi_hashes['baseline'] == phi_hashes['candidate']
    report = {'stem': stem, 'label': label, 'archived_inputs_verified': True,
              'logged_results_identical': True, 'three_phi_files_bitwise_equal': True,
              'phi_sha256': phi_hashes, 'backends': reports}
    (HERE / 'validation/scan_mutant_comparison.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': True, 'stem': stem,
                      'wall_seconds': {k: v['wall_seconds'] for k, v in reports.items()}}), flush=True)


if __name__ == '__main__':
    main()
