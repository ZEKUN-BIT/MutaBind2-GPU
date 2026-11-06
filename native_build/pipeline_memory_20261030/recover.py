"""Continue frozen preparation, including a batch interrupted before SASA."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
INSTALL = HERE.parents[1]
sys.path[:0] = [str(HERE), str(INSTALL / 'native_build/scan_pipeline_20261030'), str(INSTALL)]
import recover_energy as recovery
from memory_admission import install


def validate_prepared(workspace, source, baseline, parameters, state):
    if (state.get('protocol_sha256') != parameters.protocol_sha256
            or state.get('cap_mode') != parameters.cap_mode
            or Path(state.get('baseline', '')) != baseline):
        raise ValueError('Recovery requires the same protocol and baseline')
    if (workspace / 'job/job.input').read_bytes() != source.read_bytes():
        raise ValueError('Staged input differs from the immutable batch')
    rows, _ = recovery.task.inspect_job(source, 'single')
    pdb_name = rows[0]['PDBfile']
    if recovery.digest(workspace / 'job' / pdb_name) != recovery.digest(source.parent / pdb_name):
        raise ValueError('Staged structure differs from the immutable batch')
    recovery.task.validate_baseline(baseline, source, parameters.cap_mode, parameters.protocol_sha256)
    table = recovery.core.read_table(workspace / 'job/job.input.cleaned')
    if [r['Result_Id'] for r in table] != [r['Result_Id'] for r in rows]:
        raise ValueError('Prepared mutation IDs differ from the immutable batch')
    timings = state.get('timings_seconds', {})
    if not {'foldx', 'mutant_vmd', 'mutant_charmm', 'energy_decks', 'surface_decks'} <= set(timings):
        raise ValueError('Preparation checkpoint is incomplete')
    expected_wt = recovery.task._expected_residues(source, table)
    for row in table:
        mutation = row['Mutation_cleaned']
        stem = row['PDBid'].lower() + '_' + mutation.lower()
        expected = dict(expected_wt)
        expected[('CH' + str(ord(mutation[1].upper()) - ord('A') + 1), mutation[2:-1])] = mutation[-1].upper()
        recovery.task.audit_geometry(workspace, stem, expected, parameters.cap_mode)
        if not {'minimize_' + stem, 'geometry_' + stem} <= set(timings):
            raise ValueError('Minimization checkpoint is incomplete')
        lines = (workspace / 'jobout' / f'min_{stem}.log').read_text().splitlines()
        energies = [line.split() for line in lines if line.startswith('ENERGY:')]
        if not energies or energies[-1][1] != '100':
            raise ValueError('Minimization did not finish')
        recovery.task.read_foldx_difference(workspace / 'jobout' / f'Dif_job_{mutation.upper()}.fxout')
    checkpoint = workspace / 'memory_recovery_checkpoint.json'
    evidence = recovery.read(checkpoint)
    if evidence:
        if evidence['protocol_sha256'] != parameters.protocol_sha256:
            raise ValueError('Recovery checkpoint protocol changed')
        for name, expected in evidence['frozen_files'].items():
            if recovery.digest(workspace / name) != expected:
                raise ValueError('Frozen recovery artifact changed: ' + name)
    for name in ('error', 'failed_stage', 'finished_utc'):
        state.pop(name, None)
    return table


def freeze(workspace, source, baseline, parameters, state):
    table = validate_prepared(workspace, source, baseline, parameters, state)
    checkpoint = workspace / 'memory_recovery_checkpoint.json'
    if not checkpoint.exists():
        paths = [p for folder in ('job', 'sinputfiles')
                 for p in (workspace / folder).rglob('*') if p.is_file()]
        for pattern in ('*_min.crd', '*_charmm.crd', '*_charmm.psf', 'min_*.log', 'Dif_job_*.fxout'):
            paths.extend((workspace / 'jobout').glob(pattern))
        recovery.core.atomic_json(checkpoint, dict(created_utc=recovery.task.now(),
            protocol_sha256=parameters.protocol_sha256,
            frozen_files={str(p.relative_to(workspace)): recovery.digest(p) for p in sorted(set(paths))}))
    return table


if __name__ == '__main__':
    install()
    recovery.validate_prepared = validate_prepared
    recovery.main()
