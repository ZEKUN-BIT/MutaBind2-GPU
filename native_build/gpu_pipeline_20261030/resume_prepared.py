"""Resume frozen mutant geometries donated by the preceding PB protocol."""
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
INSTALL=HERE.parents[1]
sys.path[:0]=[str(INSTALL/'native_build/scan_pipeline_20261030'),str(INSTALL)]
import recover_energy as recovery


def validate(workspace,source,baseline,parameters,state):
    evidence=recovery.read(workspace/'gpu_migration.json')
    if not evidence or evidence['new_protocol']!=parameters.protocol_sha256:
        raise ValueError('Missing protocol migration evidence')
    if state.get('protocol_sha256')!=parameters.protocol_sha256 or Path(state.get('baseline',''))!=baseline:
        raise ValueError('Migrated worker protocol/baseline mismatch')
    if state.get('cap_mode')!=parameters.cap_mode or state.get('pb_runs'):
        raise ValueError('Only a fully prepared batch without energy work may migrate')
    if (workspace/'job/job.input').read_bytes()!=source.read_bytes():
        raise ValueError('Migrated input differs from the new immutable batch')
    for name,expected in evidence['scientific_files'].items():
        if recovery.digest(workspace/name)!=expected:
            raise ValueError('Donated artifact changed: '+name)
    recovery.task.validate_baseline(baseline,source,parameters.cap_mode,parameters.protocol_sha256)
    table=recovery.core.read_table(workspace/'job/job.input.cleaned')
    rows,_=recovery.task.inspect_job(source,'single')
    if [r['Result_Id'] for r in rows]!=[r['Result_Id'] for r in table]:
        raise ValueError('Migrated mutation IDs differ from the plan')
    expected_wt=recovery.task._expected_residues(source,table)
    timings=state.get('timings_seconds',{})
    if not {'foldx','mutant_vmd','mutant_charmm','surface','energy_decks'}<=set(timings):
        raise ValueError('Incomplete preparation checkpoint')
    for row in table:
        mutation=row['Mutation_cleaned']
        stem=row['PDBid'].lower()+'_'+mutation.lower()
        expected=dict(expected_wt)
        expected[('CH'+str(ord(mutation[1].upper())-ord('A')+1),mutation[2:-1])]=mutation[-1].upper()
        recovery.task.audit_geometry(workspace,stem,expected,parameters.cap_mode)
        if not {'minimize_'+stem,'geometry_'+stem}<=set(timings):
            raise ValueError('Missing frozen minimization')
        lines=(workspace/'jobout'/f'min_{stem}.log').read_text().splitlines()
        energies=[line.split() for line in lines if line.startswith('ENERGY:')]
        if not energies or energies[-1][1]!='100':
            raise ValueError('Donated minimization did not finish')
        recovery.task.read_foldx_difference(workspace/'jobout'/f'Dif_job_{mutation.upper()}.fxout')
    return table


if __name__=='__main__':
    recovery.validate_prepared=validate
    recovery.main()
