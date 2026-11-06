"""Recover the interrupted, fully prepared batch using original engine functions."""
import argparse
import json
import os
from pathlib import Path
import signal
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1])]
from pipeline import core, digest, plan, read, task
from pipeline_worker import wait_for_energy
from mutabind2_tool.core.runtime import configured_database, load_config


def validate_prepared(workspace, source, baseline, parameters, state):
    if (state.get('protocol_sha256') != parameters.protocol_sha256
            or state.get('cap_mode') != parameters.cap_mode
            or Path(state.get('baseline', '')) != baseline
            or not state.get('pb_runs')):
        raise ValueError('Recovery requires a prepared batch with the same protocol and baseline')
    if (workspace / 'job/job.input').read_bytes() != source.read_bytes():
        raise ValueError('Staged input differs from the immutable batch')
    rows, _ = task.inspect_job(source, 'single')
    pdb_name = rows[0]['PDBfile']
    if digest(workspace / 'job' / pdb_name) != digest(source.parent / pdb_name):
        raise ValueError('Staged structure differs from the immutable batch')
    task.validate_baseline(baseline, source, parameters.cap_mode, parameters.protocol_sha256)
    table = core.read_table(workspace / 'job/job.input.cleaned')
    if [r['Result_Id'] for r in table] != [r['Result_Id'] for r in rows]:
        raise ValueError('Prepared mutation IDs differ from the immutable batch')
    required = {'foldx', 'mutant_vmd', 'mutant_charmm', 'surface', 'energy_decks'}
    if not required <= set(state.get('timings_seconds', {})):
        raise ValueError('Preparation checkpoint is incomplete')
    expected_wt = task._expected_residues(source, table)
    for row in table:
        mutation = row['Mutation_cleaned']
        stem = row['PDBid'].lower() + '_' + mutation.lower()
        expected = dict(expected_wt)
        expected[('CH' + str(ord(mutation[1].upper()) - ord('A') + 1), mutation[2:-1])] = mutation[-1].upper()
        task.audit_geometry(workspace, stem, expected, parameters.cap_mode)
        timings = state['timings_seconds']
        if not {'minimize_' + stem, 'geometry_' + stem} <= set(timings):
            raise ValueError('Minimization checkpoint is incomplete')
        lines = (workspace / 'jobout' / f'min_{stem}.log').read_text().splitlines()
        energies = [line.split() for line in lines if line.startswith('ENERGY:')]
        if not energies or energies[-1][1] != '100':
            raise ValueError('Minimization did not finish')
        task.read_foldx_difference(workspace / 'jobout' / f'Dif_job_{mutation.upper()}.fxout')
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('module', choices=['mutabind2_tool.scan.task'])
    for name in ['source', 'workspace', 'parameters', 'config', 'baseline']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    workspace, source, baseline = args.workspace.resolve(), args.source.resolve(), args.baseline.resolve()
    parameters = argparse.Namespace(**read(args.parameters))
    config = load_config(args.config)
    output = args.config.resolve().parent
    _, protocol = core.protocol(plan.load(output), config, output / 'provean')
    if protocol != parameters.protocol_sha256:
        raise ValueError('Scientific protocol changed before recovery')
    if any(p.name.isdigit() and int(p.name) != os.getpid()
           and core.engine_alive(int(p.name), workspace) for p in Path('/proc').iterdir()):
        raise ValueError('Another worker still owns this workspace')
    state = read(workspace / 'status.json')
    record = workspace / 'energy_recovery.json'
    evidence = dict(started_utc=task.now(), worker_sha256=digest(Path(__file__)),
                    previous_state=json.loads(json.dumps(state)),
                    previous_status_sha256=digest(workspace / 'status.json'))
    core.atomic_json(record, evidence)
    def interrupted(signum, frame):
        raise InterruptedError(f'Recovery interrupted by signal {signum}')
    signal.signal(signal.SIGTERM, interrupted)
    try:
        table = validate_prepared(workspace, source, baseline, parameters, state)
        os.environ.update(BLAST_DB_NR=configured_database(config) or '', RGL_USE_NULL='TRUE',
                          PYTHONNOUSERSITE='1', MUTABIND2_CONFIG=str(workspace / 'tools.runtime.json'),
                          PATH=str(workspace / 'bin') + os.pathsep + os.environ.get('PATH', ''))
        engine = task.load_engine(workspace, parameters.foldx_workers, parameters.threads,
                                  parameters.no_contributions)
        state.update(stage='recovering_energy', pid=os.getpid(), energy_recovery=str(record))
        backend = task.install_backends(engine, workspace, config, state, parameters.surface_backend,
            parameters.pb_executable, parameters.vdw_executable, parameters.pb_threads)
        reused = []
        for row in table:
            stem = row['PDBid'].lower() + '_' + row['Mutation_cleaned'].lower()
            if backend.valid_cached(workspace, stem, stem, parameters.pb_executable,
                                    parameters.pb_threads, parameters.vdw_executable):
                reused.append(stem)
        evidence.update(validated_energies=reused)
        core.atomic_json(record, evidence)
        state['recovered_energy_count'] = len(reused)
        def step(name, function, *values):
            state.update(stage=name, stage_started_utc=task.now())
            core.atomic_json(workspace / 'status.json', state)
            started = time.monotonic()
            value = function(*values)
            state.setdefault('timings_seconds', {})[name] = round(time.monotonic() - started, 3)
            core.atomic_json(workspace / 'status.json', state)
            return value
        # Recompute the inexpensive surface/contact outputs; retain frozen geometries.
        step('surface', engine.surface)
        step('reuse_wt_contacts', task.reuse_wt_contacts, engine, baseline, table)
        original = task.run_pb
        def gated(workspace, state, operation, function, *values):
            wait_for_energy(workspace, state, operation)
            return original(workspace, state, operation, function, *values)
        task.run_pb = gated
        task.finish_pb(engine, workspace, state, table, step, parameters.mutated_partner,
                       baseline=baseline, backend=backend, executable=parameters.pb_executable,
                       threads=parameters.pb_threads, vdw_executable=parameters.vdw_executable)
        task._features_and_predict(engine, workspace, state, table, step, parameters.provean_mode,
            parameters.provean_cache, parameters.no_contributions, parameters.identities)
        state.update(stage='complete', completed_utc=task.now(), predictions=len(table))
        core.atomic_json(workspace / 'status.json', state)
    except BaseException as exc:
        state.update(stage='failed', pid=os.getpid(), failed_stage=state.get('stage'), error=str(exc))
        core.atomic_json(workspace / 'status.json', state)
        raise


if __name__ == '__main__':
    main()
