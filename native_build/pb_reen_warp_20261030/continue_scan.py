"""Continue this scan with the validated backend and retain both provenances."""
import argparse
import csv
import fcntl
import json
from pathlib import Path
import subprocess
import time

INSTALL = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OLD = INSTALL / 'projects/7mkj_128_255_20261030'
NEW = INSTALL / 'projects/7mkj_128_255_20261030_optimized'
COMBINED = INSTALL / 'projects/7mkj_128_255_20261030_combined'
BACKEND = INSTALL / 'bin/charmm-pb-gpu-warp'


def read_json(path):
    return json.loads(path.read_text()) if path.exists() else {}


def read_table(path):
    if not path.exists():
        return []
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write_json(path, value):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')
    temp.replace(path)


def write_table(path, fields, rows):
    temp = path.with_suffix('.tmp')
    with temp.open('w') as stream:
        writer = csv.DictWriter(stream, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)


def collect():
    rows = read_table(OLD / 'predictions.tsv')
    expected = {rid: (chain, mutation) for rid, chain, mutation in read_json(OLD / 'scan.json')['mutation_identity']}
    if len(rows) != len(expected) or {r['Result_Id'] for r in rows} != set(expected):
        raise ValueError('Original prediction table does not cover the immutable input')
    completed = {}
    for scan in (OLD, NEW):
        protocol = read_json(scan / 'protocol.json').get('sha256', '')
        for row in read_table(scan / 'completed_predictions.tsv'):
            rid = row['Result_Id']
            if rid in completed or expected.get(rid) != (row['Chain'], row['Mutation']):
                raise ValueError('Duplicate or mismatched completed mutation: ' + rid)
            if row['Status'] != 'complete' or not protocol:
                raise ValueError('Completed prediction lacks status or protocol')
            completed[rid] = dict(row, Source_scan=str(scan), Protocol_sha256=protocol)
    merged = [completed.get(r['Result_Id'], dict(r, Source_scan='', Protocol_sha256='')) for r in rows]
    state = read_json(NEW / 'batch_status.json') or read_json(OLD / 'batch_status.json')
    if state.get('stage') == 'running' and state.get('workspace'):
        native = read_json(Path(state['workspace']) / 'status.json')
        state.update(native_stage=native.get('stage'), current_mutation=native.get('pb_current_mutation'),
                     current_energies_complete=sum(k.startswith('mutant_') and v.get('outcome') == 'complete'
                                                  for k, v in native.get('pb_runs', {}).items()))
    state.update(total=len(expected), completed=len(completed), pending=len(expected) - len(completed),
                 original_completed=len(read_table(OLD / 'completed_predictions.tsv')),
                 optimized_scan=str(NEW), combined_results=str(COMBINED / 'completed_predictions.tsv'))
    if len(completed) == len(expected):
        state['stage'] = 'complete'
    return merged, state


def aggregate():
    rows, state = collect()
    fields = list(rows[0])
    write_table(COMBINED / 'predictions.tsv', fields, rows)
    write_table(COMBINED / 'completed_predictions.tsv', fields, [r for r in rows if r['Status'] == 'complete'])
    write_json(COMBINED / 'status.json', state)
    return state


def run():
    if read_json(OLD / 'batch_status.json').get('stage') != 'paused':
        raise ValueError('Original scan must finish its current batch and pause first')
    check = read_json(HERE / 'validation/scan_mutant_comparison.json')
    if not check.get('three_phi_files_bitwise_equal') or not check.get('logged_results_identical'):
        raise ValueError('Current-scan numerical validation has not passed')
    manifest = read_json(OLD / 'scan.json')
    completed = {r['Result_Id'] for r in read_table(OLD / 'completed_predictions.tsv')}
    pending = [dict(id=rid, chain=chain, mutation=mutation, note='continued_from_original_scan')
               for rid, chain, mutation in manifest['mutation_identity'] if rid not in completed]
    input_path = COMBINED / 'remaining_mutations.tsv'
    if not (NEW / 'scan.json').exists():
        write_table(input_path, ['id', 'chain', 'mutation', 'note'], pending)
        command = [str(INSTALL / 'mutabind2'), 'scan', 'plan', '--pdb', manifest['source'],
                   '--partner1', ','.join(manifest['partner1']), '--partner2', ','.join(manifest['partner2']),
                   '--mutations', str(input_path), '--output', str(NEW),
                   '--pb-executable', str(BACKEND), '--vdw-executable', manifest['vdw_executable'],
                   '--provean-cache', manifest['provean_cache'], '--config', str(OLD / 'config.snapshot.json')]
        for key in ['fragment_mode', 'gap_threshold', 'cap_mode', 'surface_backend', 'provean_mode',
                    'max_per_batch', 'workers', 'foldx_workers', 'backend']:
            command += ['--' + key.replace('_', '-'), str(manifest[key])]
        subprocess.run(command, cwd=INSTALL, check=True)
    updated = read_json(NEW / 'scan.json')
    if updated['plan']['prepared_pdb_sha256'] != manifest['plan']['prepared_pdb_sha256']:
        raise ValueError('Prepared structure changed during continuation')
    expected = {r['id']: (r['chain'], r['mutation']) for r in pending}
    if {rid: (chain, mutation) for rid, chain, mutation in updated['mutation_identity']} != expected:
        raise ValueError('Continuation plan does not exactly cover the remaining requests')
    write_json(COMBINED / 'continuation.json', {
        'original_scan': str(OLD), 'optimized_scan': str(NEW), 'backend': str(BACKEND),
        'original_completed_ids': sorted(completed, key=int), 'remaining_count': len(pending),
        'numerical_evidence': str(HERE / 'validation/scan_mutant_comparison.json'),
        'policy': 'Retain original completed results; recompute WT with the new backend; preserve request IDs and per-result protocol hashes.'})
    command = [str(INSTALL / 'mutabind2'), 'scan', 'run', '--resume', '--output', str(NEW), '--threads', '8']
    with (NEW / 'scheduler.log').open('a') as log:
        process = subprocess.Popen(command, cwd=INSTALL, stdout=log, stderr=subprocess.STDOUT)
        while process.poll() is None:
            aggregate()
            time.sleep(10)
    state = aggregate()
    if process.returncode or state['stage'] != 'complete':
        raise RuntimeError('Optimized scan stopped; inspect its scheduler.log and status')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['run', 'status', 'aggregate'])
    action = parser.parse_args().action
    if action == 'status':
        print(json.dumps(collect()[1], ensure_ascii=False, indent=2))
    else:
        COMBINED.mkdir(exist_ok=True)
        with (COMBINED / 'continuation.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            if action == 'run':
                run()
            else:
                print(json.dumps(aggregate(), ensure_ascii=False, indent=2))
