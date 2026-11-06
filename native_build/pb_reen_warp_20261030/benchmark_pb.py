"""Sequential, resource-bounded CHARMM benchmarks on frozen coordinates.

Never stages structures or changes production templates. Each case is immutable.
"""
from pathlib import Path
import argparse
import fcntl
import hashlib
import json
import os
import re
import resource
import signal
import subprocess
import time

TASK = Path('/home/zek/Documents/Codex/2026-10-19/https-github-com-minghuilab-mutabind2-https/outputs/predictions/7mkj_rpod')
PROJECT = TASK.parents[1] / 'mutabind2'
ROOT = Path('/home/zek/mutabind2_runs/7mkj_rpod/optimization_01')
GIB = 1024**3


def available():
    return int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines()
                    if x.startswith('MemAvailable:'))) * 1024


def parse_log(path):
    text = Path(path).read_text()
    solver_text = text
    if re.search(r'^MB2_PB_STAGE AB_BEGIN$', text, re.MULTILINE):
        sections = []
        for component in ('AB', 'A', 'B'):
            matches = re.findall(r'^MB2_PB_STAGE '+component+r'_BEGIN\n(.*?)^MB2_PB_STAGE '+component+r'_END$',
                                 text, re.MULTILINE | re.DOTALL)
            if len(matches) != 1:
                raise ValueError(f'Missing or duplicate PB component {component}')
            # READ PHI also prints its cached box; compare the target SOLVE box.
            start = matches[0].index('Maximum # iterations')
            sections.append(matches[0][start:])
        solver_text = '\n'.join(sections)
    values = {}
    for key in ('VDWENER', 'GAB', 'GA', 'GB', 'DELTAGELEC'):
        matches = re.findall(r'Parameter:\s+'+key+r'\s+<-\s+"([^"]+)"', text)
        if matches:
            values[key] = float(matches[-1])
    return dict(values=values, normal=bool(re.search(r'^\s+NORMAL TERMINATION', text, re.MULTILINE)),
                converged='calculation not converged' not in text.lower(),
                iterations=list(map(int, re.findall(r'Number of iterations:\s+(\d+)', solver_text))),
                energies=list(map(float, re.findall(r'Electrostatic energy \[KCAL/MOL\]\s*=\s*([-+\d.Ee]+)', solver_text))),
                grids=[list(map(int, x)) for x in re.findall(r'NUMBER OF GRID POINTS:\s+(\d+)\s+(\d+)\s+(\d+)', solver_text)],
                boxes=re.findall(r'Box in [XYZ] from[^\n]+', solver_text))


def make_deck(preamble, resources, mode, osor=False):
    preamble = re.sub(r'^set pathpara .*$', 'set pathpara '+str(resources), preamble,
                      flags=re.MULTILINE).split('stream @pathpara/vdw_pb.str')[0]
    # CHARMM 51 ignores MAYER_REEN's local seed unless global RNG seeds are
    # explicit. Its default clock seed otherwise changes the molecular surface.
    preamble += '\nrandom clcg iseed 314159265 22222222 33333333 44444444\n'
    stream = (resources/'vdw_pb.str').read_text()
    if mode == 'vdw':
        stream = stream.split('! Inititalzation')[0] + '\nstop\n'
    elif mode == 'pb':
        stream = stream[stream.index('! Inititalzation'):]
        stream = stream.split('! molecular surface')[0] + '\nstop\n'
    if osor:
        stream = re.sub(r'^SOLVE ', 'SOLVE OSOR ', stream, flags=re.MULTILINE)
    return preamble + stream


def run_case(label, deck_text, data, stem, executable, timeout=1800, reference=None, threads=1, require_gpu=False):
    ROOT.mkdir(parents=True, exist_ok=True)
    case = ROOT/label
    case.mkdir(exist_ok=False)
    (case/'data').symlink_to(data.resolve(), target_is_directory=True)
    deck, log, timing = (case/name for name in ('input.inp', 'output.log', 'time.txt'))
    deck.write_text(deck_text)
    source_paths = [data/f'{stem}_charmm.psf', data/f'{stem}_min.crd']
    hashes = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    result = dict(label=label, executable=str(executable), source_hashes=hashes,
                  timeout_seconds=timeout, address_limit_gib=20, memory_reserve_gib=10,
                  threads=threads, require_gpu=require_gpu)
    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (20*GIB, 20*GIB))
    with (TASK/'pb_global.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if available() < 30*GIB:
            raise RuntimeError('Requires 30 GiB available memory before benchmark')
        env = os.environ.copy()
        env.update(MUTABIND2_CHARMM_THREADS=str(threads), OMP_NUM_THREADS=str(threads),
                   OMP_PROC_BIND='close', OMP_PLACES='cores')
        argv = ['/usr/bin/time', '-v', '-o', str(timing), str(executable),
                f'path={case / "data"}', f'protname={stem}']
        start = time.monotonic()
        print(json.dumps({'starting':label,'available_gib':available()/GIB}), flush=True)
        reason = None
        with deck.open('rb') as stdin, log.open('wb') as stdout:
            proc = subprocess.Popen(argv, stdin=stdin, stdout=stdout, stderr=subprocess.STDOUT,
                                    cwd=case, env=env, preexec_fn=limits, start_new_session=True)
            try:
                while proc.poll() is None:
                    if time.monotonic()-start > timeout:
                        reason = 'timeout'
                        break
                    if available() < 10*GIB:
                        reason = 'memory reserve reached'
                        break
                    time.sleep(1)
            finally:
                if proc.poll() is None:
                    os.killpg(proc.pid, signal.SIGTERM)
                    try:
                        proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(proc.pid, signal.SIGKILL)
                        proc.wait()
        result.update(returncode=proc.returncode, stopped_reason=reason, wall_seconds=time.monotonic()-start)
    result.update(parse_log(log))
    raw_log = log.read_text()
    result['gpu_solves'] = len(re.findall(r'^\s*PBEQ_CUDA_SUCCESS\b', raw_log, re.MULTILINE))
    result['stage_seconds'] = {}
    for stage, seconds in re.findall(r'PBEQ_STAGE\s+(\w+)\s+seconds=\s*([\d.Ee+-]+)', raw_log):
        result['stage_seconds'].setdefault(stage, []).append(float(seconds))
    match = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', timing.read_text())
    result['max_rss_gib'] = int(match[1])/1024**2 if match else None
    result['sources_unchanged'] = all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())
    if reference:
        ref = parse_log(reference)
        result['reference'] = str(reference)
        result['comparison'] = {
            'values_difference':{k:v-ref['values'][k] for k,v in result['values'].items() if k in ref['values']},
            'energies_difference':[a-b for a,b in zip(result['energies'],ref['energies'])],
            'grids_equal':result['grids']==ref['grids'],
            'boxes_equal':result['boxes']==ref['boxes'],
            'iterations_equal':result['iterations']==ref['iterations']}
    (case/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result), flush=True)
    missing_solves = bool(re.search(r'^PBEQ\s*$', deck_text, re.MULTILINE)) and len(result['iterations']) != 3
    if (result['returncode'] or not result['normal'] or not result['converged']
            or not result['sources_unchanged'] or missing_solves or (require_gpu and result['gpu_solves'] != 3)):
        raise RuntimeError(f'Benchmark failed: {label}; see {case}')
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--label', required=True)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--stem', required=True)
    parser.add_argument('--template', type=Path, required=True)
    parser.add_argument('--resources', type=Path, required=True)
    parser.add_argument('--mode', choices=['full','pb','vdw'], default='pb')
    parser.add_argument('--osor', action='store_true')
    parser.add_argument('--executable', type=Path, default=PROJECT/'tools/charmm')
    parser.add_argument('--reference', type=Path)
    parser.add_argument('--timeout', type=int, default=1800)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--require-gpu', action='store_true')
    args = parser.parse_args()
    deck = make_deck(args.template.read_text(), args.resources, args.mode, args.osor)
    run_case(args.label, deck, args.data.resolve(), args.stem, args.executable.resolve(),
             args.timeout, args.reference, args.threads, args.require_gpu)


if __name__ == '__main__':
    main()
