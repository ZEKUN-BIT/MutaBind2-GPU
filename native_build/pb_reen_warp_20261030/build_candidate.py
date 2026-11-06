"""Build an isolated REEN kernel while reusing the validated Fortran/SOR objects."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import resource
import shlex
import subprocess
import time

ROOT = Path(__file__).resolve().parent
ORIGINAL = Path('/home/zek/mutabind2_runs/7mkj_rpod/optimization_cuda_reen')
previous = json.loads((ORIGINAL / 'build_manifest.json').read_text())
report = {'started_utc': datetime.now(timezone.utc).isoformat(), 'commands': []}


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (8 * 1024**3, 8 * 1024**3))
    os.nice(10)


def run(stage, argv):
    available = int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                         if line.startswith('MemAvailable:'))) * 1024
    if available < 18 * 1024**3:
        raise RuntimeError('Insufficient available memory for bounded build')
    started = time.monotonic()
    with (ROOT / (stage + '.log')).open('w') as log:
        result = subprocess.run(argv, cwd=ORIGINAL, stdout=log, stderr=subprocess.STDOUT,
                                preexec_fn=limits)
    report['commands'].append({'stage': stage, 'argv': list(map(str, argv)),
                               'returncode': result.returncode,
                               'seconds': time.monotonic() - started})
    (ROOT / 'build_manifest.json').write_text(json.dumps(report, indent=2) + '\n')
    if result.returncode:
        raise RuntimeError(f'{stage} failed; inspect {ROOT / (stage + ".log")}')


old_link = next(record['argv'] for record in reversed(previous['commands'])
                if record['stage'] == 'charmm_link')
inputs = [Path(value) for value in old_link if value.endswith('.o')]
inputs.append(ORIGINAL / 'lib/libchmm.so')
report['reused_input_sha256'] = {str(path): digest(path) for path in inputs}
obj = ROOT / 'pb_reen_cuda.o'
library = ROOT / 'lib/libchmm.so'
run('cuda_compile', ['/usr/local/cuda/bin/nvcc', *previous['cuda_flags'], '-c',
                     str(ROOT / 'pb_reen_cuda.cu'), '-o', str(obj)])
run('standalone_link', ['/usr/local/cuda/bin/nvcc', '--shared', '--cudart=shared',
                        str(obj), '-o', str(ROOT / 'lib/libreen.so')])
mapping = {str(ORIGINAL / 'pb_reen_cuda.o'): str(obj),
           str(ORIGINAL / 'lib/libchmm.so'): str(library)}
argv = [mapping.get(value, '-Wl,--dependency-file=' + str(ROOT / 'link.d')
                    if value.startswith('-Wl,--dependency-file=') else value)
        for value in old_link]
if str(obj) not in argv or str(library) not in argv:
    raise ValueError('Candidate object or library missing from link command')
run('charmm_link', argv)
executable = Path('/home/zek/mutabind2-tool/toolchain/external/charmm-51/bin/charmm')
ldpath = str(ROOT / 'lib') + ':/usr/local/cuda/lib64:/home/zek/miniforge3/lib'
wrapper = ROOT / 'bin/charmm'
wrapper.write_text('#!/bin/sh\nset -eu\nexport LD_LIBRARY_PATH=' + shlex.quote(ldpath)
                   + '\nexport OMP_NUM_THREADS=1\nexport MUTABIND2_PB_BACKEND=cuda\nexec '
                   + shlex.quote(str(executable)) + ' "$@"\n')
wrapper.chmod(0o755)
loaded = subprocess.run(['ldd', '-r', str(executable)], text=True, capture_output=True,
                        env=dict(os.environ, LD_LIBRARY_PATH=ldpath), check=True)
loader = loaded.stdout + loaded.stderr
(ROOT / 'ldd.txt').write_text(loader)
if str(library) not in loader or 'not found' in loader or 'undefined symbol' in loader:
    raise ValueError('Candidate shared library failed loader validation')
if any(digest(path) != checksum for path, checksum in report['reused_input_sha256'].items()):
    raise ValueError('A reused build input changed')
report.update(stage='complete', source_sha256=digest(ROOT / 'pb_reen_cuda.cu'),
              library_sha256=digest(library), wrapper=str(wrapper))
(ROOT / 'build_manifest.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({key: report[key] for key in ('stage', 'library_sha256', 'wrapper')}, indent=2))
