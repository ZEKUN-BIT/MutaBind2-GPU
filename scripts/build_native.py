#!/usr/bin/env python3
"""Build isolated PB CUDA or compact VDW libraries from a compatible CHARMM build.

Only the output directory is modified. Input source hashes must match the
validated CHARMM 51 patch base; a different release needs porting and validation.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import shlex
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parent.parent


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (8 * 1024**3, 8 * 1024**3))
    os.nice(10)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kind', choices=['pb', 'vdw'], required=True)
    parser.add_argument('--source', type=Path, required=True, help='licensed CHARMM source root')
    parser.add_argument('--build', type=Path, required=True, help='completed CMake Makefiles build directory')
    parser.add_argument('--executable', type=Path, required=True, help='matching CHARMM ELF executable')
    parser.add_argument('--output', type=Path, required=True, help='new empty lowercase ASCII directory')
    parser.add_argument('--cuda', type=Path, default=Path('/usr/local/cuda'))
    parser.add_argument('--arch', default='sm_120', help='CUDA architecture, e.g. sm_120; match your GPU')
    parser.add_argument('--runtime-lib', type=Path, action='append', default=[], help='additional shared-library directory; repeatable')
    args = parser.parse_args()
    source, build, executable, output = (p.expanduser().resolve() for p in
                                        (args.source, args.build, args.executable, args.output))
    cuda = args.cuda.resolve()
    if not re.fullmatch(r'/[a-z0-9_./-]+', str(output)):
        parser.error('Output must be lowercase ASCII without spaces')
    if not re.fullmatch(r'(?:sm|compute)_[0-9]+[af]?', args.arch):
        parser.error('Invalid CUDA architecture')
    if output.exists() and any(output.iterdir()):
        parser.error('Output must be new or empty; failed build evidence is retained')
    pb = ROOT/'native/pb'
    reen = ROOT/'native/reen'
    vdw = ROOT/'native/vdw'
    expected = {source/'source/misc/pbeq.F90': json.loads((pb/'bridge_source_manifest.json').read_text())['original_sha256']}
    if args.kind == 'vdw':
        manifest = json.loads((vdw/'source_manifest.json').read_text())['files']
        expected = {source/'source'/directory/name: manifest[name]['original_sha256']
                    for name, directory in [('datstr.F90','util'), ('intere.F90','energy')]}
    for path, checksum in expected.items():
        if digest(path) != checksum:
            raise ValueError(f'Unsupported CHARMM patch base: {path}; port and validate this release first')
    link = shlex.split((build/'CMakeFiles/chmm.dir/link.txt').read_text())
    compiler = link[0]
    variables = dict(re.findall(r'^(Fortran_\w+) = (.*)$',
                               (build/'CMakeFiles/chmm.dir/flags.make').read_text(), re.M))
    flags = sum((shlex.split(variables[key]) for key in
                 ('Fortran_DEFINES','Fortran_INCLUDES','Fortran_FLAGS')), [])
    flags = [f for f in flags if f not in ('-O2','-O3') and not f.startswith('-J')]
    flags = ['-I'+str(output/'modules'), '-I'+str(build), '-J'+str(output/'modules'),
             *flags, '-O3', '-fno-fast-math', '-ffp-contract=off']
    runtime = [output/'lib'] + ([cuda/'lib64'] if args.kind=='pb' else []) + [p.resolve() for p in args.runtime_lib]
    for path in [executable, *runtime]:
        if not re.fullmatch(r'/[A-Za-z0-9_./-]+', str(path)):
            raise ValueError('Executable and runtime library paths must be ASCII without spaces')
    original = {str(p):digest(p) for p in [*expected, executable, build/'libchmm.so']}
    for directory in ('lib','bin','modules'):
        (output/directory).mkdir(parents=True, exist_ok=True)
    report = {'started_utc':datetime.now(timezone.utc).isoformat(), 'kind':args.kind,
              'input_hashes':original, 'commands':[], 'stage':'preparing'}

    def save():
        (output/'build_manifest.json').write_text(json.dumps(report, indent=2)+'\n')

    def run(label, argv, compile=False):
        available = int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines()
                             if x.startswith('MemAvailable:')))*1024
        if compile and available < 18 * 1024**3:
            raise RuntimeError('Build needs 18 GiB available (8 GiB allowance + 10 GiB reserve)')
        record = {'stage':label, 'argv':list(map(str,argv))}; report['commands'].append(record); save()
        started = time.monotonic()
        with (output/(label+'.log')).open('w') as log:
            # Fortran searches its current directory before -I. Compile outside
            # the original build so INTERE sees our new DATSTR module interface.
            proc = subprocess.run(argv, cwd=build if label=='link' else output, stdout=log, stderr=subprocess.STDOUT,
                                  preexec_fn=limits if compile else None)
        record.update(returncode=proc.returncode, seconds=round(time.monotonic()-started,3)); save()
        if proc.returncode:
            raise RuntimeError(f'{label} failed; see {output/(label+".log")}')

    try:
        cuda_objects = []
        if args.kind == 'pb':
            prepared = output/'pbeq.F90'
            shutil.copyfile(source/'source/misc/pbeq.F90',prepared)
            for n, patch in enumerate((pb/'pbeq_cuda_bridge.patch',reen/'reen_bridge.patch')):
                run(f'patch_{n}', ['patch','--batch','--forward','--fuzz=0',str(prepared),str(patch)])
            if digest(prepared) != json.loads((reen/'source_manifest.json').read_text())['candidate_sha256']:
                raise ValueError('Patched PB source differs from validated candidate')
            replacements = {'pbeq.F90.o':output/'pbeq.F90.o'}
            run('fortran_pb',[compiler,*flags,'-c',str(prepared),'-o',str(replacements['pbeq.F90.o'])],True)
            for native in (pb/'pbeq_cuda.cu',pb/'pb_boundary_cuda.cu',reen/'pb_reen_cuda.cu'):
                obj = output/(native.stem+'.o'); cuda_objects.append(obj)
                run(native.stem,[str(cuda/'bin/nvcc'),'--gpu-architecture='+args.arch,'--fmad=false',
                    '--ftz=false','-O3','-std=c++17','-Xcompiler=-fPIC','-Xcompiler=-ffp-contract=off',
                    '-c',str(native),'-o',str(obj)],True)
        else:
            replacements = {}
            for name in ('datstr.F90','intere.F90'):
                prepared = output/name
                directory = 'util' if name == 'datstr.F90' else 'energy'
                shutil.copyfile(source/'source'/directory/name, prepared)
                run('patch_'+name, ['patch','--batch','--forward','--fuzz=0',
                                   str(prepared),str(vdw/(name+'.patch'))])
                if digest(prepared) != manifest[name]['candidate_sha256']:
                    raise ValueError('VDW patch source checksum mismatch')
                obj = output/(name+'.o'); replacements[name+'.o'] = obj
                run(name,[compiler,*flags,'-c',str(prepared),'-o',str(obj)],True)
        library = output/'lib/libchmm.so'
        argv=[]; replaced=set()
        for i, value in enumerate(link):
            if i and link[i-1]=='-o': value=str(library)
            elif Path(value).name in replacements and not value.startswith('-'):
                replaced.add(Path(value).name); value=str(replacements[Path(value).name])
            elif value.startswith('-Wl,--dependency-file='): value='-Wl,--dependency-file='+str(output/'link.d')
            argv.append(value)
        if replaced != set(replacements): raise ValueError('Missing CHARMM objects in CMake link command')
        if args.kind=='pb': argv += [*map(str,cuda_objects),'-L'+str(cuda/'lib64'),'-lcudart','-lstdc++']
        run('link',argv,True)
        ldpath = ':'.join(map(str,runtime))
        wrapper = output/'bin/charmm'
        switch = 'MUTABIND2_PB_BACKEND=cuda' if args.kind=='pb' else 'MUTABIND2_VDW_COMPACT=1'
        wrapper.write_text(f'#!/bin/sh\nset -eu\nexport LD_LIBRARY_PATH="{ldpath}"\n'
                           f'export OMP_NUM_THREADS=1\nexport {switch}\nexec "{executable}" "$@"\n')
        wrapper.chmod(0o755)
        loaded = subprocess.run(['ldd','-r',str(executable)], env=dict(os.environ,LD_LIBRARY_PATH=ldpath),
                                text=True,capture_output=True,check=True)
        loader = loaded.stdout+loaded.stderr; (output/'ldd.txt').write_text(loader)
        if str(library) not in loader or 'not found' in loader or 'undefined symbol' in loader:
            raise ValueError('Rebuilt library failed loader/symbol validation')
        if any(digest(p)!=h for p,h in original.items()): raise ValueError('Original inputs changed during build')
        report.update(stage='complete',library_sha256=digest(library),wrapper=str(wrapper),
                      native_validation='Run CPU/GPU energy comparison before using in production')
        save(); print(json.dumps(report,indent=2))
    except BaseException as exc:
        report.update(stage='failed',error=str(exc)); save(); raise


if __name__=='__main__':
    main()
