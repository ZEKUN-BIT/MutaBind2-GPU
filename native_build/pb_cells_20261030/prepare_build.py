"""Build the spatial-candidate bridge without changing installed scientific code."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import resource
import shlex
import subprocess
import time

ROOT = Path(__file__).resolve().parent
INSTALL = ROOT.parents[1]
WARP = INSTALL / 'native_build/pb_reen_warp_20261030'
OLD = Path('/home/zek/mutabind2_runs/7mkj_rpod/optimization_cuda_reen')


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    for name in ['lib', 'bin', 'modules']:
        (ROOT / name).mkdir(exist_ok=True)
    source = (INSTALL / 'native_build/pb/pbeq.F90').read_text()
    start = source.index('  SUBROUTINE MAYER_REEN(')
    end = source.index('  END SUBROUTINE MAYER_REEN', start)
    part = source[start:end]
    part = part.replace('    use clcg_mod,only:clcginit',
        '    use mb_reen_cell_index,only:build_reen_cells\n    use clcg_mod,only:clcginit', 1)
    part = part.replace('    LOGICAL :: QCUDA_REEN,QAXIS_SAFE',
        '    LOGICAL :: QCUDA_REEN,QAXIS_SAFE,QCELL_INDEX\n'
        '    INTEGER,ALLOCATABLE :: CELL_OFFSETS(:),CELL_CANDIDATES(:),ATOM_CELL(:)\n'
        '    INTEGER :: CELL_J,CELL_LIMIT,CELL_BEGIN', 1)
    needle = '          CUDA_PREP_ATOMS: DO IL=1,NTPRP'
    part = part.replace(needle,
        '          CALL BUILD_REEN_CELLS(NTPRP,LSTPRP,X,Y,Z,PBRAD,RADW, &\n'
        '               CELL_OFFSETS,CELL_CANDIDATES,ATOM_CELL,QCELL_INDEX)\n'
        '          IF(QCELL_INDEX.AND.PRNLEV.GE.2) &\n'
        "               WRITE(OUTU,'(A,I0)') 'PBEQ_REEN_CELL_CANDIDATES entries=',SIZE(CELL_CANDIDATES)\n" + needle, 1)
    needle = '             DO JL=1,NTPRP\n                IF(IL.EQ.JL) CYCLE'
    assert part.count(needle) == 1
    part = part.replace(needle,
        '             CELL_LIMIT=NTPRP\n'
        '             CELL_BEGIN=1\n'
        '             IF(QCELL_INDEX) THEN\n'
        '                CELL_BEGIN=CELL_OFFSETS(ATOM_CELL(IL))\n'
        '                CELL_LIMIT=CELL_OFFSETS(ATOM_CELL(IL)+1)-CELL_BEGIN\n'
        '             ENDIF\n'
        '             DO CELL_J=1,CELL_LIMIT\n'
        '                JL=CELL_J\n'
        '                IF(QCELL_INDEX) JL=CELL_CANDIDATES(CELL_BEGIN+CELL_J-1)\n'
        '                IF(IL.EQ.JL) CYCLE', 1)
    part += ''
    needle = '    IF(ALLOCATED(CUDA_ML)) DEALLOCATE(CUDA_ML)'
    part = part.replace(needle,
        '    IF(ALLOCATED(CELL_OFFSETS)) DEALLOCATE(CELL_OFFSETS)\n'
        '    IF(ALLOCATED(CELL_CANDIDATES)) DEALLOCATE(CELL_CANDIDATES)\n'
        '    IF(ALLOCATED(ATOM_CELL)) DEALLOCATE(ATOM_CELL)\n' + needle, 1)
    candidate = ROOT / 'pbeq_cells.F90'
    candidate.write_text(source[:start] + part + source[end:])
    old = json.loads((OLD / 'build_manifest.json').read_text())
    deployed = json.loads((WARP / 'build_manifest.json').read_text())
    link = next(c['argv'] for c in deployed['commands'] if c['stage'] == 'charmm_link')
    protected = [Path(v) for v in link if v.endswith('.o')] + [WARP / 'lib/libchmm.so']
    report = dict(started_utc=datetime.now(timezone.utc).isoformat(), commands=[],
                  protected_inputs={str(p):digest(p) for p in protected},
                  source_sha256=digest(candidate), cells_sha256=digest(ROOT/'reen_cells.f90'))
    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (8*1024**3, 8*1024**3))
        os.nice(10)
    def run(stage, argv):
        available = int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines()
                             if x.startswith('MemAvailable:'))) * 1024
        if available < 18*1024**3:
            raise RuntimeError('Need 18 GiB available for a bounded build')
        started = time.monotonic()
        with (ROOT / (stage+'.log')).open('w') as log:
            result = subprocess.run(argv, cwd=ROOT, preexec_fn=limits, stdout=log,
                                    stderr=subprocess.STDOUT)
        report['commands'].append(dict(stage=stage, argv=argv, returncode=result.returncode,
                                       seconds=time.monotonic()-started))
        (ROOT/'build_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
        if result.returncode:
            raise RuntimeError('Build failed: '+stage)
    run('cells_compile', ['/usr/bin/gfortran','-O3','-fPIC','-fno-fast-math','-ffp-contract=off',
        '-J'+str(ROOT/'modules'),'-c',str(ROOT/'reen_cells.f90'),'-o',str(ROOT/'reen_cells.o')])
    flags = [('-J'+str(ROOT/'modules')) if f.startswith('-J') else f for f in old['fortran_flags']]
    run('fortran_compile', ['/usr/bin/gfortran',*flags,'-I'+str(ROOT/'modules'),'-c',
                            str(candidate),'-o',str(ROOT/'pbeq.F90.o')])
    mapping = {str(OLD/'pbeq.F90.o'):str(ROOT/'pbeq.F90.o'),
               str(WARP/'lib/libchmm.so'):str(ROOT/'lib/libchmm.so')}
    assert str(OLD/'pbeq.F90.o') in link
    argv = [mapping.get(v, '-Wl,--dependency-file='+str(ROOT/'link.d')
                        if v.startswith('-Wl,--dependency-file=') else v) for v in link]
    argv.append(str(ROOT/'reen_cells.o'))
    run('charmm_link', argv)
    executable=INSTALL/'toolchain/external/charmm-51/bin/charmm'
    ldpath=str(ROOT/'lib')+':/usr/local/cuda/lib64:/home/zek/miniforge3/lib'
    wrapper=ROOT/'bin/charmm'
    wrapper.write_text('#!/bin/sh\nset -eu\nexport LD_LIBRARY_PATH='+shlex.quote(ldpath)
        +'\nexport OMP_NUM_THREADS=1\nexport MUTABIND2_PB_BACKEND=cuda\nexec '
        +shlex.quote(str(executable))+' "$@"\n')
    wrapper.chmod(0o755)
    result=subprocess.run(['ldd','-r',str(executable)],env=dict(os.environ,LD_LIBRARY_PATH=ldpath),
                          text=True,capture_output=True,check=True)
    loader=result.stdout+result.stderr
    (ROOT/'ldd.txt').write_text(loader)
    assert str(ROOT/'lib/libchmm.so') in loader and 'not found' not in loader and 'undefined symbol' not in loader
    assert all(digest(Path(p))==h for p,h in report['protected_inputs'].items())
    report.update(stage='complete',library_sha256=digest(ROOT/'lib/libchmm.so'))
    (ROOT/'build_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(stage='complete',wrapper=str(wrapper))))


if __name__ == '__main__':
    main()
