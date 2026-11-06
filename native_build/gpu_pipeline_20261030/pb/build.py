from pathlib import Path
import hashlib,json,os,resource,shlex,subprocess,time
R=Path(__file__).resolve().parent
P=R.parents[2]
CELLS=P/'native_build/pb_cells_20261030'
OLD=Path('/home/zek/mutabind2_runs/7mkj_rpod/optimization_cuda_reen')
for d in ('lib','bin','modules'): (R/d).mkdir(exist_ok=True)
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
src=(CELLS/'pbeq_cells.F90').read_text(); start=src.index('  SUBROUTINE MAYER_MSTEP('); end=src.index('  END SUBROUTINE MAYER_MSTEP',start)+len('  END SUBROUTINE MAYER_MSTEP')
part=src[start:end]
# Independent reference uses the exact unmodified original loop.
ref=part.replace('MAYER_MSTEP','mstep_reference').replace('    use chm_kinds','').replace('    use number','').replace('    use consta','').replace('CHM_REAL4','4').replace('chm_real','8')
ref=ref.replace('    PB_DETAIL_START=PB_WALLTIME()','').replace("    CALL PB_STAGE_TIME('mstep_reference',PB_DETAIL_START)",'')
(R/'mstep_reference.f90').write_text(ref+'\n')
part=part.replace('    use chm_kinds','    use iso_c_binding, only:c_int,c_double,c_float,c_char,c_null_char\n    use psf, only:natom\n    use stream, only:outu\n    use chm_kinds',1)
interface='''    INTEGER :: CUDA_STATUS,ENV_STATUS,CUDA_J
    CHARACTER(LEN=16) :: CUDA_BACKEND,MSTEP_BACKEND
    CHARACTER(KIND=C_CHAR) :: CUDA_ERROR(256)
    INTERFACE
       FUNCTION CHARMM_PB_MSTEP_CUDA(NX,NY,NZ,NA,NS,H,RW,RI,EP,TX,TY,TZ,BX,BY,BZ, &
            QR,QK,SEL,XS,YS,ZS,RAD,G,GX,GY,GZ,ERR,ERRLEN) &
            BIND(C,NAME='charmm_pb_mstep_cuda') RESULT(STATUS)
          IMPORT C_INT,C_DOUBLE,C_FLOAT,C_CHAR
          INTEGER(C_INT),VALUE :: NX,NY,NZ,NA,NS,QR,QK,ERRLEN
          REAL(C_DOUBLE),VALUE :: H,RW,RI,EP,TX,TY,TZ,BX,BY,BZ
          INTEGER(C_INT),INTENT(IN) :: SEL(*)
          REAL(C_DOUBLE),INTENT(IN) :: XS(*),YS(*),ZS(*),RAD(*)
          REAL(C_FLOAT) :: G(*),GX(*),GY(*),GZ(*)
          CHARACTER(KIND=C_CHAR),INTENT(OUT) :: ERR(*)
          INTEGER(C_INT) :: STATUS
       END FUNCTION
    END INTERFACE
'''
part=part.replace('    REAL(KIND=8) :: PB_DETAIL_START',interface+'    REAL(KIND=8) :: PB_DETAIL_START',1)
bridge='''
    IF(CHM_REAL4.NE.C_FLOAT.OR.CHM_REAL.NE.C_DOUBLE.OR.KIND(NTPRP).NE.C_INT) ERROR STOP 101
    CUDA_BACKEND='cuda'
    CALL GET_ENVIRONMENT_VARIABLE('MUTABIND2_PB_BACKEND',CUDA_BACKEND,STATUS=ENV_STATUS)
    IF(ENV_STATUS.NE.0) CUDA_BACKEND='cuda'
    MSTEP_BACKEND='cuda'
    CALL GET_ENVIRONMENT_VARIABLE('MUTABIND2_PB_MSTEP',MSTEP_BACKEND,STATUS=ENV_STATUS)
    IF(ENV_STATUS.NE.0) MSTEP_BACKEND='cuda'
    IF(TRIM(CUDA_BACKEND).EQ.'cuda'.AND.TRIM(MSTEP_BACKEND).NE.'cpu') THEN
       CUDA_ERROR=C_NULL_CHAR
       CUDA_STATUS=CHARMM_PB_MSTEP_CUDA(NCLX,NCLY,NCLZ,NATOM,NTPRP,DCEL,RADW,RADI,EPSP, &
            TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,MERGE(1,0,QREEN),MERGE(1,0,QFKAP), &
            LSTPRP,X,Y,Z,PBRAD,RMAY,RMAYX,RMAYY,RMAYZ,CUDA_ERROR,256)
       IF(CUDA_STATUS.EQ.0) THEN
          WRITE(OUTU,'(A)') ' PBEQ_MSTEP_CUDA_SUCCESS'
          CALL PB_STAGE_TIME('MAYER_MSTEP',PB_DETAIL_START)
          RETURN
       ENDIF
       WRITE(OUTU,'(A,I0,1X,256A)') ' PBEQ_MSTEP_CUDA_STATUS ',CUDA_STATUS, &
            (CUDA_ERROR(CUDA_J),CUDA_J=1,256)
       IF(CUDA_STATUS.NE.2) ERROR STOP 104
    ENDIF
'''
part=part.replace('    PB_DETAIL_START=PB_WALLTIME()','    PB_DETAIL_START=PB_WALLTIME()'+bridge,1)
(R/'pbeq_mstep.F90').write_text(src[:start]+part+src[end:])
link=next(c['argv'] for c in json.loads((CELLS/'build_manifest.json').read_text())['commands'] if c['stage']=='charmm_link')
protected={p:sha(p) for p in link if p.endswith('.o')}
report={'commands':[],'protected_inputs':protected,'sources':{str(p):sha(p) for p in (R/'pbeq_mstep.F90',R/'mstep_cuda.cu',R/'mstep_reference.f90')}}
def limits():
    resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3));os.nice(10)
def run(stage,argv):
    available=int(next(l.split()[1] for l in Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:')))*1024
    assert available>=18*1024**3
    t=time.monotonic()
    with (R/(stage+'.log')).open('w') as f:p=subprocess.run(argv,cwd=R,preexec_fn=limits,stdout=f,stderr=subprocess.STDOUT)
    report['commands'].append(dict(stage=stage,argv=argv,returncode=p.returncode,seconds=time.monotonic()-t))
    (R/'build_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    if p.returncode:raise RuntimeError(stage)
flags=['--gpu-architecture=sm_120','--fmad=false','--ftz=false','-O3','-std=c++17','-Xcompiler=-fPIC','-Xcompiler=-ffp-contract=off']
run('cuda_compile',['/usr/local/cuda/bin/nvcc',*flags,'-c',str(R/'mstep_cuda.cu'),'-o',str(R/'mstep_cuda.o')])
run('standalone_link',['/usr/local/cuda/bin/nvcc','-shared',str(R/'mstep_cuda.o'),'-o',str(R/'libmstep.so')])
run('reference_compile',['gfortran','-O2','-fPIC','-ffp-contract=off','-shared',str(R/'mstep_reference.f90'),'-o',str(R/'libreference.so')])
ff=json.loads((OLD/'build_manifest.json').read_text())['fortran_flags']
ff=[('-J'+str(R/'modules')) if v.startswith('-J') else v for v in ff]
run('fortran_compile',['gfortran',*ff,'-I'+str(CELLS/'modules'),'-c',str(R/'pbeq_mstep.F90'),'-o',str(R/'pbeq.F90.o')])
mapping={str(CELLS/'pbeq.F90.o'):str(R/'pbeq.F90.o'),str(CELLS/'lib/libchmm.so'):str(R/'lib/libchmm.so')}
assert str(CELLS/'pbeq.F90.o') in link
link=[mapping.get(v,'-Wl,--dependency-file='+str(R/'link.d') if v.startswith('-Wl,--dependency-file=') else v) for v in link]
link.append(str(R/'mstep_cuda.o'));run('charmm_link',link)
ldpath=str(R/'lib')+':/usr/local/cuda/lib64:/home/zek/miniforge3/lib'
exe=P/'toolchain/external/charmm-51/bin/charmm'
wrapper=R/'bin/charmm';wrapper.write_text('#!/bin/sh\nset -eu\nexport LD_LIBRARY_PATH='+shlex.quote(ldpath)+'\nexport OMP_NUM_THREADS=1\nexport MUTABIND2_PB_BACKEND=cuda\nexec '+shlex.quote(str(exe))+' "$@"\n');wrapper.chmod(0o755)
p=subprocess.run(['ldd','-r',str(exe)],env=dict(os.environ,LD_LIBRARY_PATH=ldpath),capture_output=True,text=True,check=True)
(R/'ldd.txt').write_text(p.stdout+p.stderr)
assert str(R/'lib/libchmm.so') in p.stdout and 'undefined symbol' not in p.stdout+p.stderr and 'not found' not in p.stdout+p.stderr
assert all(sha(p)==h for p,h in protected.items())
report.update(stage='complete',library_sha256=sha(R/'lib/libchmm.so'))
(R/'build_manifest.json').write_text(json.dumps(report,indent=2)+'\n');print('MSTEP build complete',flush=True)
