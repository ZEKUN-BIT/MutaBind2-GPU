from pathlib import Path
import argparse,hashlib,json,os,re,resource,shlex,subprocess,time
R=Path(__file__).resolve().parent;P=R.parents[2];SRC=Path('/home/zek/下载/charmm/source');BUILD=Path('/home/zek/下载/charmm/build/cmake')
a=argparse.ArgumentParser();a.add_argument('--audit',action='store_true');args=a.parse_args()
OUT=R/('audit' if args.audit else 'production')
for d in ('lib','bin','modules'):(OUT/d).mkdir(parents=True,exist_ok=True)
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
# Original NBONDA continues to establish all metadata and group geometry. Only
# the pair loop is bypassed after conservative, explicit domain checks.
src=(SRC/'nbonds/nbonda.F90').read_text()
src=src.replace('  use chm_kinds','  use iso_c_binding, only:c_int,c_char\n  use ieee_arithmetic, only:ieee_is_finite\n  use chm_kinds',1)
interface='''  logical :: mb_dense,mb_done
  integer :: mb_status,mb_env_status,mb_env_length
  character(len=16) :: mb_env
  character(kind=c_char) :: mb_error(256)
  real(chm_real) :: mb_diameter2
  interface
     function mb_dense_cuda(n,capacity,iblo,inb,rows,pairs,count,err,elen) &
          bind(C,name='charmm_vdw_dense_cuda') result(status)
       import c_int,c_char
       integer(c_int),value :: n,capacity,elen
       integer(c_int),intent(in) :: iblo(*),inb(*)
       integer(c_int) :: rows(*),pairs(*),count
       character(kind=c_char) :: err(*)
       integer(c_int) :: status
     end function
  end interface
'''
src=src.replace('  implicit none','  implicit none\n'+interface,1)
start=src.index('#if KEY_EXPAND == 1 && KEY_FOURD == 1',src.index('  NGAT=0'))
end=src.index('!=======================================================================',start)
bridge='''  call get_environment_variable('MUTABIND2_VDW_GPU',mb_env,length=mb_env_length,status=mb_env_status)
  mb_dense=mb_env_status==0.and.mb_env_length==1.and.mb_env=='1'
  mb_done=.false.
  if(mb_dense) then
     mb_dense=NATOM>0.and.NBTHOL==0.and..not.LEXS.and..not.QGMREM
     mb_dense=mb_dense.and.all(IMOVE(1:NATOM)<=0).and.all(IMOVEG(1:NGRP)<=0)
     mb_dense=mb_dense.and..not.any(ISDRUDE(1:NATOM))
#if KEY_MTS==1
     mb_dense=mb_dense.and..not.QTBMTS
#endif
#if KEY_PERT==1
     mb_dense=mb_dense.and..not.QPERT
#endif
#if KEY_TSM==1
     mb_dense=mb_dense.and..not.LTSM
#endif
#if KEY_FOURD==1
     mb_dense=mb_dense.and..not.DIM4
#endif
#if KEY_REPLICA==1
     mb_dense=mb_dense.and..not.QREP
#endif
#if KEY_PBOUND==1
     mb_dense=mb_dense.and..not.QBOUN
#endif
#if KEY_PARALLEL==1
     mb_dense=mb_dense.and.NUMNOD==1
#endif
     mb_dense=mb_dense.and.all(ieee_is_finite(X(1:NATOM))).and. &
          all(ieee_is_finite(Y(1:NATOM))).and.all(ieee_is_finite(Z(1:NATOM)))
     ! A factor-of-two margin ensures every group pair has RSDISP=1,
     ! including rounding in group centers. Short-cutoff cases use NBONDA.
     mb_dense=mb_dense.and.maxval(abs(X(1:NATOM)))<1.0e6_chm_real.and. &
          maxval(abs(Y(1:NATOM)))<1.0e6_chm_real.and.maxval(abs(Z(1:NATOM)))<1.0e6_chm_real
     mb_diameter2=(maxval(X(1:NATOM))-minval(X(1:NATOM)))**2+ &
          (maxval(Y(1:NATOM))-minval(Y(1:NATOM)))**2+ &
          (maxval(Z(1:NATOM))-minval(Z(1:NATOM)))**2
     mb_dense=mb_dense.and.ieee_is_finite(mb_diameter2).and.ieee_is_finite(CTNBSQ).and. &
          CUTNB>0.and.mb_diameter2<0.25_chm_real*CTNBSQ
     if(mb_dense) then
        mb_status=mb_dense_cuda(NATOM,MAXJNB,IBLO14,INB14,INBLO,JNB,NNNB,mb_error,256)
        if(mb_status==3) return ! Preserve the original allocation retry contract.
        if(mb_status==0) then
           mb_done=.true.
           RSDISP(1:NGRP)=1
           if(PRNLEV>=2) write(OUTU,'(A,I0)') ' MBVDW_DENSE_CUDA_SUCCESS pairs=',NNNB
        else if(mb_status==2) then
           NNNB=0 ! The original builder overwrites any declined candidate list.
        else
           write(OUTU,'(A,256A)') ' MBVDW_DENSE_CUDA_ERROR ',mb_error
           error stop 105
        endif
     endif
  endif
  if(.not.mb_done) then
'''
src=src[:start]+bridge+src[start:end]+'  endif\n'+src[end:]
(R/'nbonda.F90').write_text(src)
src=(P/'toolchain/external/charmm-vdw-compact/source/intere.F90').read_text()
start=src.index('SUBROUTINE INTER2(');end=src.index('END SUBROUTINE INTER2',start)
part=src[start:end]
part=part.replace('  use chm_kinds','  use iso_c_binding, only:c_int,c_char\n  use chm_kinds',1)
interface='''  integer :: mb_gpu_status,mb_gpu_env_status
  logical :: mb_gpu_selected
  character(len=16) :: mb_gpu_env
  character(kind=c_char) :: mb_gpu_error(256)
  interface
     function mb_select_cuda(n,capacity,rows,pairs,sel,other,drows,dest,err,elen) &
          bind(C,name='charmm_vdw_select_cuda') result(status)
       import c_int,c_char
       integer(c_int),value :: n,capacity,elen
       integer(c_int),intent(in) :: rows(*),pairs(*),sel(*),other(*)
       integer(c_int) :: drows(*),dest(*)
       character(kind=c_char) :: err(*)
       integer(c_int) :: status
     end function
  end interface
'''
part=part.replace('  implicit none','  implicit none\n'+interface,1)
start2=part.index('     ! Inter3g added by Lazaridis')
end2=part.index('#if KEY_GRID==1\n  EndIf',start2)
bridge='''     mb_gpu_selected=.false.
     call get_environment_variable('MUTABIND2_VDW_GPU',mb_gpu_env,status=mb_gpu_env_status)
     if(mb_compact.and.mb_gpu_env_status==0.and.mb_gpu_env=='1') then
        mb_gpu_status=mb_select_cuda(NATOM,size(bnbndc%JNB),bnbnd%INBLO,bnbnd%JNB, &
             ISLCT,JSLCT,bnbndc%INBLO,bnbndc%JNB,mb_gpu_error,256)
        if(mb_gpu_status==0) then
           mb_gpu_selected=.true.
           bnbndc%INBLOG(1:NGRP)=0
           write(OUTU,'(A,I0)') ' MBVDW_SELECT_CUDA_SUCCESS pairs=',bnbndc%INBLO(NATOM)
        else
           write(OUTU,'(A,I0,1X,256A)') ' MBVDW_SELECT_CUDA_ERROR ',mb_gpu_status,mb_gpu_error
           error stop 106
        endif
     endif
     if(.not.mb_gpu_selected) then
'''
part=part[:start2]+bridge+part[start2:end2]+'     endif\n'+part[end2:]
(R/'intere.F90').write_text(src[:start]+part+src[end:])
(R/'datstr.F90').write_bytes((P/'toolchain/external/charmm-vdw-compact/source/datstr.F90').read_bytes())
report={'commands':[],'audit':args.audit,'sources':{str(p):sha(p) for p in R.glob('*') if p.suffix in ('.cu','.F90')}}
def limits():
    resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3));os.nice(10)
def run(stage,argv):
    available=int(next(l.split()[1] for l in Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:')))*1024
    assert available>=18*1024**3
    t=time.monotonic()
    with (OUT/(stage+'.log')).open('w') as f:p=subprocess.run(argv,cwd=OUT,preexec_fn=limits,stdout=f,stderr=subprocess.STDOUT)
    report['commands'].append(dict(stage=stage,argv=argv,returncode=p.returncode,seconds=time.monotonic()-t))
    (OUT/'build_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    if p.returncode:raise RuntimeError(stage)
run('cuda_compile',['/usr/local/cuda/bin/nvcc','--gpu-architecture=sm_120','--fmad=false','-O3','-std=c++17','-Xcompiler=-fPIC','-c',str(R/'pairs_cuda.cu'),'-o',str(OUT/'pairs_cuda.o')])
run('standalone_link',['/usr/local/cuda/bin/nvcc','-shared',str(OUT/'pairs_cuda.o'),'-o',str(OUT/'libpairs.so')])
fm=(BUILD/'CMakeFiles/chmm.dir/flags.make').read_text();flags=[]
for key in ('Fortran_DEFINES','Fortran_INCLUDES','Fortran_FLAGS'):flags+=shlex.split(re.search(r'^'+key+r' = (.*)$',fm,re.M)[1])
flags=['-I'+str(OUT/'modules'),'-I'+str(BUILD),'-J'+str(OUT/'modules'),'-I'+str(SRC/'nbonds'),*flags]
if args.audit:flags.append('-DMUTABIND2_VDW_AUDIT')
for name in ('datstr.F90','intere.F90','nbonda.F90'):
    run('compile_'+name,['gfortran',*flags,'-c',str(R/name),'-o',str(OUT/(name+'.o'))])
link=shlex.split((BUILD/'CMakeFiles/chmm.dir/link.txt').read_text());protected={}
for i,v in enumerate(link):
    if v in ('CMakeFiles/chmm.dir/datstr.F90.o','CMakeFiles/chmm.dir/intere.F90.o','CMakeFiles/chmm.dir/nbonda.F90.o'):link[i]=str(OUT/Path(v).name)
    elif v.startswith('CMakeFiles/') or v in ('libcharmm_c.a','libcharmm_cxx.a'):link[i]=str(BUILD/v);protected[link[i]]=sha(link[i])
    elif v=='libchmm.so':link[i]=str(OUT/'lib/libchmm.so')
    elif v.startswith('-Wl,--dependency-file='):link[i]='-Wl,--dependency-file='+str(OUT/'link.d')
link += [str(OUT/'pairs_cuda.o'),'-L/usr/local/cuda/lib64','-lcudart','-lstdc++']
if args.audit:link.append('/usr/lib/x86_64-linux-gnu/libcrypto.so.3')
report['protected_inputs']=protected;run('link',link)
ldpath=str(OUT/'lib')+':/usr/local/cuda/lib64:/home/zek/miniforge3/lib';exe=P/'toolchain/external/charmm-51/bin/charmm'
for mode in ('gpu','cpu'):
    wrapper=OUT/'bin'/('charmm_'+mode)
    wrapper.write_text('#!/bin/sh\nset -eu\nexport LD_LIBRARY_PATH='+shlex.quote(ldpath)+'\nexport OMP_NUM_THREADS=1\nexport MUTABIND2_VDW_COMPACT=1\nexport MUTABIND2_VDW_GPU='+('1' if mode=='gpu' else '0')+'\nexec '+shlex.quote(str(exe))+' "$@"\n');wrapper.chmod(0o755)
p=subprocess.run(['ldd','-r',str(exe)],env=dict(os.environ,LD_LIBRARY_PATH=ldpath),capture_output=True,text=True,check=True)
(OUT/'ldd.txt').write_text(p.stdout+p.stderr)
assert str(OUT/'lib/libchmm.so') in p.stdout and 'undefined symbol' not in p.stdout+p.stderr and 'not found' not in p.stdout+p.stderr
assert all(sha(p)==h for p,h in protected.items())
report.update(stage='complete',library_sha256=sha(OUT/'lib/libchmm.so'))
(OUT/'build_manifest.json').write_text(json.dumps(report,indent=2)+'\n');print('VDW build complete',args.audit,flush=True)
