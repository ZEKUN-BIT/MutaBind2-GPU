"""Bitwise tests against extracted, unmodified Fortran MAYER_MSTEP loop."""
import ctypes as C,fcntl,json,os,resource,sys
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parent;P=R.parents[2]
I=C.c_int;D=C.c_double;F=C.c_float
ip=C.POINTER(I);dp=C.POINTER(D);fp=C.POINTER(F)
gpu=C.CDLL(str(R/'libmstep.so')).charmm_pb_mstep_cuda
gpu.argtypes=[I]*5+[D]*10+[I]*2+[ip]+[dp]*4+[fp]*4+[C.c_char_p,I];gpu.restype=I
cpu=C.CDLL(str(R/'libreference.so')).mstep_reference_
cpu.argtypes=[ip,ip]+[dp]*7+[ip]*3+[dp]*7+[fp]*4+[ip]*2

def call_gpu(xyz,rads,sel,grids,dims,h,rw,ri,eps,shift,flags):
    err=C.create_string_buffer(256)
    code=gpu(*dims,len(rads),len(sel),h,rw,ri,eps,*shift,*flags,sel.ctypes.data_as(ip),
        *(a.ctypes.data_as(dp) for a in (*xyz,rads)),*(a.ctypes.data_as(fp) for a in grids),err,256)
    return code,err.value.decode()
def call_cpu(xyz,rads,sel,grids,dims,h,rw,ri,eps,shift,flags):
    ni=I(len(sel));dd=[D(a) for a in (rw,ri,eps,h,*shift)];ii=[I(a) for a in (*dims,*flags)]
    cpu(C.byref(ni),sel.ctypes.data_as(ip),*(a.ctypes.data_as(dp) for a in (*xyz,rads)),
        *(C.byref(a) for a in dd[:3]),*(C.byref(a) for a in ii[:3]),*(C.byref(a) for a in dd[3:]),
        *(a.ctypes.data_as(fp) for a in grids),*(C.byref(a) for a in ii[3:]))

def main():
    rng=np.random.default_rng(301026);cases=[]
    for k in range(96):
        dims=(int(rng.integers(1,24)),int(rng.integers(1,25)),int(rng.integers(1,23)))
        h=[.25,.5,1.,1.3][k%4]; rw=[0.,1.4,2.][k%3];ri=[0.,1.1][(k//3)%2];eps=[0.,1.,2.,4.][k%4]
        flags=(k%2,(k//2)%2);n=30
        xyz=[np.ascontiguousarray(rng.uniform(-10,35,n)) for _ in range(3)]
        rads=np.ascontiguousarray(rng.uniform(-1,4,n));sel=rng.integers(1,n+1,45,dtype=np.int32)
        shift=tuple(rng.uniform(-2,2,6));size=np.prod(dims)
        grids=[np.ascontiguousarray(rng.choice([0.,-0.,1.,-1.,eps,-eps,80.,-80.],size).astype(np.float32)) for _ in range(4)]
        cases.append((xyz,rads,sel,grids,dims,h,rw,ri,eps,shift,flags))
    # Exact radius boundaries, coincident atoms, different atom order and shells.
    for ri in (0.,1.):
        for flags in ((0,0),(0,1),(1,0),(1,1)):
            for rev in (False,True):
                xyz=[np.array([2.,2.,3.,-1.,20.],dtype=np.float64) for _ in range(3)]
                rads=np.array([1.,2.,.5,0.,1.]);sel=np.arange(1,6,dtype=np.int32)
                if rev:sel=sel[::-1].copy()
                grids=[np.full(9**3,v,dtype=np.float32) for v in (1.,80.,80.,80.)]
                cases.append((xyz,rads,sel,grids,(9,9,9),.5,1.4,ri,1.,(0.,)*6,flags))
    for k,args in enumerate(cases):
        xyz,rads,sel,base,*params=args;a=[g.copy() for g in base];b=[g.copy() for g in base]
        call_cpu(xyz,rads,sel,a,*params);code,msg=call_gpu(xyz,rads,sel,b,*params)
        assert code==0,(k,code,msg)
        for j,(old,new) in enumerate(zip(a,b)):
            diff=np.flatnonzero(old.view(np.uint32)!=new.view(np.uint32))
            assert not len(diff),(k,j,diff[:10],old[diff[:10]],new[diff[:10]])
    # Rejection must preserve every host byte. Test negative radii layer, precision,
    # zero spacing, invalid selections and overflow without dereferencing invalid data.
    args=cases[-1];xyz,rads,sel,base,dims,h,rw,ri,eps,shift,flags=args
    invalid=[dict(h=0.),dict(ri=-1.),dict(eps=1.1),dict(sel=np.array([0],dtype=np.int32)),
             dict(shift=(float('nan'),0.,0.,0.,0.,0.)),dict(dims=(INT_MAX:=2147483647,3,5))]
    values=dict(xyz=xyz,rads=rads,sel=sel,dims=dims,h=h,rw=rw,ri=ri,eps=eps,shift=shift,flags=flags)
    for change in invalid:
        grids=[g.copy() for g in base];code,msg=call_gpu(**(values|change),grids=grids)
        assert code==2,(change,code,msg)
        assert all(a.tobytes()==b.tobytes() for a,b in zip(base,grids))
    result=dict(passed=True,bitwise_cases=len(cases),transactional_rejections=len(invalid),reference='unmodified Fortran MAYER_MSTEP loop')
    (R/'test_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':
    with (P/'runtime/energy.lock').open('a') as lock:
        print('Waiting for production energy lock',flush=True);fcntl.flock(lock,fcntl.LOCK_EX)
        available=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))*1024
        assert available>=30*1024**3
        resource.setrlimit(resource.RLIMIT_AS,(20*1024**3,20*1024**3))
        main()
