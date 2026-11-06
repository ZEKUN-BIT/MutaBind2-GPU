"""Test stable integer lists against independent Python references."""
from pathlib import Path
import ctypes as C,fcntl,json,resource
import numpy as np
R=Path(__file__).resolve().parent;P=R.parents[2]
I=C.c_int;ip=C.POINTER(I);lib=C.CDLL(str(R/'audit/libpairs.so'))
dense=lib.charmm_vdw_dense_cuda;dense.argtypes=[I,I]+[ip]*5+[C.c_char_p,I]
select=lib.charmm_vdw_select_cuda;select.argtypes=[I,I]+[ip]*6+[C.c_char_p,I]
def ptr(v):return v.ctypes.data_as(ip)
def main():
    rng=np.random.default_rng(301027);count=0
    for n in (1,2,17,255,256,257,1000):
        exclusions=[];iblo=[];expected=[];rows=[]
        for i in range(1,n+1):
            row={j:(-j if rng.random()<.3 else j) for j in range(i+1,n+1) if rng.random()<.05}
            exclusions.extend(row.values());iblo.append(len(exclusions))
            expected.extend(row.get(j,j) for j in range(i+1,n+1) if row.get(j,0)<=0);rows.append(len(expected))
        ex=np.array(exclusions or [0],dtype=np.int32);iblo=np.array(iblo,dtype=np.int32)
        out=np.full(max(1,len(expected)),123,dtype=np.int32);actualrows=np.full(n,123,dtype=np.int32);total=I();err=C.create_string_buffer(256)
        code=dense(n,len(out),ptr(iblo),ptr(ex),ptr(actualrows),ptr(out),C.byref(total),err,256)
        assert code==0,(code,err.value);assert total.value==len(expected)
        assert actualrows.tolist()==rows and out[:len(expected)].tolist()==expected
        for mode in ('partition','overlap','partial','all','none'):
            sel=rng.integers(0,2,n,dtype=np.int32);other=1-sel
            if mode=='overlap':other=np.ones(n,dtype=np.int32)
            if mode=='partial':other=rng.integers(0,2,n,dtype=np.int32)
            if mode=='all':sel[:]=1;other[:]=1
            if mode=='none':sel[:]=0;other[:]=0
            wanted=[];wr=[];prev=0
            for i,end in enumerate(rows):
                wanted.extend(j for j in expected[prev:end] if (sel[i]==1 and other[abs(j)-1]==1) or (other[i]==1 and sel[abs(j)-1]==1))
                wr.append(len(wanted));prev=end
            got=np.full(max(1,len(wanted)),777,dtype=np.int32);gr=np.full(n,777,dtype=np.int32)
            code=select(n,len(got),ptr(actualrows),ptr(out),ptr(sel),ptr(other),ptr(gr),ptr(got),err,256)
            assert code==0,(mode,code,err.value)
            assert gr.tolist()==wr and got[:len(wanted)].tolist()==wanted,(n,mode)
            count+=1
        if expected:
            out.fill(777);actualrows.fill(777)
            code=dense(n,len(expected)-1,ptr(iblo),ptr(ex),ptr(actualrows),ptr(out),C.byref(total),err,256)
            assert code==3 and total.value==len(expected)
            assert np.all(out==777) and np.all(actualrows==777)
    # Cross a 64 MiB chunk boundary in both kernels, with positive exclusions
    # and negative 1-4 entries in every row; verify each row independently.
    n=6000;ex=[];ib=[];wantrows=[];total_pairs=0
    for i in range(1,n+1):
        if i+1<=n:ex.append(i+1)
        if i+2<=n:ex.append(-(i+2))
        ib.append(len(ex));total_pairs+=max(0,n-i-1);wantrows.append(total_pairs)
    out=np.empty(total_pairs,dtype=np.int32);rr=np.empty(n,dtype=np.int32);total=I();err=C.create_string_buffer(256)
    code=dense(n,total_pairs,ptr(np.array(ib,dtype=np.int32)),ptr(np.array(ex,dtype=np.int32)),ptr(rr),ptr(out),C.byref(total),err,256)
    assert code==0,(code,err.value)
    assert total.value==total_pairs and rr.tolist()==wantrows and total_pairs>16*1024**2
    sel=np.arange(n,dtype=np.int32)%2;other=1-sel
    dst=np.empty(total_pairs,dtype=np.int32);dr=np.empty(n,dtype=np.int32)
    code=select(n,total_pairs,ptr(rr),ptr(out),ptr(sel),ptr(other),ptr(dr),ptr(dst),err,256)
    assert code==0,(code,err.value)
    prev=0;dp=0
    for i,end in enumerate(wantrows,1):
        expected=np.arange(i+2,n+1,dtype=np.int32)
        if expected.size:expected[0]*=-1
        assert np.array_equal(out[prev:end],expected)
        wanted=expected[other[np.abs(expected)-1]==1] if sel[i-1] else expected[sel[np.abs(expected)-1]==1]
        assert dr[i-1]==dp+len(wanted) and np.array_equal(dst[dp:dr[i-1]],wanted)
        prev=end;dp=dr[i-1]
    # Exclusion ordering rejection must happen before writes.
    rr=np.full(3,777,dtype=np.int32);oo=np.full(3,777,dtype=np.int32)
    code=dense(3,3,ptr(np.array([2,2,2],dtype=np.int32)),ptr(np.array([3,2],dtype=np.int32)),ptr(rr),ptr(oo),C.byref(I()),C.create_string_buffer(256),256)
    assert code==2 and np.all(rr==777) and np.all(oo==777)
    report=dict(passed=True,dense_cases=8,stable_selection_cases=count+1,streaming_pairs=total_pairs,capacity_and_rejection_checks=True)
    (R/'test_results.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':
    with (P/'runtime/energy.lock').open('a') as lock:
        print('Waiting for production energy lock',flush=True);fcntl.flock(lock,fcntl.LOCK_EX)
        avail=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))*1024
        assert avail>=30*1024**3
        resource.setrlimit(resource.RLIMIT_AS,(20*1024**3,20*1024**3));main()
