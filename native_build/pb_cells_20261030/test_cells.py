"""Test the actual compiled Fortran index for conservative, ordered candidates."""
from pathlib import Path
import ctypes
import json
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parent


def main():
    wrapper = ROOT / 'test_wrapper.f90'
    wrapper.write_text('''subroutine test_cells(n,m,selected,x,y,z,radius,water,table,counts,success) bind(C)
  use iso_c_binding
  use mb_reen_cell_index
  implicit none
  integer(c_int),value :: n,m
  integer(c_int),intent(in) :: selected(n)
  real(c_double),intent(in) :: x(m),y(m),z(m),radius(m)
  real(c_double),value :: water
  integer(c_int),intent(out) :: table(n,n),counts(n),success
  integer,allocatable :: offsets(:),candidates(:),atom_cell(:)
  logical :: ok
  integer :: i,k,c
  table=-1;counts=0
  call build_reen_cells(n,selected,x,y,z,radius,water,offsets,candidates,atom_cell,ok)
  success=merge(1,0,ok)
  if(ok) then
    do i=1,n
      c=atom_cell(i)
      k=offsets(c+1)-offsets(c)
      counts(i)=k
      table(1:k,i)=candidates(offsets(c):offsets(c+1)-1)
    enddo
  endif
end subroutine test_cells
''')
    subprocess.run(['/usr/bin/gfortran', '-shared', '-fPIC', '-fcheck=all',
                    '-I'+str(ROOT/'modules'), str(wrapper), str(ROOT/'reen_cells.o'),
                    '-o', str(ROOT/'test_cells.so')], check=True)
    routine = ctypes.CDLL(str(ROOT/'test_cells.so')).test_cells
    intarray = np.ctypeslib.ndpointer(dtype=np.int32)
    realarray = np.ctypeslib.ndpointer(dtype=np.float64)
    routine.argtypes = [ctypes.c_int, ctypes.c_int, intarray, realarray, realarray,
                        realarray, realarray, ctypes.c_double, intarray, intarray,
                        ctypes.POINTER(ctypes.c_int)]
    rng = np.random.default_rng(4011030)
    results = []

    def check(name, points, radii, water=1.4, selected=None, expect=True):
        m = len(points)
        if selected is None:
            selected = np.arange(1,m+1,dtype=np.int32)
        selected = np.asarray(selected,dtype=np.int32)
        n = len(selected)
        table = np.empty((n,n),dtype=np.int32,order='F')
        counts = np.empty(n,dtype=np.int32)
        flag = ctypes.c_int()
        routine(n,m,selected,*(np.ascontiguousarray(points[:,i],dtype=np.float64) for i in range(3)),
                np.ascontiguousarray(radii,dtype=np.float64), water, table, counts, ctypes.byref(flag))
        assert bool(flag.value)==expect, (name,'unexpected fallback/activation')
        retained = 0
        if expect:
            p = points[selected-1]
            r = radii[selected-1]+water
            for i in range(n):
                candidates = table[:counts[i],i]-1
                assert len(candidates)>0 and np.all(np.diff(candidates)>0), (name,'order')
                assert i in candidates, (name,'self')
                # Include tangent and contained spheres, even when SEGPAR would
                # eventually reject them; use a conservative tolerance here.
                distances = np.sqrt(np.sum((p-p[i])**2,axis=1))
                required = np.flatnonzero(distances <= r+r[i]+1e-9)
                assert np.all(np.isin(required,candidates)), (name,i,'missing intersection')
                retained += len(required)-1
        results.append(dict(name=name,active=bool(flag.value),atoms=n,intersections=retained))

    for seed in range(12):
        n = 256+seed*32
        pts = rng.uniform(-25,25,(n,3))
        rad = rng.uniform(.1,3,n)
        check('random_dense_'+str(seed),pts,rad)
        permutation = rng.permutation(n).astype(np.int32)+1
        check('permuted_'+str(seed),pts,rad,selected=permutation)
    pts = rng.uniform(-12,12,(1024,3))
    check('noncontiguous_selection',pts,np.full(1024,2.),selected=rng.permutation(1024)[:512]+1)
    # Cell faces, translated coordinate origins, tiny and large radii.
    for shift in (0.,-999900.,999900.):
        pts = rng.uniform(-32,15,(512,3))
        pts[:80] = np.resize(np.array([-32.,-16.,0.]),(80,3))
        pts[:40,0] = np.nextafter(pts[:40,0],np.inf)
        check('faces_shift_'+str(shift),pts+shift,np.full(512,2.6))
    for scale in (1e-4,1000.):
        check('radius_scale_'+str(scale),rng.uniform(-5,5,(512,3))*scale,
              np.full(512,scale),water=0.)
    check('over_512_intersections',rng.uniform(-.1,.1,(768,3)),np.full(768,2.))
    check('coincident_contained',np.zeros((256,3)),np.linspace(.01,3,256))
    check('small_fallback',np.zeros((255,3)),np.ones(255),expect=False)
    check('sparse_fallback',np.column_stack((np.arange(256)*10.,np.zeros((256,2)))),
          np.ones(256),water=0.,expect=False)
    check('grid_size_fallback',rng.uniform(-1e6,1e6,(256,3)),np.full(256,.01),expect=False)
    for value in (float('nan'),float('inf'),-2.,1e5):
        radii = np.ones(256);radii[5]=value
        check('invalid_radius_'+str(value),np.zeros((256,3)),radii,expect=False)
    for value in (float('nan'),float('inf'),1e7):
        pts = np.zeros((256,3));pts[5,0]=value
        check('invalid_coordinate_'+str(value),pts,np.ones(256),expect=False)
    for value in (-1.,float('nan'),float('inf')):
        check('invalid_water_'+str(value),np.zeros((256,3)),np.ones(256),water=value,expect=False)
    check('invalid_index',np.zeros((256,3)),np.ones(256),selected=np.zeros(256),expect=False)
    report = dict(passed=True,cases=len(results),tests=results)
    (ROOT/'cell_test_results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(passed=True,cases=len(results))))


if __name__=='__main__':
    main()
