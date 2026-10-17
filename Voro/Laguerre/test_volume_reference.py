import ctypes,json
import numpy as np
from scipy.spatial import HalfspaceIntersection,ConvexHull
from pathlib import Path
lib=ctypes.CDLL(str(Path('libgate.so').resolve()));array=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');lib.cell.argtypes=[ctypes.c_int,ctypes.c_double,array,array,ctypes.POINTER(ctypes.c_double)];lib.cell.restype=ctypes.c_int
rng=np.random.default_rng(20261017);rows=[]
for case in range(64):
    n=(8,12,20,32)[case%4];direction=rng.normal(size=(n,3));direction/=np.linalg.norm(direction,axis=1)[:,None];xyz=direction*rng.uniform(.8,1.5,size=(n,1));w=rng.uniform(0,.1,size=n);w0=.04
    data=np.ascontiguousarray(np.column_stack([xyz,w]));output=np.zeros(1536);volume=ctypes.c_double();count=lib.cell(n,w0,data,output,ctypes.byref(volume))
    assert count>3,(case,count)
    bound=.5*(np.sum(xyz*xyz,axis=1)+w0-w);half=np.vstack([np.column_stack([xyz,-bound]),np.column_stack([np.vstack([np.eye(3),-np.eye(3)]),np.full(6,-2.)])])
    reference=HalfspaceIntersection(half,np.zeros(3)).intersections;expected=ConvexHull(reference).volume
    error=abs(volume.value-expected)/expected;violation=float(np.max(output[:3*count].reshape(-1,3)@half[:,:3].T+half[:,3]));assert error<1e-10 and violation<1e-10,(case,error,violation)
    data[:,3]+=2;shifted=ctypes.c_double();nc=lib.cell(n,w0+2,data,output,ctypes.byref(shifted));assert nc>3 and abs(shifted.value-volume.value)<1e-10
    permuted=np.ascontiguousarray(data[rng.permutation(n)]);nv=lib.cell(n,w0+2,permuted,output,ctypes.byref(shifted));assert nv>3 and abs(shifted.value-volume.value)<1e-10
    rows.append(dict(case=case,neighbors=n,vertices=count,relative_volume_error=error,max_halfspace_violation=violation))
Path('geometry_reference_report.json').write_text(json.dumps(dict(passed=True,cases=rows,max_relative_volume_error=max(r['relative_volume_error'] for r in rows),scope='3D convex geometry only; no fluid evolution. Non-simple cut-through-vertex configurations explicitly return -1.'),indent=2))
print('PASS',len(rows),'random weighted cells, gauge and permutation checks')
