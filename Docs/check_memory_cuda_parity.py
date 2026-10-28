"""Execute the real quadratic-row header on CUDA, compare host QR and Python SVD.

Selected actual t6 interior/boundary rows and invalid-input guards. This checks
device arithmetic only, not production GPU force scatter or MPI evolution.
"""
import argparse,ctypes,json,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from smooth_owner_memory import build_row
root=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--previous-library',type=Path)
args=parser.parse_args()
lib=ctypes.CDLL(str(root/'memory_cuda_probe.so'))
ptr=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
previous=ctypes.CDLL(str(args.previous_library)) if args.previous_library else None
functions=[lib.rm_host,lib.rm_device]+([previous.rm_host] if previous else [])
for fn in functions:
    fn.argtypes=[ctypes.c_int,ptr,ptr,ctypes.c_int,ctypes.c_double,ptr,ptr,ptr]
    fn.restype=ctypes.c_int
def call(fn,dx,dy,coverage):
    dx=np.ascontiguousarray(dx);dy=np.ascontiguousarray(dy);n=len(dx)
    out=np.zeros(n);work=np.zeros(7*n);meta=np.zeros(2)
    code=fn(n,dx,dy,0,coverage,out,work,meta)
    return code,out,meta
s=np.loadtxt('/gpfs/kjhan/LagForce/Codex_t6_entropy_force_probes_20261018/adaptive/raw_faces.dat.hgstate.rank0000')
s=s[np.argsort(s[:,0])];x=s[:,3:5];tree=cKDTree(x%6,boxsize=6)
r=np.linalg.norm(x-3,axis=1)
interior=np.flatnonzero(r<.6);edge=np.flatnonzero(np.any((x<.2)|(x>5.8),axis=1))
owners=np.unique(np.concatenate((np.linspace(0,len(x)-1,128,dtype=int),
    interior[np.linspace(0,len(interior)-1,64,dtype=int)],edge[np.linspace(0,len(edge)-1,64,dtype=int)])))
maxdiff=0.;maxsvd=0.;maxradius=0.;maxmoment=0.;periodic=0;previous_equal=True
for i in owners:
    distance,ids=tree.query(x[i]%6,k=128);delta=x[ids]-x[i]
    periodic+=bool(np.any(abs(delta)>3));delta-=6*np.rint(delta/6)
    cpu,a,ma=call(lib.rm_host,delta[:,0],delta[:,1],min(distance[-1],3.))
    gpu,b,mb=call(lib.rm_device,delta[:,0],delta[:,1],min(distance[-1],3.))
    assert cpu==gpu==0,(int(i),cpu,gpu)
    if previous:
        oldcode,oldrow,oldmeta=call(previous.rm_host,delta[:,0],delta[:,1],min(distance[-1],3.))
        previous_equal &= oldcode==cpu and np.array_equal(oldrow,a) and np.array_equal(oldmeta,ma)
    ii,ref,_=build_row(x,i,tree);lookup=dict(zip(ii,ref));expected=np.array([lookup.get(j,0.) for j in ids])
    maxdiff=max(maxdiff,float(np.max(abs(a-b))));maxsvd=max(maxsvd,float(np.max(abs(b-expected))))
    maxradius=max(maxradius,float(abs(ma[0]-mb[0])));maxmoment=max(maxmoment,float(mb[1]))
q=np.linspace(-1.,1.,64);q[0]=0.;nan=q.copy();nan[1]=np.nan
guards={}
for name,dx,dy,cov,want in [('collinear',q,np.zeros(64),2.,3),
    ('uncovered',q,np.zeros(64),.001,2),('nonfinite',nan,np.zeros(64),2.,1)]:
    codes=[call(fn,dx,dy,cov)[0] for fn in (lib.rm_host,lib.rm_device)]
    guards[name]=codes;assert codes==[want,want]
out={'scope':__doc__,'actual_rows':len(owners),'periodic_candidate_rows':periodic,
    'max_CPU_GPU_coefficient_difference':maxdiff,'max_GPU_SVD_coefficient_difference':maxsvd,
    'max_CPU_GPU_radius_difference':maxradius,'max_GPU_moment':maxmoment,
    'host_device_guard_codes':guards,
    'previous_host_rows_and_meta_bitwise_equal':bool(previous_equal) if previous else None,
    'header_sha256':hashlib.sha256((root/'reversible_memory.h').read_bytes()).hexdigest(),
    'limitations':'Single CUDA thread per tested row, no performance claim. No production force assembly, moving stage, MPI or long-time CUDA trajectory validation.'}
print(json.dumps(out,indent=2))
assert max(maxdiff,maxsvd,maxradius,maxmoment)<1e-12
assert previous_equal
