"""CUDA gradient/residual arithmetic on saved actual t6 rows, no evolution."""
import ctypes,hashlib,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
root=Path('/gpfs/kjhan/LagForce/Codex_memory_gradient_cuda_20261018')
lib=ctypes.CDLL(str(root/'probe.so'))
ptr=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
for fn in (lib.gradient_host,lib.gradient_device):
    fn.argtypes=[ctypes.c_int,ptr,ptr,ctypes.c_double,ptr,ptr,ptr];fn.restype=ctypes.c_int
def call(fn,delta,coverage):
    n=len(delta);xx=np.ascontiguousarray(delta[:,0]);yy=np.ascontiguousarray(delta[:,1])
    out=np.zeros(3*n);work=np.zeros(7*n);meta=np.zeros(2)
    code=fn(n,xx,yy,coverage,out,work,meta)
    return code,out.reshape(3,n),meta
s=np.loadtxt('/gpfs/kjhan/LagForce/Codex_t6_entropy_force_probes_20261018/adaptive/raw_faces.dat.hgstate.rank0000')
s=s[np.argsort(s[:,0])];x=s[:,3:5];tree=cKDTree(x%6,boxsize=6)
r=np.linalg.norm(x-3,axis=1);inner=np.flatnonzero(r<.6)
edge=np.flatnonzero(np.any((x<.2)|(x>5.8),axis=1))
owners=np.unique(np.r_[np.linspace(0,len(x)-1,128,dtype=int),
    inner[np.linspace(0,len(inner)-1,64,dtype=int)],edge[np.linspace(0,len(edge)-1,64,dtype=int)]])
maxcpu=0.;maxsvd=0.;moment=0.;maxradius=0.;periodic=0
for owner in owners:
    dist,ids=tree.query(x[owner]%6,k=128);delta=x[ids]-x[owner]
    periodic+=bool(np.any(abs(delta)>3));delta-=6*np.rint(delta/6)
    ca,a,ma=call(lib.gradient_host,delta,min(dist[-1],3.))
    cb,b,mb=call(lib.gradient_device,delta,min(dist[-1],3.))
    assert ca==cb==0,(int(owner),ca,cb)
    scale=np.array([1.,mb[0],mb[0]])[:,None]
    maxcpu=max(maxcpu,float(np.max(abs(a-b)*scale)))
    maxradius=max(maxradius,abs(float(ma[0]-mb[0])))
    q=delta/mb[0];P=np.column_stack((np.ones(len(ids)),q,q[:,0]**2,q[:,0]*q[:,1],q[:,1]**2))
    sqrtw=np.maximum(1-np.sum(q*q,axis=1),0)**2
    u,sv,vt=np.linalg.svd(sqrtw[:,None]*P,full_matrices=False)
    for k in range(3):
        expected=sqrtw*(u@(vt[:,k]/sv))
        if k==0:
            expected=-expected;expected[0]+=1
            target=np.zeros(6)
        else:
            target=np.zeros(6);target[k]=1
        normalized=b[k]*scale[k,0]
        maxsvd=max(maxsvd,float(np.max(abs(normalized-expected))))
        moment=max(moment,float(np.max(abs(P.T@normalized-target))))
q=np.linspace(-1.,1.,64);q[0]=0
guards={}
for name,coverage,want in [('collinear',2.,3),('uncovered',.001,2),('nonfinite',2.,1)]:
    delta=np.column_stack((q,np.zeros(64)))
    if name=='nonfinite':delta[1,0]=np.nan
    codes=[call(fn,delta,coverage)[0] for fn in (lib.gradient_host,lib.gradient_device)]
    guards[name]=codes;assert codes==[want,want]
assert max(maxcpu,maxsvd,moment,maxradius)<1e-12
out={'rows':len(owners),'periodic_rows':periodic,'max_dimensionless_CPU_GPU_error':maxcpu,
     'max_dimensionless_GPU_SVD_error':maxsvd,'max_polynomial_moment_error':moment,
     'max_radius_difference':maxradius,'guards':guards,
     'header_sha256':hashlib.sha256((root/'reversible_memory.h').read_bytes()).hexdigest(),
     'scope':'Actual CUDA execution of residual and gradient rows. One thread per row, not performance, production transport, MPI or long-time physics validation.'}
(root/'analysis.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
