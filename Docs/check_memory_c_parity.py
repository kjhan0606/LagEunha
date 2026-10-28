"""Full-domain C QR vs independent Python SVD row comparison, no evolution."""
import ctypes,json,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from smooth_owner_memory import build_row
here=Path(__file__).resolve().parent
lib=ctypes.CDLL(str(here/'memory_probe.so'));function=lib.rm_probe
pointer=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
function.argtypes=[ctypes.c_int,pointer,pointer,ctypes.c_int,ctypes.c_double,ctypes.c_double,pointer,pointer,pointer]
function.restype=ctypes.c_int
path=Path('/gpfs/kjhan/LagForce/Codex_t6_entropy_force_probes_20261018/adaptive/raw_faces.dat.hgstate.rank0000')
s=np.loadtxt(path);s=s[np.argsort(s[:,0])];x=s[:,3:5];tree=cKDTree(x%6,boxsize=6)
maximum=0.;rmax=0.;moment=0.;periodic=0
def call(dx,dy,coverage,owner=0):
    dx=np.ascontiguousarray(dx,dtype=np.float64);dy=np.ascontiguousarray(dy,dtype=np.float64)
    n=len(dx);out=np.zeros(n);work=np.zeros(7*n);meta=np.zeros(2)
    code=function(n,dx,dy,owner,coverage,12.,out,work,meta)
    return code,out,meta
for i in range(len(x)):
    distance,ids=tree.query(x[i]%6,k=128)
    offset=x[ids]-x[i];periodic+=bool(np.any(abs(offset)>3));offset-=6*np.rint(offset/6)
    code,actual,meta=call(offset[:,0],offset[:,1],min(distance[-1],3.))
    assert code==0,(i,code)
    ii,reference,info=build_row(x,i,tree)
    lookup={int(k):v for k,v in zip(ii,reference)}
    expected=np.array([lookup.get(int(k),0.) for k in ids])
    maximum=max(maximum,float(np.max(abs(actual-expected))))
    rmax=max(rmax,abs(meta[0]-info['radius']));moment=max(moment,float(meta[1]))
assert maximum<1e-12 and rmax<1e-12 and moment<1e-12
q=np.linspace(-1,1,64);q[0]=0
collinear=call(q,np.zeros(64),2.)[0]
uncovered=call(q,np.zeros(64),.001)[0]
bad=q.copy();bad[1]=np.nan
nonfinite=call(bad,np.zeros(64),2.)[0]
assert collinear==3 and uncovered==2 and nonfinite==1
print(json.dumps({'scope':__doc__,'rows':len(x),'periodic_rows':periodic,
                  'max_C_QR_vs_Python_SVD_coefficient_error':maximum,
                  'max_radius_difference':rmax,'max_C_polynomial_moment':moment,
                  'rejection_codes':{'collinear':collinear,'uncovered':uncovered,'nonfinite':nonfinite},
                  'header_sha256':hashlib.sha256((here/'reversible_memory.h').read_bytes()).hexdigest(),
                  'limitations':'C row function only, no force scatter, RK/restart, OpenMP, MPI or CUDA execution.'},indent=2))
