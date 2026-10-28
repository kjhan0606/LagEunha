"""Static full active-domain Bv/-B^Tz backend and frozen-work comparison.

Independent SVD rows are shared across backends to isolate assembly errors.
This does not validate production neighbor search, row build, MPI or RK steps.
"""
import ctypes,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from smooth_owner_memory import build_row
here=Path(__file__).resolve().parent
root=Path('/gpfs/kjhan/LagForce/Codex_t6_entropy_force_probes_20261018/adaptive')
s=np.loadtxt(root/'raw_faces.dat.hgstate.rank0000');s=s[np.argsort(s[:,0])]
rhs=np.loadtxt(root/'raw_faces.dat.rhs.rank0000');rhs=rhs[np.argsort(rhs[:,0])]
x=s[:,3:5];m=s[:,5];U=s[:,8];volume=rhs[:,6];n=len(s)
assert np.array_equal(s[:,0],rhs[:,0]) and np.array_equal(s[:,0],np.arange(n))
frozen=np.linalg.norm(x-3,axis=1)>=2.5;owners=np.flatnonzero(~frozen);nr=len(owners)
cache=here/'operator.npz'
if cache.exists():
    data=np.load(cache);assert np.array_equal(data['owners'],owners)
    ptr=data['ptr'];col=data['col'];b=data['b']
else:
    tree=cKDTree(x%6,boxsize=6);ptr=[0];col=[];b=[];gamma=float(np.float32(1.666667))
    for i in owners:
        ids,row,meta=build_row(x,i,tree)
        assert not np.any(abs(x[ids]-x[i])>3) # torque uses unwrapped active support
        scale=np.sqrt(m[i])*np.sqrt(gamma*(gamma-1)*U[i]/m[i])/np.sqrt(volume[i])
        col.extend(ids);b.extend(scale*row);ptr.append(len(col))
    ptr=np.array(ptr,dtype=np.int32);col=np.array(col,dtype=np.int32);b=np.array(b)
    np.savez(cache,owners=owners,ptr=ptr,col=col,b=b)
assert ptr[0]==0 and ptr[-1]==len(col)==len(b) and np.all(np.diff(ptr)>0)
assert np.all((col>=0)&(col<n)) and np.isfinite(b).all()
B=csr_matrix((b,col,ptr),shape=(nr,n))
lib=ctypes.CDLL(str(here/'memory_assembly.so'));fn=lib.memory_assemble
ip=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
dp=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
fn.argtypes=[ctypes.c_int,ctypes.c_int,ip,ip,dp,dp,dp,dp,dp,ctypes.c_int];fn.restype=ctypes.c_int
rng=np.random.default_rng(20261018);z=np.ascontiguousarray(rng.normal(size=(nr,2)))
cases={}
for case,v in [('random_with_moving_boundary',rng.normal(size=(n,2))),('actual_saved_velocity',s[:,6:8])]:
    v=np.ascontiguousarray(v);expected_dz=B@v;expected_force=-(B.T@z)
    rows=[]
    for mode in [0,1,2,2,2]:
        dz=np.zeros((nr,2));force=np.zeros((n,2))
        assert fn(n,nr,ptr,col,b,v,z,dz,force,mode)==0
        assert np.isfinite(dz).all() and np.isfinite(force).all()
        ed=float(np.max(abs(dz-expected_dz))/max(np.max(abs(expected_dz)),1e-300))
        ef=float(np.max(abs(force-expected_force))/max(np.max(abs(expected_force)),1e-300))
        aux=float(np.sum(z*dz));free=float(np.sum(force[~frozen]*v[~frozen]));external=-float(np.sum(force[frozen]*v[frozen]))
        den=float(np.sum(abs(force*v))+np.sum(abs(z*dz)))
        closure=abs(free+aux-external)/den
        fnorm=float(np.linalg.norm(force,axis=1).sum());xx=x-3
        torque=float(np.sum(xx[:,0]*force[:,1]-xx[:,1]*force[:,0]))
        row={'backend':['serial','openmp4','cuda_atomic'][mode],
             'dz_relative_max_error':ed,'force_relative_max_error':ef,
             'work_relative_closure':closure,'net_force_relative':float(np.linalg.norm(force.sum(axis=0))/fnorm),
             'net_torque_relative':abs(torque)/(3*fnorm),'external_boundary_power':external}
        assert max(ed,ef,closure,row['net_force_relative'],row['net_torque_relative'])<1e-12
        rows.append(row)
    cases[case]=rows
print(json.dumps({'scope':__doc__,'particles':n,'active_rows':nr,'nonzero_coefficients':len(b),
    'cases':cases,'limitations':'A standalone assembler; live simulator remains host assembly. CPU/GPU share precomputed SVD coefficients. Separate row parity evidence is required. No production CUDA integration, MPI exchange or long-time stability claim.'},indent=2))
