"""Full-domain smooth row coverage and frozen-boundary work audit.

All rows are constructed; active auxiliary modes are owned at r<2.5.
Frozen support velocities are retained in z_dot and reaction work recorded.
"""
import json,sys,time,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_reversible_memory_moving_20261018')
from smooth_owner_memory import build_row
root=Path('/gpfs/kjhan/LagForce/Codex_t6_entropy_force_probes_20261018/adaptive')
state_path=root/'raw_faces.dat.hgstate.rank0000'
s=np.loadtxt(state_path);s=s[np.argsort(s[:,0])]
rhs=np.loadtxt(root/'raw_faces.dat.rhs.rank0000');rhs=rhs[np.argsort(rhs[:,0])]
x=s[:,3:5];m=s[:,5];U=s[:,8];volume=rhs[:,6];n=len(s)
assert np.array_equal(s[:,0],np.arange(n))
frozen=np.linalg.norm(x-3,axis=1)>=2.5;tree=cKDTree(x%6,boxsize=6)
rng=np.random.default_rng(20261018);v=rng.normal(size=(n,2));z=rng.normal(size=(n,2));z[frozen]=0
force=np.zeros((n,2));dz=np.zeros((n,2));cross_rows=0;periodic_rows=0;stats=[]
gamma=float(np.float32(1.666667));start=time.monotonic()
for i in range(n):
    ids,r,meta=build_row(x,i,tree);stats.append(meta)
    wrapped=np.any(abs(x[ids]-x[i])>3)
    periodic_rows+=bool(wrapped)
    if frozen[i]:continue
    # Active circular domain is safely separated from the box boundary.
    assert not wrapped
    cross_rows+=bool(np.any(frozen[ids]))
    sound=np.sqrt(gamma*(gamma-1)*U[i]/m[i]);b=np.sqrt(m[i])*sound/np.sqrt(volume[i])*r
    dz[i]=b@v[ids];np.add.at(force,ids,-b[:,None]*z[i])
    if i and i%8192==0:print('active rows reached',i,file=sys.stderr,flush=True)
aux=float(np.sum(z*dz));moving_work=float(np.sum(force[~frozen]*v[~frozen]))
frozen_work=float(np.sum(force[frozen]*v[frozen]));external_work=-frozen_work
den=float(np.sum(abs(force*v))+np.sum(abs(z*dz)))
closure=abs(moving_work+aux-external_work)/den
fnorm=float(np.linalg.norm(force,axis=1).sum());xx=x-3
torque=float(np.sum(xx[:,0]*force[:,1]-xx[:,1]*force[:,0]))
summary={key:{'min':min(row[key] for row in stats),'max':max(row[key] for row in stats)} for key in stats[0]}
assert summary['moment_max']['max']<1e-12 and closure<1e-12
assert np.linalg.norm(force.sum(axis=0))/fnorm<1e-12 and abs(torque)/(3*fnorm)<1e-12
out={'scope':__doc__,'cells':n,'active_owners':int((~frozen).sum()),'periodic_rows':periodic_rows,
     'active_rows_crossing_frozen_boundary':cross_rows,'row_summary':summary,
     'moving_kinetic_power':moving_work,'auxiliary_power':aux,'external_boundary_power':external_work,
     'boundary_work_relative_closure':closure,'total_force_relative':float(np.linalg.norm(force.sum(axis=0))/fnorm),
     'total_torque_relative':float(abs(torque)/(3*fnorm)),'seconds':time.monotonic()-start,
     'state_sha256':hashlib.sha256(state_path.read_bytes()).hexdigest(),
     'limitations':'Static full-domain Python operator. Nonzero random boundary velocities intentionally exercise external work. Production frozen-state integration, restart and parallel implementations remain unverified.'}
print(json.dumps(out,indent=2))
