"""Reference-field shear contribution to departure energy; no hydro force inferred."""
import argparse,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read,sha
from audit_canonical_pair_history import gravity_coefficient
p=argparse.ArgumentParser(description=__doc__);p.add_argument('snapshot',type=Path)
args=p.parse_args();a,t=read(args.snapshot)
ic,_=read(Path('/gpfs/kjhan/LagForce/Codex_canonical_initial_20261017/run/khout.000000.dat'))
maskid=np.uint64(72057594037927935)
assert np.array_equal(a['id']&maskid,ic['id']&maskid)
x=np.column_stack((a['x']-3,a['y']-3));r=np.linalg.norm(x,axis=1)
omega=np.sqrt(gravity_coefficient(r));q=r/.5
M=omega**2*r**3
dm=np.zeros_like(r)
lo=q<=.5;hi=(q>.5)&(q<1)
dm[lo]=32*q[lo]**2-192*q[lo]**4+192*q[lo]**5
dm[hi]=64*q[hi]**2-192*q[hi]**3+192*q[hi]**4-64*q[hi]**5
good=r>1e-12
dlog=np.zeros_like(r);dlog[good]=.5*(q[good]*dm[good]/M[good]-3)
vc=omega[:,None]*np.column_stack((-x[:,1],x[:,0]))
u=np.column_stack((a['vx'],a['vy']))-vc
ur=np.zeros_like(r);ut=np.zeros_like(r)
ur[good]=np.sum(u[good]*x[good],axis=1)/r[good]
ut[good]=(x[good,0]*u[good,1]-x[good,1]*u[good,0])/r[good]
work=-a['mass']*omega*dlog*ur*ut
assert np.isfinite(work).all()
r0=np.linalg.norm(np.column_stack((ic['x']-3,ic['y']-3)),axis=1)
masks={'all_active':r<2.5,'inner_background':(r0<.5)&(ic['den']<.5),
 'inner_shear_band':(r0>=.3)&(r0<.5),'inner_disk_edge':(r0<.6)&(ic['den']>.5),
 'disk_core':(r0>=.6)&(r0<1.9)&(ic['den']>.5)}
out={'snapshot':str(args.snapshot),'snapshot_sha256':sha(args.snapshot),'time':t,
 'regions':{name:{'cells':int(m.sum()),'shear_power':float(work[m&(r<2.5)].sum())} for name,m in masks.items()},
 'scope':'Fixed-mass canonical spline potential, rc=.5, center(3,3), no steering. '
 'For u=v-vc, d(.5*m*u^2)/dt=u·Fhydro+u·Faux-m*r*Omega_prime*ur*ut. '
 'Only the final shear term is measured here. Hydro work and HG split remain separate; not a full derivative or stability proof.'}
print(json.dumps(out,indent=2))
