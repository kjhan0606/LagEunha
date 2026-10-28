"""Initial pressure derivative and induced pressure-force derivative."""
import json
from pathlib import Path
import numpy as np
p=Path('/gpfs/kjhan/LagForce/Codex_pair_hg_modes_t0_20261018')
a=np.loadtxt(p/'raw_faces.dat.rhs.rank0000');a=a[np.argsort(a[:,0])]
st=np.loadtxt(p/'raw_faces.dat.hgstate.rank0000');st=st[np.argsort(st[:,0])]
f=np.loadtxt(p/'raw_faces.dat.common.rank0000');hg=np.load(p/'instantaneous_actual.npz')
assert np.array_equal(a[:,0],st[:,0]);n=len(a);gamma=float(np.float32(1.666667))
P=a[:,5];V=a[:,6];m=a[:,7];vel=a[:,3:5];r=np.linalg.norm(a[:,1:3]-3,axis=1)
comp=-gamma*P*a[:,12]/V
tcf=(gamma-1)*a[:,11]/V
hgp=(gamma-1)*hg['heat']/V
# Independent derivative from total hydro-energy and force, before HG.
Udot=a[:,10]-np.sum(vel*a[:,8:10],axis=1)
direct=(gamma-1)*Udot/V-P*a[:,12]/V
i,j=f[:,:2].astype(int).T;dx=f[:,2:4];A=f[:,4];d=np.linalg.norm(dx,axis=1);e=dx/d[:,None]
t=f[:,5:7]-.5*dx;t-=np.sum(t*e,axis=1)[:,None]*e
def force(pres):
 pair=-A[:,None]*(.5*(pres[i]+pres[j])[:,None]*e+((pres[j]-pres[i])/d)[:,None]*t)
 return np.column_stack([np.bincount(i,weights=pair[:,k],minlength=n) for k in range(2)])
parts={'compression':comp,'TCF':tcf,'HG':hgp};jerks={k:force(v)/m[:,None] for k,v in parts.items()}
out={'scope':'Instantaneous t0, fixed mass and ideal gas. Pressure-only acceleration derivative from changing pressure, evaluated at fixed geometry. Not total jerk or finite-time growth.',
 'Pdot_formula':'-gamma P Vdot/V + (gamma-1)(H_TCF+H_HG)/V',
 'hydro_Pdot_identity_max_abs_r_lt_1p9':float(np.max(abs(direct-comp-tcf)[r<1.9])),
 'initial_pressure_range_r_lt_1p9':[float(P[r<1.9].min()),float(P[r<1.9].max())],
 'regions':{}}
for name,mask in [('inner',r<.5),('edge',(r>=.5)&(r<.6)),('core',(r>=.6)&(r<1.9))]:
 rms=lambda v:float(np.sqrt(np.average(v[mask]**2,weights=m[mask])))
 out['regions'][name]={'cells':int(mask.sum()),'relative_Pdot_RMS':{k:rms(v/P) for k,v in parts.items()},
 'pressure_acceleration_derivative_RMS':{k:rms(np.linalg.norm(v,axis=1)) for k,v in jerks.items()},
 'initial_pressure_acceleration_RMS':rms(np.linalg.norm(force(P)/m[:,None],axis=1))}
with (p/'initial_pressure_seed.json').open('x') as fp:json.dump(out,fp,indent=2)
print(json.dumps(out,indent=2))
