"""Uniform-density tangent check of weight-swept mass transport at W=0.

Uses actual t0 faces and velocities, but a manufactured uniform density to
isolate the kinematic identity. Not an evolved candidate or actual mass reset.
"""
import json
from pathlib import Path
import numpy as np

root=Path('/gpfs/kjhan/LagForce/Codex_pair_hg_modes_t0_20261018')
a=np.loadtxt(root/'raw_faces.dat.rhs.rank0000');a=a[np.argsort(a[:,0])]
f=np.loadtxt(root/'raw_faces.dat.common.rank0000');i,j=f[:,:2].astype(int).T
x=a[:,1:3]-3;r=np.linalg.norm(x,axis=1);active=r<1.9
d=np.linalg.norm(f[:,2:4],axis=1);coef=f[:,4]/(2*d)
# Bounded smooth relative rates, with several independent amplitudes.
z=np.sin(7*x[:,0])*np.cos(5*x[:,1])*np.exp(-r*r)
base=a[:,12]/a[:,6]
out={'scope':__doc__,'identity':'Vdot=Bv+L Wdot; mdot=rho L Wdot implies rhodot/rho=-Bv/V at uniform density.',
     'expected_circular_density_rate_RMS_r_lt_0p5':float(np.sqrt(np.average(base[r<.5]**2,weights=a[r<.5,7]))),
     'cases':[]}
for amplitude in [0,.001,.01,.1]:
    wdot=amplitude*z
    vw=np.bincount(i,weights=coef*(wdot[i]-wdot[j]),minlength=len(a))
    # rho=1 => m=V for this manufactured tangent state only.
    density_rate=(vw-(a[:,12]+vw))/a[:,6]
    err=float(abs(density_rate[active]+base[active]).max())
    assert err<1e-12,err
    out['cases'].append({'weight_rate_amplitude':amplitude,
       'weight_volume_rate_RMS':float(np.sqrt(np.mean((vw[active]/a[active,6])**2))),
       'density_rate_identity_max_error':err})
with (root/'weight_only_density_response.json').open('x') as fp:json.dump(out,fp,indent=2);fp.write('\n')
print(json.dumps(out,indent=2))
