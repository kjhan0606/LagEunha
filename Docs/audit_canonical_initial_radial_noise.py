"""Read existing t0 geometry/RHS on lageunha; no state modification."""
import json
from pathlib import Path
import numpy as np
from audit_canonical_pair_history import gravity_coefficient

root=Path('/gpfs/kjhan/LagForce/Codex_pair_hg_modes_t0_20261018')
state=np.loadtxt(root/'raw_faces.dat.hgstate.rank0000');state=state[np.argsort(state[:,0])]
rhs=np.loadtxt(root/'raw_faces.dat.rhs.rank0000');rhs=rhs[np.argsort(rhs[:,0])]
assert np.array_equal(state[:,0],rhs[:,0])
r=np.linalg.norm(state[:,3:5]-3,axis=1);mass=state[:,5]
omega=np.sqrt(gravity_coefficient(r))
# Logarithmic symmetric difference; two step sizes check derivative resolution.
def shear(e):
    return abs((np.sqrt(gravity_coefficient(r*np.exp(e)))-np.sqrt(gravity_coefficient(r*np.exp(-e))))/(2*e))
s=shear(1e-5);s2=shear(5e-6)
assert np.max(abs(s-s2))<1e-7
vd=rhs[:,12]/rhs[:,6]
out={'scope':'Instantaneous t0 radial profile; association is not causal isolation.',
     'definitions':{'h':'sqrt(cell geometric area)','shear':'abs(r dOmega/dr)',
                    'volume_rate':'geometric Vdot/V from original C RHS'},
     'shear_derivative_max_difference':float(np.max(abs(s-s2))), 'bins':[]}
edges=[0,.1,.2,.3,.4,.5,.6,.8,1,1.5,1.9]
for lo,hi in zip(edges[:-1],edges[1:]):
    mask=(r>=lo)&(r<hi)
    if not mask.any():continue
    w=mass[mask]
    rms=lambda a:float(np.sqrt(np.average(a[mask]**2,weights=w)))
    out['bins'].append({'r_interval':[lo,hi],'cells':int(mask.sum()),
       'density_range':[float(state[mask,10].min()),float(state[mask,10].max())],
       'h_over_r_mass_RMS':rms(np.sqrt(rhs[:,6])/np.maximum(r,1e-30)),
       'omega_mass_RMS':rms(omega),'shear_mass_RMS':rms(s),
       'volume_rate_mass_RMS':rms(vd),
       'volume_rate_to_shear_RMS_ratio':rms(vd)/rms(s) if rms(s)>0 else None})
with (root/'radial_noise.json').open('x') as f:json.dump(out,f,indent=2);f.write('\n')
print(json.dumps(out,indent=2))
