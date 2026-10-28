"""Normalize stored candidate velocity moment defects without rerunning RHS."""
import json,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');sys.path.insert(0,str(base))
from run_canonical_av_branches_20261018 import read
root=base/'Codex_contact_prior_complete_rhs_20261018'
s=np.load(root/'rhs.npz')
h=np.load(base/'Codex_continuous_initial_rhs_quadratic_split_20261018/rhs_r0.55_f8.npz')
assert np.array_equal(s['ids'],h['ids'])
raw=np.column_stack([s['mass'],h['state'][:-4].reshape(-1,20)[:,3:]])
m=s['mass'];defect=s['moment_defects'];active=~h['frozen']
ic,_=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat')
width=np.sqrt(ic['volume'][s['ids']]);speed=np.sqrt((raw[:,13]+raw[:,16])/m)
U=raw[:,17]-.5*(raw[:,13]+raw[:,16]);assert np.all(speed>0) and np.all(U>0)
fields=dict(p=defect[:,:2]/(m*speed)[:,None],
 D=defect[:,2:6]/(m*width*speed)[:,None],
 Pi=defect[:,6:10]/(m*speed**2)[:,None],
 implied_thermal_fraction=(defect[:,6]+defect[:,9])/(2*U))
metrics={}
for key,value in fields.items():
    a=value[active].reshape(int(active.sum()),-1)
    metrics[key]=dict(max_abs=float(abs(a).max()),mass_weighted_RMS=float(np.sqrt(np.average(np.sum(a*a,axis=1),weights=m[active]))))
pressure=s['pressure_raw_rates'];remaining=s['transport_kinematic_raw_rates']
actual_hydro=np.column_stack([np.zeros(len(m)),s['rates'][:-4].reshape(-1,20)[:,3:]])
# Include the measured finite mass residual when reconstructing the whole RHS.
actual_hydro[:,0]=s['raw_derivative'][:,0]
closure=float(abs(pressure+remaining-actual_hydro).max())
print(json.dumps(dict(active_cells=int(active.sum()),normalized_moment_budget_defects=metrics,
 pressure_plus_remaining_hydro_max_defect=closure,
 scope='Momentum scale=m*v_rms, D scale=m*sqrt(area)*v_rms, Pi scale=m*v_rms^2; thermal fraction uses stored positive U. Active cells only. Frozen gravity arrays are evaluated diagnostics, not applied evolution sources.'),indent=2))
