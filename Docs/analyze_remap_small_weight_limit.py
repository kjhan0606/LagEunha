"""Finite displacement consistency, not time convergence or zero-weight solver reduction."""
from pathlib import Path
import json,numpy as np
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_remap_small_weight_limit_20261018'
assert (root/'result.json').exists()
ref=np.load(root/'reference_q8.npz')['state'];m0=ref[:,0];v0=ref[:,6:8]/m0[:,None]
U0=(ref[:,15]-.5*(ref[:,12]+ref[:,14])).sum()
rows=[]
for f in (1.,.25,.0625,.015625,.00390625):
    folder=base/'Codex_checkpoint_weight_remap_20261018' if f>=.25 else root
    q=np.load(folder/f'state_w{f:g}_q8.npz')['state'];m=q[:,0];v=q[:,6:8]/m[:,None]
    internal=.5*(q[:,12]+q[:,14])-.5*m*np.sum(v*v,axis=1)
    row={'fraction':f,'relative_mass_change_L1':float(abs(m-m0).sum()/m0.sum()),'mean_velocity_change_mass_RMS':float(np.sqrt(np.average(np.sum((v-v0)**2,axis=1),weights=m0))),'internal_K_over_original_U':float(internal.sum()/U0),'minimum_U':float(np.min(q[:,15]-.5*(q[:,12]+q[:,14])))}
    rows.append(row)
for prev,row in zip(rows,rows[1:]):
    row['displacement_scaling_exponents']={key:float(np.log(prev[key]/row[key])/np.log(prev['fraction']/row['fraction'])) for key in ('relative_mass_change_L1','mean_velocity_change_mass_RMS','internal_K_over_original_U')}
out={'scope':__doc__,'rows':rows,'limits':['Conditional lift, fixed spatial resolution, single remap.','Displacement scaling does not prove repeated-remap accumulated error vanishes as timestep decreases.','Stationary-face trace mismatch is not by itself rejection evidence for all affine schemes.']}
with (root/'scaling.json').open('x') as fp:json.dump(out,fp,indent=2);fp.write('\n')
print(json.dumps(out,indent=2))
