"""Accepted HG budget over completed onset branches; no new evolution."""
import json,re,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
root=Path('/gpfs/kjhan/LagForce/Codex_onset_velocity_tangent_20261018')
def departure(a):
 x=np.column_stack((a['x']-3,a['y']-3));r=np.linalg.norm(x,axis=1)
 vc=np.sqrt(gravity_coefficient(r))[:,None]*np.column_stack((-x[:,1],x[:,0]))
 return .5*a['mass']*((np.column_stack((a['vx'],a['vy']))-vc)**2).sum(axis=1),r
results={}
for name in ['reference','plus','minus']:
 d=root/name;assert (d/'exit.code').read_text().strip()=='0'
 a,t0=read(d/'khout.000184.dat');z,t1=read(d/'khout.000224.dat')
 assert np.array_equal(a['id'],z['id']) and np.array_equal(a['mass'],z['mass'])
 rows={key:[] for key in ['ACCEPTED_STATE','HG_PERTURBATION','HG_TRANSFER']}
 for line in (d/'log').read_text().splitlines():
  for key in rows:
   if line.startswith('['+key+']'):
    rows[key].append(dict(re.findall(r'(\w+)=([^\s]+)',line)))
 assert len(rows['ACCEPTED_STATE'])==40
 assert all(s['finite']=='1' and s['retries']=='0' for s in rows['ACCEPTED_STATE'])
 assert len(rows['HG_PERTURBATION'])==320 and len(rows['HG_TRANSFER'])==320
 assert all(s['rank']=='0' for s in rows['HG_PERTURBATION'])
 e0,r0=departure(a);e1,r1=departure(z)
 hg=sum(float(s['dEpert']) for s in rows['HG_PERTURBATION'])
 delta=float((e1-e0).sum())
 inner_hg=sum(float(s['dEpert']) for s in rows['HG_PERTURBATION'] if s['bin']=='0')
 mask=r0<.5
 results[name]={
  't0':t0,'t1':t1,'steps':40,
  'global_departure_energy_start':float(e0.sum()),
  'global_departure_energy_end':float(e1.sum()),
  'global_departure_delta':delta,'global_HG_delta':hg,
  'global_RK_delta_by_subtraction':delta-hg,
  'HG_thermal_energy_added':sum(float(s['dU']) for s in rows['HG_TRANSFER']),
  'inner_start_material_delta':float((e1-e0)[mask].sum()),
  'inner_HG_instantaneous_bin_delta':inner_hg,
  'inner_endpoint_membership_changes':int(np.count_nonzero(mask!=(r1<.5))),
  'HG_identity_abs_sum':sum(float(s['identity_abs_sum']) for s in rows['HG_PERTURBATION'])}
out={'cases':results,'scope':'Global HG exact finite increments summed over 80 split applications, all 40 steps accepted without retries. RK is total endpoint change minus HG, not independent force integration. Inner HG uses instantaneous radius while inner endpoint difference follows initial IDs; do not infer an exact regional budget without intermediate membership tracking. Pressure/TCF/shear separation not available in these logs.'}
with (root/'operator_budget.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
