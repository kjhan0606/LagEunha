"""Integrate accepted RK work and HG increments from the continued run."""
import json,re,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
root=Path('/gpfs/kjhan/LagForce/Codex_orbit_work_interval_20261018')
assert (root/'exit.code').read_text().strip()=='0'
def energy(a):
 x=np.column_stack((a['x']-3,a['y']-3));r=np.linalg.norm(x,axis=1)
 vc=np.sqrt(gravity_coefficient(r))[:,None]*np.column_stack((-x[:,1],x[:,0]))
 return float(np.sum(.5*a['mass']*np.sum((np.column_stack((a['vx'],a['vy']))-vc)**2,axis=1)))
a,t0=read(root/'khout.000225.dat');z,t1=read(root/'khout.000265.dat')
assert np.array_equal(a['id'],z['id']) and np.array_equal(a['mass'],z['mass'])
rows={k:[] for k in ['ACCEPTED_STATE','HG_PERTURBATION','ORBIT_WORK_STAGE']}
for line in (root/'log').read_text().splitlines():
 for k in rows:
  if line.startswith('['+k+']'):rows[k].append(dict(re.findall(r'(\w+)=([^\s]+)',line)))
assert len(rows['ACCEPTED_STATE'])==40
assert all(s['finite']=='1' and s['retries']=='0' for s in rows['ACCEPTED_STATE'])
assert len(rows['ORBIT_WORK_STAGE'])==160 and len(rows['HG_PERTURBATION'])==320
assert [int(s['stage']) for s in rows['ORBIT_WORK_STAGE']]==[1,2,3,4]*40
assert abs(t1-t0-.1)<1e-12
keys=['hydro_minus_TCF','TCF','shear','total']
increments=[]
for i in range(40):
 stages=rows['ORBIT_WORK_STAGE'][4*i:4*i+4]
 inc={k:sum(float(s[k])*float(s['wdt']) for s in stages) for k in keys}
 inc['HG']=sum(float(s['dEpert']) for s in rows['HG_PERTURBATION'][8*i:8*i+8])
 inc['step']=i+1;increments.append(inc)
totals={k:sum(s[k] for s in increments) for k in keys+['HG']}
e0=energy(a);e1=energy(z);res=e1-e0-totals['total']-totals['HG']
out={'start':t0,'end':t1,'accepted_steps':40,'initial_departure_energy':e0,'final_departure_energy':e1,
 'departure_energy_ratio':e1/e0,'departure_delta':e1-e0,'integrals':totals,'closure_residual':res,
 'closure_relative_to_net_change':res/(e1-e0),
 'step_sign_counts':{k:{'positive':sum(s[k]>0 for s in increments),'negative':sum(s[k]<0 for s in increments)} for k in keys+['HG']},
 'step_increments':increments,
 'scope':'Global accepted trajectory work, no retries. RK rates integrated with stage weights; finite-step closure retained. Hydro-minus-TCF not assumed pure pressure. No half-dt or independent long-interval passivity check. Historical t6-10 failure remains unresolved.'}
with (root/'analysis.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps({k:v for k,v in out.items() if k!='step_increments'},indent=2))
