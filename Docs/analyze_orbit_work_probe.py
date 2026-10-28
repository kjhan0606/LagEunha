import json,re,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
root=Path('/gpfs/kjhan/LagForce/Codex_orbit_work_20261018')
def energy(a):
 x=np.column_stack((a['x']-3,a['y']-3));r=np.linalg.norm(x,axis=1)
 vc=np.sqrt(gravity_coefficient(r))[:,None]*np.column_stack((-x[:,1],x[:,0]))
 return float((.5*a['mass']*((np.column_stack((a['vx'],a['vy']))-vc)**2).sum(axis=1)).sum())
states={};out={}
for name,steps in [('off',1),('on',1),('half',2)]:
 d=root/name;assert (d/'exit.code').read_text().strip()=='0'
 a,t0=read(d/'khout.000224.dat');z,t1=read(d/f'khout.{224+steps:06d}.dat');states[name]=z
 assert np.array_equal(a['id'],z['id']) and np.array_equal(a['mass'],z['mass'])
 rows={k:[] for k in ['ACCEPTED_STATE','HG_PERTURBATION','ORBIT_WORK_STAGE']}
 for l in (d/'log').read_text().splitlines():
  for k in rows:
   if l.startswith('['+k+']'):rows[k].append(dict(re.findall(r'(\w+)=([^\s]+)',l)))
 assert len(rows['ACCEPTED_STATE'])==steps
 assert all(s['finite']=='1' and s['retries']=='0' for s in rows['ACCEPTED_STATE'])
 assert len(rows['HG_PERTURBATION'])==steps*8
 assert len(rows['ORBIT_WORK_STAGE'])==(0 if name=='off' else steps*4)
 delta=energy(z)-energy(a);hg=sum(float(s['dEpert']) for s in rows['HG_PERTURBATION'])
 result={'start':t0,'end':t1,'departure_delta':delta,'HG_delta':hg,'RK_delta_by_subtraction':delta-hg}
 if name!='off':
  assert [int(s['stage']) for s in rows['ORBIT_WORK_STAGE']]==[1,2,3,4]*steps
  integrals={k:sum(float(s[k])*float(s['wdt']) for s in rows['ORBIT_WORK_STAGE']) for k in ['hydro_minus_TCF','TCF','shear','total']}
  result.update(stage_integrals=integrals,closure_residual=delta-hg-integrals['total'])
 out[name]=result
fields=['x','y','vx','vy','mass','pressure','ie','den','volume','csound']
# Report every reader field rather than assuming a naming convention.
common=states['off'].dtype.names
equal={f:bool(np.array_equal(states['off'][f],states['on'][f])) for f in common}
assert all(equal.values())
out['same_binary_on_off_fields_bitwise_equal']=equal
out['scope']='One full step versus two half steps from same preserved checkpoint; stage quadrature closure, not long physical accuracy. Hydro minus TCF includes any other active hydro forces. No retry exercise yet.'
with (root/'analysis.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
