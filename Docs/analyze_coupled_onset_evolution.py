import json,re,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
root=Path('/gpfs/kjhan/LagForce/Codex_coupled_onset_seed_20261018');eps=.001;gam=float(np.float32(1.666667));initial={};final={};times={}
for name in ['reference','plus','minus']:
 d=root/name;assert (d/'evolution_exit.code').read_text().strip()=='0'
 initial[name],t0=read(d/'khout.000184.dat');final[name],times[name]=read(d/'khout.000224.dat')
 states=[dict(re.findall(r'(\w+)=([^\s]+)',l)) for l in (d/'evolution.log').read_text().splitlines() if l.startswith('[ACCEPTED_STATE]')]
 assert len(states)==40 and all(s['finite']=='1' and s['retries']=='0' for s in states)
 assert np.array_equal(initial[name]['id'],final[name]['id'])
 assert np.array_equal(initial[name]['mass'],final[name]['mass'])
assert max(times.values())==min(times.values())
def xy(a):return np.column_stack((a['x'],a['y']))
def vel(a):return np.column_stack((a['vx'],a['vy']))
def departure(a):
 x=xy(a)-3;om=np.sqrt(gravity_coefficient(np.linalg.norm(x,axis=1)))
 return vel(a)-om[:,None]*np.column_stack((-x[:,1],x[:,0]))
def metrics(states):
 a=states['reference'];p=states['plus'];n=states['minus'];m=a['mass']
 assert np.array_equal(p['id'],a['id']) and np.array_equal(n['id'],a['id'])
 dx=(xy(p)-xy(n))/(2*eps);dv=(vel(p)-vel(n))/(2*eps);du=(departure(p)-departure(n))/(2*eps);dp=(p['pressure']-n['pressure'])/(2*eps)
 omega2=gravity_coefficient(np.linalg.norm(xy(a)-3,axis=1))
 kin=float(np.sum(m[:,None]*dv**2));relative=float(np.sum(m[:,None]*du**2))
 ac=float(np.sum(a['volume']*dp**2/(gam*a['pressure'])))
 position=float(np.sum(m[:,None]*omega2[:,None]*dx**2))
 even=.5*(vel(p)+vel(n))-vel(a)
 return {'velocity':kin,'circular_departure_velocity':relative,'acoustic':ac,'omega_weighted_displacement':position,
 'velocity_plus_acoustic':kin+ac,'departure_plus_acoustic_plus_displacement':relative+ac+position,
 'even_velocity_relative_to_odd':float(np.sqrt(np.sum(m[:,None]*even**2)/kin)/eps),
 'reference_global_departure_energy':float(.5*np.sum(m[:,None]*departure(a)**2))}
start=metrics(initial);end=metrics(final)
out={'start_time':t0,'end_time':times['reference'],'epsilon':eps,'initial_metrics':start,'final_metrics':end,
 'ratios':{k:end[k]/start[k] for k in start if k!='even_velocity_relative_to_odd'},
 'scope':'40-step MPI1 OMP4 central coupled position/velocity response, isentropic initial geometry. Positive norms are diagnostic and not conserved perturbation energies; displacement weighting uses local Omega^2. No initial entropy perturbation or half-amplitude/timestep convergence. Historical long-time failure unresolved.'}
with (root/'evolution_analysis.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
