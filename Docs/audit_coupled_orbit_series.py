"""Gravity-only tangent reference over two time units, at two resolutions."""
import json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
b=Path('/gpfs/kjhan/LagForce');src=b/'Codex_coupled_onset_seed_20261018';root=b/'Codex_coupled_orbit_series_20261018';root.mkdir(exist_ok=False)
states=[]
for name in ['reference','plus','minus']:
 a,t0=read(src/name/'khout.000184.dat');states.append(a)
x=np.array([np.column_stack((a['x']-3,a['y']-3)) for a in states]);v=np.array([np.column_stack((a['vx'],a['vy'])) for a in states])
sel=np.any(x[1]!=x[2],axis=1)|np.any(v[1]!=v[2],axis=1)
# Unselected particles have exactly identical initial x/v in all branches;
# independent gravity gives them exactly zero central response at all times.
assert np.array_equal(x[0,~sel],x[1,~sel]) and np.array_equal(v[0,~sel],v[1,~sel])
x0=x[:,sel];v0=v[:,sel];m=states[0]['mass'][sel];eps=.001
assert np.max(np.linalg.norm(x0,axis=2))<.7
def acc(x):return -gravity_coefficient(np.linalg.norm(x,axis=2))[...,None]*x
def metric(x,v):
 om2=gravity_coefficient(np.linalg.norm(x,axis=2));vc=np.sqrt(om2)[...,None]*np.stack((-x[:,:,1],x[:,:,0]),axis=2)
 dx=(x[1]-x[2])/(2*eps);dv=(v[1]-v[2])/(2*eps);du=((v[1]-vc[1])-(v[2]-vc[2]))/(2*eps)
 return {'velocity':float(np.sum(m[:,None]*dv**2)),'circular_departure_velocity':float(np.sum(m[:,None]*du**2)),
 'omega_weighted_displacement':float(np.sum(m[:,None]*om2[0,:,None]*dx**2))}
initial=metric(x0,v0);results={}
for per_sample in [200,400]:
 dt=.1/per_sample;x=x0.copy();v=v0.copy();frames=[]
 for frame in range(20):
  for step in range(per_sample):
   a=v;b1=acc(x);c=v+dt*b1/2;d=acc(x+dt*a/2);e=v+dt*d/2;f=acc(x+dt*c/2);g=v+dt*f;h=acc(x+dt*e)
   x=x+dt*(a+2*c+2*e+g)/6;v=v+dt*(b1+2*d+2*f+h)/6
  met=metric(x,v);frames.append({'duration':.1*(frame+1),'time':t0+.1*(frame+1),'ratios':{k:met[k]/initial[k] for k in met}})
 results[str(per_sample)]={'dt':dt,'frames':frames}
out={'initial_time':t0,'selected_particles':int(sel.sum()),'initial_metrics':initial,'resolutions':results,
 'maximum_absolute_ratio_difference':{k:max(abs(a['ratios'][k]-c['ratios'][k]) for a,c in zip(results['200']['frames'],results['400']['frames'])) for k in initial},
 'scope':'Independent gravity only; restricted to initially nonzero x/v seed support because all other particles have identically zero response. No fluid feedback, acoustic metric or physical stability pass. Compare matching durations with hydro checkpoints.'}
with (root/'analysis.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps({'selected_particles':out['selected_particles'],'max_ratio_difference':out['maximum_absolute_ratio_difference'],'final':results['400']['frames'][-1]},indent=2))
