"""Exact-time gravity controls for complete live hydro checkpoints."""
import json,sys,math
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_coupled_onset_to_t4_20261018';src=b/'Codex_coupled_onset_seed_20261018';outdir=root/'orbit_comparisons';outdir.mkdir(exist_ok=True)
states=[]
for name in ['reference','plus','minus']:
 a,t0=read(src/name/'khout.000184.dat');states.append(a)
x=np.array([np.column_stack((a['x']-3,a['y']-3)) for a in states]);v=np.array([np.column_stack((a['vx'],a['vy'])) for a in states])
sel=np.any(x[1]!=x[2],axis=1)|np.any(v[1]!=v[2],axis=1)
x0=x[:,sel];v0=v[:,sel];m=states[0]['mass'][sel];eps=.001
def acceleration(x):return -gravity_coefficient(np.linalg.norm(x,axis=2))[...,None]*x
def metrics(x,v):
 om2=gravity_coefficient(np.linalg.norm(x,axis=2));vc=np.sqrt(om2)[...,None]*np.stack((-x[:,:,1],x[:,:,0]),axis=2)
 dx=(x[1]-x[2])/(2*eps);dv=(v[1]-v[2])/(2*eps);du=((v[1]-vc[1])-(v[2]-vc[2]))/(2*eps)
 return {'velocity':float(np.sum(m[:,None]*dv**2)),'circular_departure_velocity':float(np.sum(m[:,None]*du**2)),
 'omega_weighted_displacement':float(np.sum(m[:,None]*om2[0,:,None]*dx**2))}
start=metrics(x0,v0);summaries=[]
for file in sorted((root/'response_frames').iterdir()):
 if not file.name.endswith('.json'):continue
 dest=outdir/file.name
 if dest.exists():summaries.append(json.loads(dest.read_text()));continue
 row=json.loads(file.read_text());duration=row['time']-t0;controls={}
 for cap in [.0005,.00025]:
  n=math.ceil(duration/cap);dt=duration/n;x=x0.copy();v=v0.copy()
  for step in range(n):
   a=v;bb=acceleration(x);c=v+dt*bb/2;d=acceleration(x+dt*a/2);e=v+dt*d/2;f=acceleration(x+dt*c/2);g=v+dt*f;h=acceleration(x+dt*e)
   x=x+dt*(a+2*c+2*e+g)/6;v=v+dt*(bb+2*d+2*f+h)/6
  met=metrics(x,v);controls[str(cap)]={k:met[k]/start[k] for k in met}
 fine=controls['0.00025'];coarse=controls['0.0005']
 result={'time':row['time'],'duration':duration,'snapshot':row['snapshot'],
  'hydro_ratios':{k:row['ratios'][k] for k in fine},'gravity_ratios':fine,
  'hydro_over_gravity_squared_norm':{k:row['ratios'][k]/fine[k] for k in fine},
  'max_gravity_ratio_change_with_half_dt':max(abs(fine[k]-coarse[k]) for k in fine),
  'hydro_even_odd_velocity':row['metrics']['even_velocity_relative_to_odd'],
  'scope':'Matching actual checkpoint duration, no interpolation. Norm quotients are comparisons, not additive causal energy fractions or eigenmode growth. Hydro amplitude/timestep convergence pending.'}
 with dest.open('x') as f:json.dump(result,f,indent=2)
 summaries.append(result)
print(json.dumps(summaries,indent=2))
