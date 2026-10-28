"""Analyze complete common checkpoints; safe while evolution is running."""
import json,sys,struct
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read,sha
from audit_canonical_pair_history import gravity_coefficient
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_coupled_onset_to_t4_20261018';src=b/'Codex_coupled_onset_seed_20261018'
outdir=root/'response_frames';outdir.mkdir(exist_ok=True)
names=['reference','plus','minus'];eps=.001;gamma=float(np.float32(1.666667));initial={}
for name in names:initial[name],t0=read(src/name/'khout.000184.dat')
mass=initial['reference']['mass'];ids=initial['reference']['id']
def metrics(states):
 x=np.array([np.column_stack((states[n]['x']-3,states[n]['y']-3)) for n in names])
 v=np.array([np.column_stack((states[n]['vx'],states[n]['vy'])) for n in names])
 om2=gravity_coefficient(np.linalg.norm(x,axis=2));vc=np.sqrt(om2)[...,None]*np.stack((-x[:,:,1],x[:,:,0]),axis=2)
 dx=(x[1]-x[2])/(2*eps);dv=(v[1]-v[2])/(2*eps);du=((v[1]-vc[1])-(v[2]-vc[2]))/(2*eps)
 p=states['reference']['pressure'];V=states['reference']['volume'];dp=(states['plus']['pressure']-states['minus']['pressure'])/(2*eps)
 def norm(z):return float(np.sum(mass[:,None]*z*z))
 kin=norm(dv);relative=norm(du);ac=float(np.sum(V*dp*dp/(gamma*p)))
 pos=float(np.sum(mass[:,None]*om2[0,:,None]*dx*dx))
 return {'velocity':kin,'circular_departure_velocity':relative,'acoustic':ac,'omega_weighted_displacement':pos,
 'departure_plus_acoustic':relative+ac,'even_velocity_relative_to_odd':float(np.sqrt(norm(.5*(v[1]+v[2])-v[0])/kin)/eps),
 'reference_global_departure_energy':.5*norm(v[0]-vc[0])}
start=metrics(initial)
common=set.intersection(*[{p.name for p in (root/n).iterdir() if p.name.startswith('khout.') and p.name.endswith('.dat')} for n in names])
frames=[];pending=[]
for filename in sorted(common):
 dest=outdir/(filename+'.json')
 if dest.exists():frames.append(json.loads(dest.read_text()));continue
 paths={n:root/n/filename for n in names};complete=True
 for p in paths.values():
  with p.open('rb') as f:header=f.read(16)
  if len(header)!=16:complete=False;break
  magic,count,size,stress=struct.unpack('<iiii',header)
  assert magic==-1 and count==len(ids) and size==624
  if p.stat().st_size!=16+count*size+52:complete=False;break
 if not complete:pending.append(filename);continue
 states={};times=[]
 for n,p in paths.items():
  states[n],t=read(p);times.append(t)
  assert np.array_equal(states[n]['id'],ids) and np.array_equal(states[n]['mass'],mass)
  assert np.all(np.isfinite(states[n]['pressure'])) and np.all(states[n]['pressure']>0)
 assert max(times)-min(times)<1e-12
 end=metrics(states)
 row={'snapshot':filename,'time':times[0],'duration_from_seed':times[0]-t0,'metrics':end,
  'ratios':{k:end[k]/start[k] for k in start if k!='even_velocity_relative_to_odd'},
  'snapshot_sha256':{n:sha(p) for n,p in paths.items()},
  'scope':'Finite +/- perturbation, not a Lyapunov exponent; orbital shear control and amplitude/timestep checks required. Even/odd ratio flags nonlinear contamination but is not a convergence proof.'}
 with dest.open('x') as f:json.dump(row,f,indent=2)
 frames.append(row)
print(json.dumps({'common_complete_frames':len(frames),'incomplete_files_skipped':pending,
 'frames':[{'time':r['time'],'snapshot':r['snapshot'],'ratios':r['ratios'],'even_odd':r['metrics']['even_velocity_relative_to_odd']} for r in frames]},indent=2))
