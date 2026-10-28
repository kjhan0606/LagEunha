import sys,json,re
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_pair_tangent_branches_20261018';meta=json.loads((root/'manifest.json').read_text())
ref,t=read(Path(meta['reference'])/'khout.002313.dat');start,t0=read(Path(meta['reference'])/'khout.002305.dat')
states={};initial={}
for name in meta['cases']:
 d=root/name;assert (d/'exit.code').read_text().strip()=='0'
 a,ta=read(d/'khout.002313.dat');assert ta==t;states[name]=a
 a,ta=read(d/'khout.002305.dat');assert ta==t0;initial[name]=a
 accepted=[dict(re.findall(r'(\w+)=([^\s]+)',l)) for l in (d/'log').read_text().splitlines() if l.startswith('[ACCEPTED_STATE]')]
 assert len(accepted)==8 and all(v['finite']=='1' and v['retries']=='0' for v in accepted)
m=start['mass'];pair=np.array(meta['pair']);outside=np.ones(len(m),bool);outside[pair]=False
def vel(a):return np.column_stack((a['vx'],a['vy']))
def pos(a):return np.column_stack((a['x'],a['y']))
def norm(v):return float(np.sqrt(np.sum(m[:,None]*v*v)))
out={'start_time':t0,'end_time':t,'duration':t-t0,'pair':meta['pair'],'scope':'Central finite-difference response of actual full solver to pair velocity perturbation; all particles free, eight steps. Velocity norm gain is not total perturbation energy or an asymptotic exponent.','amplitudes':{}}
derivatives=[]
for tag,pn,mn,eps in [('full','plus','minus',1e-6),('half','half_plus','half_minus',5e-7)]:
 dv=(vel(states[pn])-vel(states[mn]))/(2*eps);dx=(pos(states[pn])-pos(states[mn]))/(2*eps)
 dv0=(vel(initial[pn])-vel(initial[mn]))/(2*eps);derivatives.append(dv)
 nonlinear=.5*(vel(states[pn])+vel(states[mn]))-vel(ref)
 out['amplitudes'][tag]={'velocity_mass_norm_gain':norm(dv)/norm(dv0),'velocity_response_outside_pair_fraction_squared':float(np.sum(m[outside,None]*dv[outside]**2)/np.sum(m[:,None]*dv**2)),
 'position_response_mass_norm':norm(dx),'initial_velocity_response_mass_norm':norm(dv0),
 'even_velocity_deviation_relative_to_odd':norm(nonlinear)/(eps*norm(dv)),
 'initial_total_momentum_difference':np.sum(m[:,None]*(vel(initial[pn])-vel(start)),axis=0).tolist()}
out['half_vs_full_velocity_derivative_relative_difference']=norm(derivatives[1]-derivatives[0])/norm(derivatives[1])
with (root/'analysis.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
