"""Locate response growth by fixed initial radial bins at a saved checkpoint."""
import json,sys,argparse
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
ap=argparse.ArgumentParser();ap.add_argument('--step',type=int,required=True);args=ap.parse_args()
b=Path('/gpfs/kjhan/LagForce');src=b/'Codex_coupled_onset_seed_20261018';root=b/'Codex_coupled_onset_to_t4_20261018';eps=.001;g=float(np.float32(1.666667))
ic,_=read(b/'Codex_canonical_initial_20261017/run/khout.000000.dat');r0=np.hypot(ic['x']-3,ic['y']-3)
def terms(folder,step):
 states=[]
 for name in ['reference','plus','minus']:
  a,t=read(folder/name/f'khout.{step:06d}.dat');assert np.array_equal(a['id'],ic['id']);states.append(a)
 def u(a):
  x=np.column_stack((a['x']-3,a['y']-3));v=np.column_stack((a['vx'],a['vy']));om=np.sqrt(gravity_coefficient(np.linalg.norm(x,axis=1)))
  return v-om[:,None]*np.column_stack((-x[:,1],x[:,0]))
 du=(u(states[1])-u(states[2]))/(2*eps);dp=(states[1]['pressure']-states[2]['pressure'])/(2*eps);ref=states[0]
 return {'departure_velocity':ref['mass']*np.sum(du**2,axis=1),'acoustic':ref['volume']*dp**2/(g*ref['pressure'])},t
start,t0=terms(src,184);end,t=terms(root,args.step);edges=[0,.1,.2,.3,.4,.5,.6,1.,1.9];rows=[]
for lo,hi in zip(edges[:-1],edges[1:]):
 mask=(r0>=lo)&(r0<hi);row={'initial_radius':[lo,hi],'cells':int(mask.sum()),'metrics':{}}
 for k in start:
  a=float(start[k][mask].sum());z=float(end[k][mask].sum());row['metrics'][k]={'initial':a,'final':z,'ratio':z/a if a>0 else None,'fraction_of_global_final':z/float(end[k].sum())}
 rows.append(row)
out={'time':t,'seed_time':t0,'snapshot_step':args.step,'regions':rows,
 'top10_cell_fraction':{k:float(np.sort(v)[-10:].sum()/v.sum()) for k,v in end.items()},
 'scope':'Fixed initial radial material bins, central finite perturbation squared response. Different from unperturbed trajectory energy and not causal localization. Outside1.9 omitted from bin table but included in global normalization.'}
dest=root/f'regional_response_{args.step:06d}.json'
with dest.open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
