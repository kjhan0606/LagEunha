"""Log-pressure decomposition; entropy proxy is not a direct AV heat budget."""
import json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
base=Path('/gpfs/kjhan/LagForce')
ic,_=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat')
r=np.hypot(ic['x']-3,ic['y']-3)
masks={'inner_background':(r<.5)&(ic['den']<.5),
 'inner_disk_edge':(r<.6)&(ic['den']>.5),
 'disk_core':(r>=.6)&(r<1.9)&(ic['den']>.5)}
cases=json.loads((base/'Codex_t0_maxdt0025_t20_20261017/t12_comparison.json').read_text())['cases']
gamma=float(np.float32(1.666667))
out={'gamma':gamma,'identity':'ln(P/P0)=gamma*ln(rho/rho0)+ln[(P/rho^gamma)/(P0/rho0^gamma)]',
 'scope':'Initial-ID region masks; entropy proxy includes all numerical/physical entropy changes, not an AV-specific budget. RMS components are correlated and cannot be added.', 'cases':{}}
for name,info in cases.items():
 a,t=read(Path(info['snapshot'])); assert np.all(a['pressure']>0) and np.all(a['den']>0)
 lp=np.log(a['pressure']/ic['pressure']); compression=gamma*np.log(a['den']/ic['den']); entropy=lp-compression
 result={}
 for region,mask in masks.items():
  w=ic['mass'][mask];avg=lambda v:float(np.average(v[mask],weights=w))
  result[region]={key:{'mass_mean':avg(v),'mass_RMS':avg(v*v)**.5} for key,v in [('log_pressure',lp),('adiabatic_density_term',compression),('log_entropy_proxy',entropy)]}
  result[region]['cross_second_moment']=avg(compression*entropy)
 out['cases'][name]={'time':t,'regions':result}
p=base/'Codex_t0_maxdt0025_t20_20261017/t12_pressure_entropy.json'
with p.open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
