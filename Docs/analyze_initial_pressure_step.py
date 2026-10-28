"""Compare two actual first-step secants with independent t0 pressure RHS."""
import sys,json,re
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read,sha
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_initial_pressure_step_20261018';src=base/'Codex_pair_hg_modes_t0_20261018'
s=np.loadtxt(src/'raw_faces.dat.rhs.rank0000');s=s[np.argsort(s[:,0])]
hg=np.load(src/'instantaneous_actual.npz')['heat'];a,t0=read(src/'khout.000001.dat');assert t0==0
gamma=float(np.float32(1.666667));comp=-gamma*s[:,5]*s[:,12]/s[:,6]
heat=(gamma-1)*(s[:,11]+hg)/s[:,6];full=comp+heat
r=np.hypot(a['x']-3,a['y']-3)
out={'scope':'One actual accepted split-HG/RK4/HG step at two dt; compare pressure secant to independent t0 RHS. Local consistency, not long-time stability.','cases':{}}
for name in ['dt1e4','dt5e5']:
 d=root/name;assert (d/'exit.code').read_text().strip()=='0'
 states=[dict(re.findall(r'(\w+)=([^\s]+)',line)) for line in (d/'log').read_text().splitlines() if line.startswith('[ACCEPTED_STATE]')]
 assert len(states)==1 and states[0]['retries']=='0' and states[0]['finite']=='1'
 b,t=read(d/'khout.000002.dat');assert abs(t-float(states[0]['dt']))<1e-15
 assert np.array_equal(a['mass'],b['mass'])
 sec=(b['pressure']-a['pressure'])/t
 row={'time':t,'snapshot_sha256':sha(d/'khout.000002.dat'),'regions':{}}
 for region,mask in [('inner',r<.5),('edge',(r>=.5)&(r<.6)),('core',(r>=.6)&(r<1.9))]:
  avg=lambda q:float(np.average(q[mask],weights=a['mass'][mask]))
  rms=lambda q:avg(q*q)**.5
  row['regions'][region]={'measured_relative_pressure_rate_RMS':rms(sec/a['pressure']),
   'predicted_relative_pressure_rate_RMS':rms(full/a['pressure']),
   'full_prediction_error_relative_L2':rms(sec-full)/rms(full),
   'compression_only_error_relative_L2':rms(sec-comp)/rms(full),
   'heating_only_error_relative_L2':rms(sec-heat)/rms(full)}
 out['cases'][name]=row
with (root/'analysis.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
