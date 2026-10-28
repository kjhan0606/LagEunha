import json,re,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_coupled_restart_control_20261018';old=b/'Codex_coupled_onset_to_t4_20261018';out={}
for name in ['plus','minus']:
 d=root/name;assert (d/'exit.code').read_text().strip()=='0'
 a,t=read(d/'khout.000263.dat');z,tz=read(old/name/'khout.000263.dat');assert abs(t-tz)<1e-12
 states=[dict(re.findall(r'(\w+)=([^\s]+)',l)) for l in (d/'log').read_text().splitlines() if l.startswith('[ACCEPTED_STATE]')]
 assert len(states)==79 and all(s['finite']=='1' and s['retries']=='0' for s in states)
 assert np.array_equal(a['id'],z['id']) and np.array_equal(a['mass'],z['mass'])
 fields={}
 for k in a.dtype.names:
  fields[k]={'bitwise_equal':bool(np.array_equal(a[k],z[k])),'max_abs_difference':float(np.max(abs(a[k].astype(float)-z[k].astype(float))))}
 out[name]={'time':t,'fields':fields,'all_reader_fields_bitwise_equal':all(v['bitwise_equal'] for v in fields.values())}
out['scope']='Same binary/OMP4/full-amplitude seeds and dt.0025; uninterrupted79 steps versus40+39 with restart. Reader fields only, not all snapshot bytes; no unperturbed reference case.'
with (root/'analysis.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
