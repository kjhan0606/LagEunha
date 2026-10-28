import sys,json,re
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_pressure_seed_20261018';meta=json.loads((root/'manifest.json').read_text())
src=b/'Codex_entropy_multistep_20261018/off';ref,t=read(src/'khout.002313.dat');start,t0=read(src/'khout.002305.dat')
states={};initial={};g=float(np.float32(1.666667));vel=lambda a:np.column_stack((a['vx'],a['vy']))
for case in meta['cases']:
 d=root/case;assert (d/'exit.code').read_text().strip()=='0'
 a,ta=read(d/'khout.002313.dat');assert ta==t;states[case]=a
 a,ta=read(d/'khout.002305.dat');assert ta==t0;initial[case]=a
 assert np.array_equal(vel(a),vel(start)) and np.array_equal(a['mass'],start['mass'])
 stateslog=[dict(re.findall(r'(\w+)=([^\s]+)',l)) for l in (d/'log').read_text().splitlines() if l.startswith('[ACCEPTED_STATE]')]
 assert len(stateslog)==8 and all(a['finite']=='1' and a['retries']=='0' for a in stateslog)
out={'start_time':t0,'end_time':t,'scope':'Eight steps, distributed energy-neutral pressure/entropy seed. Positive kinetic/acoustic quadratic diagnostic only; no full stability claim.','cases':{}}
responses=[]
for name,p,m,e in [('full','plus','minus',1e-5),('half','half_plus','half_minus',5e-6)]:
 dp0=(initial[p]['pressure']-initial[m]['pressure'])/(2*e);dp=(states[p]['pressure']-states[m]['pressure'])/(2*e);dv=(vel(states[p])-vel(states[m]))/(2*e)
 norm0=np.sum(start['volume']*dp0**2/(g*start['pressure']));kin=np.sum(start['mass'][:,None]*dv**2);ac=np.sum(ref['volume']*dp**2/(g*ref['pressure']))
 even=.5*(states[p]['pressure']+states[m]['pressure'])-ref['pressure']
 out['cases'][name]={'kinetic_ratio':float(kin/norm0),'acoustic_ratio':float(ac/norm0),'combined_ratio':float((kin+ac)/norm0),
 'even_pressure_deviation_relative_to_odd':float(np.linalg.norm(even)/(e*np.linalg.norm(dp)))}
 responses.append((dv,dp,norm0))
dv=responses[0][0]-responses[1][0];dp=responses[0][1]-responses[1][1]
out['half_full_difference_normalized_to_initial_norm']=float(np.sqrt((np.sum(start['mass'][:,None]*dv**2)+np.sum(ref['volume']*dp**2/(g*ref['pressure'])))/responses[1][2]))
out['initial_delta_total_U']={k:v['delta_total_U'] for k,v in meta['cases'].items()}
with (root/'analysis.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
