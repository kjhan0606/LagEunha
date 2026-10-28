import sys,json,re
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_pressure_followup_20261018';old=b/'Codex_pressure_seed_20261018';g=float(np.float32(1.666667));eps=1e-5
v=lambda a:np.column_stack((a['vx'],a['vy']))
start,t0=read(b/'Codex_entropy_multistep_20261018/off/khout.002305.dat');p0,_=read(old/'plus/khout.002305.dat');m0,_=read(old/'minus/khout.002305.dat');dp0=(p0['pressure']-m0['pressure'])/(2*eps);n0=np.sum(start['volume']*dp0**2/(g*start['pressure']))
for case,meta in json.loads((root/'manifest.json').read_text()).items():
 d=root/case;assert (d/'exit.code').read_text().strip()=='0'
 states=[dict(re.findall(r'(\w+)=([^\s]+)',l)) for l in (d/'log').read_text().splitlines() if l.startswith('[ACCEPTED_STATE]')]
 assert len(states)==(32 if case.startswith('resume') else 16)
 assert all(s['finite']=='1' and s['retries']=='0' for s in states)
def response(pr,mr,rr):
 p,t=read(pr);m,tm=read(mr);ref,tr=read(rr);assert abs(t-tm)<1e-12 and abs(t-tr)<1e-12
 dv=(v(p)-v(m))/(2*eps);dp=(p['pressure']-m['pressure'])/(2*eps)
 kin=np.sum(start['mass'][:,None]*dv**2)/n0;ac=np.sum(ref['volume']*dp**2/(g*ref['pressure']))/n0
 return {'time':t,'duration':t-t0,'kinetic_ratio':float(kin),'acoustic_ratio':float(ac),'combined_ratio':float(kin+ac)},dv,dp,ref
coarse=response(old/'plus/khout.002313.dat',old/'minus/khout.002313.dat',b/'Codex_entropy_multistep_20261018/off/khout.002313.dat')
fine=response(root/'fine_plus/khout.002321.dat',root/'fine_minus/khout.002321.dat',root/'fine_reference/khout.002321.dat')
long=response(root/'resume_plus/khout.002345.dat',root/'resume_minus/khout.002345.dat',b/'Codex_pair_tangent_resume_20261018/reference/khout.002345.dat')
dv=coarse[1]-fine[1];dp=coarse[2]-fine[2];ref=fine[3]
out={'scope':'Actual full solver, distributed energy-neutral pressure seed. Positive quadratic diagnostic is not total perturbation energy. Long branches resumed; fine dt independently evolved from original perturbed checkpoint.',
 'short_coarse':coarse[0],'short_fine':fine[0],'long_coarse':long[0],
 'short_dt_response_difference_normalized_to_initial':float(np.sqrt((np.sum(start['mass'][:,None]*dv**2)+np.sum(ref['volume']*dp**2/(g*ref['pressure'])))/n0))}
with (root/'analysis.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
