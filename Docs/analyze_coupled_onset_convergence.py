import json,re,sys,argparse
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_coupled_convergence_20261018';old=b/'Codex_coupled_onset_to_t4_20261018';seed=b/'Codex_coupled_onset_seed_20261018'
parser=argparse.ArgumentParser();parser.add_argument('--amplitude-only',action='store_true')
parser.add_argument('--root',type=Path,default=root);parser.add_argument('--coarse-step',type=int,default=263)
parser.add_argument('--intermediate',action='store_true',help='Analyze a complete saved checkpoint while the runner continues.')
args=parser.parse_args();root=args.root
manifest=json.loads((root/'manifest.json').read_text());cases={};times={}
target_time=manifest.get('target_time',manifest.get('target'))
if args.intermediate:_,target_time=read(old/'reference'/f'khout.{args.coarse_step:06d}.dat')
for name,info in manifest['cases'].items():
 if args.amplitude_only and not name.startswith('amplitude_'):continue
 d=root/name
 if not args.intermediate:assert (d/'exit.code').read_text().strip()=='0'
 steps=round((target_time-info['start_time'])/info['dt']) if args.intermediate else info['expected_steps']
 a,t=read(d/f"khout.{info.get('source_step',184)+steps:06d}.dat");cases[name]=a;times[name]=t
 rows=[dict(re.findall(r'(\w+)=([^\s]+)',l)) for l in (d/'log').read_text().splitlines() if l.startswith('[ACCEPTED_STATE]')]
 if args.intermediate:
  assert len(rows)>=steps;rows=rows[:steps]
 assert len(rows)==steps and all(s['finite']=='1' and s['retries']=='0' for s in rows)
 assert abs(t-target_time)<1e-11
coarse={}
for name in ['reference','plus','minus']:
 coarse[name],tc=read(old/name/f'khout.{args.coarse_step:06d}.dat');assert abs(tc-target_time)<1e-11
ref=coarse['reference'];m=ref['mass'];g=float(np.float32(1.666667))
for a in list(cases.values())+list(coarse.values()):
 assert np.array_equal(a['id'],ref['id']) and np.array_equal(a['mass'],m)
def vectors(a):
 x=np.column_stack((a['x']-3,a['y']-3));v=np.column_stack((a['vx'],a['vy']));om=np.sqrt(gravity_coefficient(np.linalg.norm(x,axis=1)))
 return x,v,v-om[:,None]*np.column_stack((-x[:,1],x[:,0]))
def response(p,n,eps):
 pp=vectors(p);nn=vectors(n)
 out={k:(u-v)/(2*eps) for k,u,v in zip(['position','velocity','circular_departure_velocity'],pp,nn)}
 out['pressure']=(p['pressure']-n['pressure'])/(2*eps);return out
responses={'coarse':response(coarse['plus'],coarse['minus'],.001),
 'half_amplitude':response(cases['amplitude_plus'],cases['amplitude_minus'],.0005)}
if not args.amplitude_only:responses['half_dt']=response(cases['fine_plus'],cases['fine_minus'],.001)
def norm(k,v):
 if k=='pressure':return float(np.sum(ref['volume']*v*v/(g*ref['pressure'])))
 return float(np.sum(m[:,None]*v*v))
out={'target_time':target_time,'intermediate_checkpoint':args.intermediate,'field_response_comparisons':{},'cases':manifest['cases'],
 'scope':'Central response comparison at matching target time, same binary and OMP4. Coarse/reference and refined cases can have different restart histories; discrepancy may include restart effects. Half-amplitude and half-dt comparisons are separate, not a formal asymptotic order estimate.'}
for name in responses:
 if name=='coarse':continue
 out['field_response_comparisons'][name]={k:{'relative_L2_difference':float(np.sqrt(norm(k,responses[name][k]-v)/norm(k,v))),
 'squared_norm_ratio_to_coarse':norm(k,responses[name][k])/norm(k,v)} for k,v in responses['coarse'].items()}
initial={}
for name in ['plus','minus']:initial[name],_=read(seed/name/'khout.000184.dat')
r0=response(initial['plus'],initial['minus'],.001)
# Use common final weights here, explicitly distinct from evolving-weight ratios.
out['common_final_weight_growth_ratios']={name:{k:norm(k,v)/norm(k,r0[k]) for k,v in resp.items()} for name,resp in responses.items()}
filename='amplitude_analysis' if args.amplitude_only else 'analysis'
if args.intermediate:filename+=f'_checkpoint{args.coarse_step:06d}'
with (root/(filename+'.json')).open('x') as f:json.dump(out,f,indent=2)
print(json.dumps({k:v for k,v in out.items() if k!='cases'},indent=2))
