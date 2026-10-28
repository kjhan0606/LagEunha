import json,re,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_onset_velocity_tangent_20261018';v=lambda a:np.column_stack((a['vx'],a['vy']))
initial={};final={};times={}
for name in ['reference','plus','minus']:
 d=root/name;assert (d/'exit.code').read_text().strip()=='0'
 a,t0=read(d/'khout.000184.dat');initial[name]=a
 outputs=[p for p in d.iterdir() if p.name.startswith('khout.') and p.name!='khout.000184.dat'];assert len(outputs)==1
 a,t=read(outputs[0]);final[name]=a;times[name]=t
 states=[dict(re.findall(r'(\w+)=([^\s]+)',l)) for l in (d/'log').read_text().splitlines() if l.startswith('[ACCEPTED_STATE]')]
 assert len(states)==40 and all(s['finite']=='1' and s['retries']=='0' for s in states)
assert max(times.values())==min(times.values())
a=initial['reference'];ref=final['reference'];m=a['mass'];eps=.001;g=float(np.float32(1.666667))
for name in initial:
 assert np.array_equal(initial[name]['id'],a['id'])
 assert np.array_equal(final[name]['id'],a['id'])
 assert np.array_equal(final[name]['mass'],m)
 for field in ['x','y','pressure','volume','ie']:
  assert np.array_equal(initial[name][field],a[field])
dv0=(v(initial['plus'])-v(initial['minus']))/(2*eps);dv=(v(final['plus'])-v(final['minus']))/(2*eps);dp=(final['plus']['pressure']-final['minus']['pressure'])/(2*eps)
n0=np.sum(m[:,None]*dv0**2);kin=np.sum(m[:,None]*dv**2);ac=np.sum(ref['volume']*dp**2/(g*ref['pressure']))
def departure_velocity(a):
 x=np.column_stack((a['x']-3,a['y']-3));r=np.linalg.norm(x,axis=1);vc=np.sqrt(gravity_coefficient(r))[:,None]*np.column_stack((-x[:,1],x[:,0]))
 return v(a)-vc
def departure(a):
 return .5*a['mass']*np.sum(departure_velocity(a)**2,axis=1)
du=(departure_velocity(final['plus'])-departure_velocity(final['minus']))/(2*eps)
departure_norm=np.sum(m[:,None]*du**2)
ic,_=read(b/'Codex_canonical_initial_20261017/run/khout.000000.dat');mask=(np.hypot(ic['x']-3,ic['y']-3)<.5)&(ic['den']<.5)
even=.5*(v(final['plus'])+v(final['minus']))-v(ref)
out={'start_time':t0,'end_time':times['reference'],'duration':times['reference']-t0,
 'epsilon':eps,'accepted_steps_per_case':40,'dt':.0025,
 'kinetic_ratio':float(kin/n0),'acoustic_ratio':float(ac/n0),
 'initial_velocity_norm_squared':float(n0),
 'circular_departure_tangent_norm_gain':float(np.sqrt(departure_norm/n0)),
 'velocity_tangent_norm_gain':float(np.sqrt(kin/n0)),'kinetic_plus_acoustic_ratio':float((kin+ac)/n0),
 'even_velocity_relative_to_odd':float(np.sqrt(np.sum(m[:,None]*even**2)/kin)/eps),
 'reference_inner_departure_energy_start':float(departure(a)[mask].sum()),'reference_inner_departure_energy_end':float(departure(ref)[mask].sum()),
 'scope':'Actual t2 early-onset velocity-pattern perturbation, 40 steps; no half-amplitude or dt convergence yet. Growing baseline departure and growing trajectory tangent are different quantities. Velocity plus acoustic norm is a positive diagnostic, not conserved perturbation energy; it omits displacement and entropy degrees of freedom.'}
out['reference_inner_departure_energy_ratio']=out['reference_inner_departure_energy_end']/out['reference_inner_departure_energy_start']
with (root/'analysis.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
