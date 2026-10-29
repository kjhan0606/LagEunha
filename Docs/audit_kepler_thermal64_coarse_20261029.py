"""Compare uncapped and capped thermal64 evolution from identical t10 state."""
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_native_radial_20261026 import GAMMA_MINUS_ONE, mean
from audit_kepler_q1_thermal64_promotion_20261029 import payload, TAG

root=Path('/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64Coarse_H200_20261029')
half=Path('/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64Half_H200_20261029')
contract=json.loads((root/'contract.json').read_text());parent=root/'coarse/khout.020790.dat'
init=root/'init/khdiag.step020790.thermal64_initial.dat';ic=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/ic_v1/khout.000001.dat')
out=root/'coarse_half_comparison.json'
if out.exists():raise FileExistsError(out)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
if sha(parent)!=contract['source_sha256'] or (root/'coarse/exit_status.txt').read_text().strip()!='0':raise ValueError('Invalid coarse run')
rebased=json.loads((root/'init/initial_rebase.json').read_text())
if not rebased['passed'] or sha(init)!=rebased['initial_h200_sha256']:raise ValueError('Initial state not audited')
p0,_=ordered_global(parent);ids=p0['u4if']&ID_MASK;pi,_=ordered_global(ic,ids);p_init,m_init=ordered_global(init,ids);q0=payload(init,len(ids))
if m_init['t']!=10.:raise ValueError('Wrong initial time')
r0=np.hypot(pi['x']-3.,pi['y']-3.);m=p0['mass'].astype(float)
regions={'core':r0<.1,'adjacent':(r0>=.1)&(r0<.5),'disk':(r0>=.5)&(r0<2.)}
def energy(p,u):
    x=p['x']-3.;y=p['y']-3.;return float(np.sum(u+.5*m*(p['vx']**2+p['vy']**2)-m/np.sqrt(x*x+y*y+.01175**2),dtype=np.longdouble))
e0=energy(p0,p0['ie'].astype(float));cases={};states={}
for name,d in [('coarse',root/'coarse'),('half',half)]:
    files=[p for p in d.glob('khout.*.dat') if p.name!='khout.020790.dat']
    if len(files)!=1:raise ValueError('Ambiguous endpoint')
    p,meta=ordered_global(files[0],ids);q=payload(files[0],len(ids))
    if meta['t']!=10.1 or not np.array_equal(p0['mass'],p['mass']) or np.any(p['w2']!=0) or np.any(q['tag']!=TAG) or not np.all(np.isfinite(q['u'])) or not np.all(np.isfinite(q['v'])) or np.any(q['u']<=0) or np.any(q['v']<=0):raise ValueError('Invalid candidate state')
    log=(d/'log').read_text();accept=re.findall(r'\[WALL_ACCEPT\][^\n\[]*',log)
    ledger=re.findall(r'\[RK4E\][^\n\[]*',log)
    if len(accept)!=len(ledger) or not len(accept) or any(float(re.search(r'\b'+k+r'=(\S+)',l)[1])!=0 for l in ledger for k in ('floor_cum','de_cum','sfl_cum')):raise ValueError('Invalid energy ledgers')
    rows=[(float(re.search(r't_start=(\S+)',l)[1]),float(re.search(r'\bh=(\S+)',l)[1])) for l in accept]
    if rows[0][0]!=10. or abs(sum(rows[-1])-10.1)>1e-12 or any(abs(t+h-rows[k+1][0])>1e-12 for k,(t,h) in enumerate(rows[:-1])):raise ValueError('Noncontinuous step times')
    du=np.log(q['u']/q0['u'])+GAMMA_MINUS_ONE*np.log(q['v']/q0['v'])
    reg={}
    for label,mask in regions.items():
        w=m[mask];x=p['x'][mask]-3.;y=p['y'][mask]-3.;xb=p0['x'][mask]-3.;yb=p0['y'][mask]-3.
        L0=np.sum(w*(xb*p0['vy'][mask]-yb*p0['vx'][mask]),dtype=np.longdouble)
        L1=np.sum(w*(x*p['vy'][mask]-y*p['vx'][mask]),dtype=np.longdouble)
        reg[label]={'entropy_abs_mass_mean':mean(abs(du[mask]),w),'entropy_signed_mass_mean':mean(du[mask],w),
                    'density_abs_fraction_change':mean(abs(p['den'][mask].astype(float)/p0['den'][mask].astype(float)-1),w),
                    'angular_momentum_relative_change':float(L1/L0-1) if abs(L0)>0 else None}
    de=(energy(p,q['u'])-e0)/abs(e0)
    if abs(de)>=1e-6:raise ValueError('Authoritative energy gate failed')
    cases[name]={'steps':len(rows),'median_h':float(np.median([h for _,h in rows])),
                 'restart_energy_relative_change':de,'regions':reg,'endpoint_sha256':sha(files[0])}
    states[name]=p
delta={}
for label,mask in regions.items():
    a=states['coarse'];b=states['half'];w=m[mask]
    delta[label]={'position_difference_abs_mass_mean':mean(np.hypot(a['x'][mask]-b['x'][mask],a['y'][mask]-b['y'][mask]),w),
                  'velocity_difference_abs_mass_mean':mean(np.hypot(a['vx'][mask]-b['vx'][mask],a['vy'][mask]-b['vy'][mask]),w)}
result={'passed':True,'cases':cases,'endpoint_differences':delta,'initial_h200_sha256':sha(init),
        'scope':'Same-parent t10→10.1 controls using H200 zero-time promotion. Coarse/half actual ratio is not a uniform refinement; not a t0 accuracy or long-time proof.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
