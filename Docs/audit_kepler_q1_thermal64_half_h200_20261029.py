"""Same-parent arithmetic comparison; use authoritative double U/V for thermal64."""
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_native_radial_20261026 import GAMMA_MINUS_ONE, mean
from audit_kepler_q1_thermal64_promotion_20261029 import payload, TAG

root=Path('/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64Half_H200_20261029')
contract=json.loads((root/'contract.json').read_text())
ref=Path(contract['reference_half']); source=Path(contract['source'])
promotion=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/thermal64_init_t10_v1/khdiag.step020790.thermal64_initial.dat')
initial=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/ic_v1/khout.000001.dat')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
if sha(source)!=contract['source_sha256'] or sha(root/'eunha2.exe')!=contract['binary_sha256']:raise ValueError('Changed source or candidate')
if (root/'exit_status.txt').read_text().strip()!='0':raise ValueError('Native trial failed')
gate=json.loads((ref/'checkpoint_gate.json').read_text());done=json.loads((ref/'completed.json').read_text())
single=Path(gate['endpoint'])
if len(gate['checks'])!=28 or not all(gate['checks'].values()) or gate['time']!=10.1 or done['sha256']!=sha(single):raise ValueError('Invalid completed reference')
p0,m0=ordered_global(source);ids=p0['u4if']&ID_MASK
ic,_=ordered_global(initial,ids);ps,ms=ordered_global(single,ids)
pt0,mt0=ordered_global(promotion,ids);qt0=payload(promotion,len(ids))
if sha(promotion)!='0456a8c6dafb3f7ee7e681089c8c475462a4fe7e4ee062a19cf19f4f2db52381':raise ValueError('Unexpected promoted reference')
endpoints=[p for p in root.glob('khout.*.dat') if p.name!=source.name]
if len(endpoints)!=1:raise ValueError('Expected exactly one candidate snapshot')
end=endpoints[0];pt,mt=ordered_global(end,ids);qt=payload(end,len(ids))
if not (m0['t']==mt0['t']==10. and ms['t']==mt['t']==10.1):raise ValueError('Wrong physical times')
if not all(np.array_equal(p0['mass'],p['mass']) and np.all(p['w2']==0) for p in [ps,pt,pt0]):raise ValueError('Mass or weight changed')
if not all(np.all(q['tag']==TAG) and np.all(np.isfinite(q['u'])) and np.all(np.isfinite(q['v'])) and np.all(q['u']>0) and np.all(q['v']>0) for q in [qt0,qt]):raise ValueError('Invalid thermal64 payload')
if not (np.array_equal(qt['u'].astype(np.float32),pt['ie']) and np.array_equal(qt['v'].astype(np.float32),pt['volume'])):raise ValueError('Invalid thermal mirror')
if not all(np.isfinite(pt[k]).all() for k in ('x','y','vx','vy','mass','den','pressure','volume','ie','w2')):raise ValueError('Nonfinite endpoint state')
rho=pt['mass'].astype(float)/qt['v'];pressure=GAMMA_MINUS_ONE*qt['u']/qt['v']
if not (np.all(rho>0) and np.all(pressure>0) and np.isfinite(rho).all() and np.isfinite(pressure).all()
        and np.max(abs(pt['den']/rho-1))<1e-6 and np.max(abs(pt['pressure']/pressure-1))<1e-6):raise ValueError('Invalid double-state EOS')
domain=.5*1024*2.75**2*np.sin(2*np.pi/1024)
area_error=abs(float(np.sum(qt['v'],dtype=np.longdouble))-domain)
if area_error>=1e-5:raise ValueError('Double volumes do not tile bounded polygon')
x=pt['x']-3.;y=pt['y']-3.;theta=np.arctan2(y,x);nearest=np.rint(theta/(2*np.pi/1024)-.5)
apothem=2.75*np.cos(np.pi/1024);gap=np.full(len(pt),-np.inf)
for shift in (-1,0,1):
    angle=(nearest+shift+.5)*2*np.pi/1024
    gap=np.maximum(gap,x*np.cos(angle)+y*np.sin(angle)-apothem)
if np.any(gap>64*np.finfo(float).eps*np.maximum(2.75,abs(x)+abs(y))):raise ValueError('Endpoint outside fixed polygon')
def energy(p,u):
    mass=p['mass'].astype(float);xx=p['x'].astype(float)-3.;yy=p['y'].astype(float)-3.
    vx=p['vx'].astype(float);vy=p['vy'].astype(float)
    return float(np.sum(u+.5*mass*(vx*vx+vy*vy)-mass/np.sqrt(xx*xx+yy*yy+.01175**2),dtype=np.longdouble))
e0=energy(p0,p0['ie'].astype(float));e1=energy(pt,qt['u']);de=(e1-e0)/abs(e0)
if abs(de)>=1e-6:raise ValueError('Double authoritative restart energy gate failed')
def steps(d):
    log=(d/'log').read_text(errors='replace'); rows=[]
    for line in log.splitlines():
        if '[WALL_ACCEPT]' in line:rows.append((float(re.search(r't_start=(\S+)',line)[1]),float(re.search(r'\bh=(\S+)',line)[1])))
    if not rows or rows[0][0]!=10. or abs(sum(rows[-1])-10.1)>1e-12 or any(abs(t+h-rows[k+1][0])>1e-12 for k,(t,h) in enumerate(rows[:-1])):raise ValueError('Invalid accepted time sequence')
    if any(not (0<h<=contract['dt_cap']*(1+1e-12)) for _,h in rows):raise ValueError('Invalid time step cap')
    if '[WALL_EVENT_STOP]' in log or 'payload_counts_dtold_ledgers_exact=0' in log:raise ValueError('Wall failure')
    ledger=[l for l in log.splitlines() if '[RK4E]' in l]
    if len(ledger)!=len(rows) or not all(all(float(re.search(r'\b'+k+r'=(\S+)',l)[1])==0 for k in ('floor_cum','de_cum','sfl_cum')) for l in ledger):raise ValueError('Unexpected injection or ledger')
    return {'steps':len(rows),'median_h':float(np.median([h for _,h in rows]))}
step_single=steps(ref);step_double=steps(root)
r0=np.hypot(ic['x']-3.,ic['y']-3.)
regions={'core_r0_lt_0p1':r0<.1,'core_r0_0p1_to0p5':(r0>=.1)&(r0<.5),'disk_r0_0p5_to2':(r0>=.5)&(r0<2.)}
ds=np.log(ps['ie'].astype(float)/p0['ie'].astype(float))+GAMMA_MINUS_ONE*np.log(ps['volume'].astype(float)/p0['volume'].astype(float))
dt=np.log(qt['u']/qt0['u'])+GAMMA_MINUS_ONE*np.log(qt['v']/qt0['v'])
out={}
for name,mask in regions.items():
    mass=p0['mass'][mask].astype(float)
    def stats(a):return {'signed_mass_mean':mean(a,mass),'abs_mass_mean':mean(abs(a),mass),'abs_cell_p99':float(np.percentile(abs(a),99))}
    out[name]={'cells':int(mask.sum()),'single_incremental_log_entropy':stats(ds[mask]),'thermal64_incremental_log_entropy':stats(dt[mask]),
               'endpoint_log_entropy_difference':stats((np.log(qt['u']/ps['ie'].astype(float))+GAMMA_MINUS_ONE*np.log(qt['v']/ps['volume'].astype(float)))[mask])}
result={'passed':True,'single':step_single,'thermal64':step_double,'regions':out,'source_sha256':sha(source),'single_endpoint_sha256':sha(single),'thermal64_endpoint_sha256':sha(end),
        'promotion_sha256':sha(promotion),'double_restart_energy_relative_change':de,'double_polygon_area_absolute_error':area_error,
        'double_max_wall_gap':float(gap.max()),'double_eos_max_relative_errors':{'density':float(np.max(abs(pt['den']/rho-1))), 'pressure':float(np.max(abs(pt['pressure']/pressure-1)))},
        'scope':'Same-parent short t10 to t10.1 comparison; initial double volume uses independent Ada promotion and may have small cross-GPU roundoff. Neither inherited-error repair nor long-time stability nor convergence order is established.'}
output=root/'arithmetic_comparison.json'
if output.exists():raise FileExistsError(output)
with output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
