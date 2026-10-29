"""Measure float U representability of a frozen native K1 half-stage predictor."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_native_radial_20261026 import mean

root=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/rhs_t10_selected_v1')
cfg=json.loads((root/'probe.json').read_text());selected=json.loads((root/'audit.json').read_text())
probe=root/'khdiag.step020790.wall_rhs_probe.dat';initial=Path(cfg['initial']);seed=Path(cfg['endpoint'])
output=root/'frozen_k1_float_resolution.json'
if output.exists():raise FileExistsError(output)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
if not selected['passed'] or sha(seed)!=cfg['endpoint_sha256'] or sha(initial)!=cfg['initial_sha256']:raise ValueError('Invalid source and selected audit')
p0,_=ordered_global(initial);ids=p0['u4if']&ID_MASK
p,meta=ordered_global(probe,ids);s,_=ordered_global(seed,ids)
if meta['t']!=10. or not all(np.array_equal(s[k],p[k]) for k in ('x','y','vx','vy','mass','ie','w2')):raise ValueError('RHS probe changed accepted state')
dtype=np.dtype({'names':['id','ax','ay','die'],'formats':['<u8','<f8','<f8','<f4'],
                'offsets':[8,80,88,120],'itemsize':576})
raw=np.memmap(probe,mode='r',dtype=dtype,offset=16,shape=(len(p),));q=raw[np.argsort(raw['id']&ID_MASK)]
if not np.array_equal(q['id']&ID_MASK,ids):raise ValueError('Raw IDs mismatch')
udot=q['die'].astype(float)-p['mass']*(p['vx']*q['ax']+p['vy']*q['ay'])
if not np.isfinite(udot).all():raise ValueError('Nonfinite K1 thermal rate')
with (root/'rk4_stage_r0.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
if len(rows)!=8:raise ValueError('Expected selected native rows')
sample=[]
for row in rows:
    idx=int(np.searchsorted(ids,int(row['id'])));native=float(row['Udot'])
    diff=abs(float(udot[idx])-native)
    if diff>1e-24+abs(native)*2e-14:raise ValueError('Raw offset or native Udot mismatch')
    sample.append(diff)
r0=np.hypot(p0['x']-3.,p0['y']-3.);mass=p['mass'].astype(float);u=p['ie'].astype(float)
regions={'core':r0<.1,'adjacent':(r0>=.1)&(r0<.5),'disk':(r0>=.5)&(r0<2.)}
steps={'coarse_median':0.0003975641302531585,'half_cap':4.760027150041424e-05,'quarter_cap':2.380013575020712e-05}
data={}
for name,h in steps.items():
    delta=.5*h*udot; rounded=(u+delta).astype(np.float32).astype(float)-u
    stats={}
    for region,mask in regions.items():
        w=mass[mask];d=delta[mask];v=rounded[mask]
        stats[region]={'cells':int(mask.sum()),'no_float_change_mass_fraction':mean((v==0).astype(float),w),
                       'nonzero_delta_no_float_change_mass_fraction':mean(((v==0)&(d!=0)).astype(float),w),
                       'requested_abs_relative_mass_mean':mean(abs(d)/u[mask],w),
                       'rounding_abs_relative_mass_mean':mean(abs(v-d)/u[mask],w),
                       'signed_rounding_relative_mass_mean':mean((v-d)/u[mask],w)}
    data[name]={'h':h,'regions':stats}
result={'passed':True,'native_selected_udot_max_abs_difference':max(sample),'cases':data,'source_sha256':sha(seed),
        'probe_sha256':sha(probe),'scope':'Actual native saved K1 RHS at t10 and float rounding of a hypothetical first RK4 half-stage predictor. Frozen only: later RHS, accumulated rounding, trajectory change, entropy production and long-time behavior are not determined.'}
with output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
