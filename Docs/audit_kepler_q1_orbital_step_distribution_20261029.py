"""Quantify local orbital cap distribution; not an individual-step proposal."""
import hashlib
import json
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK

root=Path('/gpfs/kjhan/LagForce/Kepler_Q1_OrbitalCap006_Ada_20261029')
state=root/'to0p4125/khout.005400.dat';initial=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/ic_v1/khout.000001.dat')
contract=json.loads((root/'chain.json').read_text());output=root/'to0p4125/local_orbital_caps.json'
if output.exists():raise FileExistsError(output)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
p0,_=ordered_global(initial);p,meta=ordered_global(state,p0['u4if']&ID_MASK)
if meta['t']!=.4125 or not np.array_equal(p0['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('Invalid Q1 state')
eta=contract['eta'];eps=.01175
r=np.hypot(p['x']-3.,p['y']-3.);r0=np.hypot(p0['x']-3.,p0['y']-3.)
local=eta*(r*r+eps*eps)**.75;global_cap=float(np.min(local))
regions={}
for name,mask in {'core':r0<.1,'adjacent':(r0>=.1)&(r0<.5),'disk':(r0>=.5)&(r0<2.),'exterior':r0>=2.}.items():
    a=local[mask]
    regions[name]={'cells':int(mask.sum()),'orbital_cap_min':float(a.min()),'orbital_cap_median':float(np.median(a)),
                   'orbital_cap_p90':float(np.percentile(a,90)),'median_to_global_ratio':float(np.median(a)/global_cap)}
relative_work=float(np.sum(global_cap/local,dtype=np.longdouble)/len(local))
result={'time':meta['t'],'eta':eta,'global_orbital_cap':global_cap,'cells':len(local),'regions':regions,
        'idealized_active_site_fraction':relative_work,'idealized_inverse_fraction':1/relative_work,
        'endpoint_sha256':sha(state),
        'scope':'Instantaneous orbital-only cap distribution at one state. Idealized active-site fraction assumes independent per-cell updates and ignores hydrodynamic/face coupling, conservation, geometry rebuilds, synchronization and other timestep constraints. It is not a measured speedup, a validated individual-timestep scheme or evidence of slow pressure work.'}
with output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
