"""Saved-Q2 counterfactuals for inner-disk generator displacements."""
import hashlib
import json
from pathlib import Path
import socket
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_t3_support_20261029 import geometry, orbit
from audit_kepler_native_radial_20261026 import mean

docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
prior=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
out=docs/'KEPLER_Q2_BOUNDED_ORBIT_SUPPORT_2026-10-29.json'
if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('Run geometry on compute host')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
initial=Path(prior['initial'])
if sha(initial)!=prior['initial_sha256']:raise ValueError('Initial changed')
a,m=ordered_global(initial);ids=a['u4if']&ID_MASK
x0=np.column_stack((a['x'],a['y']));r0=np.hypot(x0[:,0]-3.,x0[:,1]-3.)
sel=np.flatnonzero((r0>=.5)&(r0<.75));weights=a['mass'][sel].astype(float)
rho0=a['den'][sel].astype(float);foreign=r0<.5
rows=[]
for old in prior['cases']:
    path=Path(old['endpoint'])
    if sha(path)!=old['endpoint_sha256']:raise ValueError('Changed endpoint')
    p,mt=ordered_global(path,ids)
    if m['t']!=0 or mt['t']!=old['time'] or not np.array_equal(a['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('Invalid state')
    xy=np.column_stack((p['x'],p['y']));exact=orbit(x0,mt['t'])
    native=p['den'][sel].astype(float)/rho0-1.
    disp=np.hypot(xy[sel,0]-exact[sel,0],xy[sel,1]-exact[sel,1])/np.sqrt(a['volume'][sel].astype(float))
    variants={'actual':np.zeros(len(a),dtype=bool), 'inner_disk_only':(r0>=.5)&(r0<.75),
              'all_inside_0p75':r0<.75,'all_inside_0p85':r0<.85}
    measures={}
    for name,mask in variants.items():
        z=xy.copy();z[mask]=exact[mask]
        area,_,_,_,polygons=geometry(z,sel,foreign,False)
        if any(np.max(np.hypot(q[:,0]-3.,q[:,1]-3.))>=2.75*np.cos(np.pi/1024) for q in polygons):
            raise ValueError('Selected polygon might touch wall')
        if name=='actual' and np.max(abs(area/p['volume'][sel].astype(float)-1))>1e-6:raise ValueError('Native areas mismatch')
        change=weights/area/rho0-1.
        measures[name]={'moved_generators':int(mask.sum()),'absolute_mass_mean':mean(abs(change),weights),
                        'signed_mass_mean':mean(change,weights),'p99_abs':float(np.percentile(abs(change),99))}
    rows.append({'time':mt['t'],'cells':len(sel),'position_error_over_sqrt_initial_area':
                 {'mass_mean':mean(disp,weights),'cell_p99':float(np.percentile(disp,99)),
                  'max':float(np.max(disp))},'actual_density_change_abs':mean(abs(native),weights),
                 'variants':measures,'endpoint_sha256':old['endpoint_sha256']})
result={'initial_sha256':prior['initial_sha256'],'cases':rows,'script_sha256':sha(Path(__file__)),
        'scope':'Each counterfactual replaces specified current generator positions with exact circular-reference positions; it is not a physical evolution, a force decomposition, a flux, or a causal explanation. Density norms are not additive.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
