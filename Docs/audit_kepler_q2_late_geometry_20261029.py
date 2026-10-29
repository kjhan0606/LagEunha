"""Decompose saved Q2 density against circular-orbit Voronoi geometry."""
import hashlib
import json
from pathlib import Path
import socket
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_q1_checkpoint_geometry_20261029 import independent
from audit_kepler_native_radial_20261026 import mean

ROOT=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
prior=json.loads((ROOT/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
out=ROOT/'KEPLER_Q2_BOUNDED_LATE_GEOMETRY_2026-10-29.json'
if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('Run geometry on compute host')
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
initial=Path(prior['initial'])
if sha(initial)!=prior['initial_sha256']:raise ValueError('Initial state changed')
a,meta=ordered_global(initial);ids=a['u4if']&ID_MASK
x=a['x']-3.;y=a['y']-3.;r=np.hypot(x,y);mass=a['mass'].astype(float)
regions={'background_0p1_to_0p5':(r>=.1)&(r<.5),'inner_disk_0p5_to_0p75':(r>=.5)&(r<.75),'disk_0p5_to_2':(r>=.5)&(r<2.)}
cases=[]
for old in prior['cases']:
    path=Path(old['endpoint'])
    if sha(path)!=old['endpoint_sha256']:raise ValueError('Endpoint changed')
    p,t=ordered_global(path,ids)
    if meta['t']!=0 or t['t']!=old['time'] or not np.array_equal(a['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('Invalid identity, time, or invariant')
    phase=t['t']*(r*r+.01175**2)**(-.75)
    xy=np.column_stack((p['x'],p['y']))
    circular=np.column_stack((3.+x*np.cos(phase)-y*np.sin(phase),3.+x*np.sin(phase)+y*np.cos(phase)))
    area,actual_check=independent(xy);circ,circ_check=independent(circular)
    vol=p['volume'].astype(float);v0=a['volume'].astype(float)
    maxerror=float(np.max(abs(vol/area-1)))
    if maxerror>=1e-6:raise ValueError('Native volume does not match independent geometry')
    metrics={'native_density_change':p['den'].astype(float)/a['den'].astype(float)-1,
             'circular_geometry_density_change':v0/circ-1,
             'native_density_relative_to_circular_geometry':circ/vol-1}
    summary={}
    for name,mask in regions.items():
        summary[name]={'cells':int(mask.sum()),**{k:{'signed_mass_mean':mean(v[mask],mass[mask]),
             'absolute_mass_mean':mean(abs(v[mask]),mass[mask]),'absolute_cell_p99':float(np.percentile(abs(v[mask]),99))} for k,v in metrics.items()}}
    cases.append({'time':t['t'],'endpoint_sha256':old['endpoint_sha256'],'native_volume_max_relative_error':maxerror,
                  'actual_geometry':actual_check,'circular_geometry':circ_check,'regions':summary})
result={'initial_sha256':prior['initial_sha256'],'cases':cases,'script_sha256':sha(Path(__file__)),
        'scope':'All-site bounded Voronoi areas, with initial sites rotated on analytic softened circular paths. The circular configuration is a kinematic counterfactual, not an evolved fluid solution; density ratios are not additive and do not isolate causes such as boundary, shear, or thermal error.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
