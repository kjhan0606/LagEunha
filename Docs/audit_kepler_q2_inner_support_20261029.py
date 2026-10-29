"""Measure how inner-background positions affect saved Q2 disk Voronoi cells."""
import hashlib
import json
from pathlib import Path
import socket
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_t3_support_20261029 import geometry, orbit
from audit_kepler_native_radial_20261026 import mean

docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
previous=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
out=docs/'KEPLER_Q2_BOUNDED_INNER_SUPPORT_2026-10-29.json'
if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('Run geometry on compute host')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
initial=Path(previous['initial'])
if sha(initial)!=previous['initial_sha256']:raise ValueError('Initial changed')
a,m=ordered_global(initial);ids=a['u4if']&ID_MASK
x0=np.column_stack((a['x'],a['y']));r0=np.hypot(x0[:,0]-3.,x0[:,1]-3.)
sel=np.flatnonzero((r0>=.5)&(r0<.75));bg=r0<.5
weights=a['mass'][sel].astype(float);rho0=a['den'][sel].astype(float)
cases=[]
for old in previous['cases']:
    path=Path(old['endpoint'])
    if sha(path)!=old['endpoint_sha256']:raise ValueError('Endpoint changed')
    p,meta=ordered_global(path,ids)
    if m['t']!=0 or meta['t']!=old['time'] or not np.array_equal(a['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('Unexpected state')
    xy=np.column_stack((p['x'],p['y']))
    area,support,edges,_,polygons=geometry(xy,sel,bg,False)
    error=float(np.max(abs(area/p['volume'][sel].astype(float)-1)))
    if error>1e-6:raise ValueError('Actual selected area differs from stored bounded volume')
    corrected=xy.copy();exact=orbit(x0,meta['t']);corrected[bg]=exact[bg]
    c_area,_,_,_,c_polygons=geometry(corrected,sel,bg,False)
    apothem=2.75*np.cos(np.pi/1024)
    if any(np.max(np.hypot(q[:,0]-3.,q[:,1]-3.))>=apothem for q in polygons+c_polygons):
        raise ValueError('A selected polygon is not guaranteed inside the wall')
    actual=p['den'][sel].astype(float)/rho0-1.
    alternative=weights/c_area/rho0-1.
    norm=weights*abs(actual)
    cases.append({'time':meta['t'],'cells':len(sel),'background_support_cells':int(support.sum()),
                  'supporting_faces':len(edges),'fraction_of_inner_disk_density_error_on_supported_cells':float(np.sum(norm[support])/np.sum(norm)),
                  'actual_absolute_density_change':mean(abs(actual),weights),
                  'background_replaced_by_circular_positions_absolute_density_change':mean(abs(alternative),weights),
                  'max_native_volume_relative_error':error,'endpoint_sha256':old['endpoint_sha256']})
result={'initial_sha256':previous['initial_sha256'],'cases':cases,'script_sha256':sha(Path(__file__)),
        'scope':'Read-only Voronoi counterfactual: background generator positions replaced by their exact circular reference while all disk/exterior positions and selected masses remain fixed. Not an evolved fluid state, a proof of causality in time, or a measurement of flux, work, or entropy.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
