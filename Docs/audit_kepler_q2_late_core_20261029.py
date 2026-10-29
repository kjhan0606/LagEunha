"""Compare saved Q2 early/late states by original radius without evolving them."""
import hashlib
import json
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_q1_restart_dt_bundle_20261029 import region_metrics
from audit_kepler_native_radial_20261026 import energy

initial=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/t001_v2/khout.000000.dat')
endpoints=[Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to30_a10_v1/khout.005247.dat'),
           Path('/gpfs/kjhan/LagForce/Kepler_Q2_WallCrossing_Replay_20261029/native_events_v2/compact_crossing/khout.011597.dat')]
out=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs/KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json')

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

a,meta0=ordered_global(initial);ids=a['u4if']&ID_MASK
if meta0['t']!=0. or len(a)!=55351:raise ValueError('Unexpected initial Q2')
r=np.hypot(a['x']-3.,a['y']-3.)
regions={'core_r_lt_0p1':r<.1,'adjacent_0p1_to_0p5':(r>=.1)&(r<.5),
         'inner_disk_0p5_to_0p75':(r>=.5)&(r<.75),'disk_0p5_to_2':(r>=.5)&(r<2.)}
cases=[]
for path,expected in zip(endpoints,[30.,48.4]):
    p,meta=ordered_global(path,ids)
    if meta['t']!=expected or not np.array_equal(a['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('Changed identity, time, mass, or weight')
    rows={name:region_metrics(a,a,p,mask,meta['t']) for name,mask in regions.items() if mask.any()}
    cases.append({'time':meta['t'],'endpoint':str(path),'endpoint_sha256':sha(path),
                  'energy_relative_to_initial':(energy(p,'baseline')-energy(a,'baseline'))/abs(energy(a,'baseline')),
                  'regions':rows})
result={'initial':str(initial),'initial_sha256':sha(initial),'cases':cases,'audit_sha256':sha(Path(__file__)),
        'scope':'Postprocess saved aligned Q2 t30 and event-continued t48.4 states using original-radius cohorts and softened circular trajectories. Time 48.4 includes a wall impact; no boundary convergence, spin, generic IC, or GIZMO long-time validation is implied. Regional energy/angular momentum are not separately conserved.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
