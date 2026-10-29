"""Compare the saved Q2 pair growth projection with local acoustic rates."""
import hashlib
import json
from pathlib import Path
import socket
import numpy as np
from scipy.spatial import cKDTree
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_t3_support_20261029 import orbit
from audit_kepler_native_radial_20261026 import mean

docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
late=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
feedback=json.loads((docs/'KEPLER_Q2_TORQUE_FEEDBACK_2026-10-29.json').read_text())
out=docs/'KEPLER_Q2_FEEDBACK_SOUND_SCALE_2026-10-29.json'
if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('Run analysis on compute host')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
initial=Path(late['initial']);a,ma=ordered_global(initial);ids=a['u4if']&ID_MASK
if sha(initial)!=late['initial_sha256'] or ma['t']!=0:raise ValueError('Changed IC')
xy0=np.column_stack((a['x'],a['y']));r0=np.hypot(xy0[:,0]-3.,xy0[:,1]-3.)
sel0=np.flatnonzero((r0>=.5)&(r0<.75))
root=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1')
path_last=root/'khout.011058.dat';p,mt=ordered_global(path_last,ids)
xy=np.column_stack((p['x'],p['y']));d,j=cKDTree(xy).query(xy[sel0],k=2,workers=1)
exact=orbit(xy0,mt['t']);ref=np.linalg.norm(exact[j[:,1]]-exact[sel0],axis=1)
sel=sel0[d[:,1]/ref<.8];mass=a['mass'][sel].astype(float)
if len(sel)!=244:raise ValueError('Selected cohort changed')
paths=[Path(late['cases'][0]['endpoint'])]+[root/('khout.'+n+'.dat') for n in ('001329','002762','004411','005767','007138','008700','011058')]
cases=[];gamma=float(np.float32(1.666667))
for path in paths:
    p,mt=ordered_global(path,ids);digest=sha(path)
    q=next((v for v in feedback['cases'] if v['time']==mt['t']),None)
    if q is None or digest!=q['source_sha256'] or not np.array_equal(a['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('Unexpected state')
    pressure=p['pressure'][sel].astype(float);density=p['den'][sel].astype(float);volume=p['volume'][sel].astype(float)
    if np.any(pressure<=0) or np.any(density<=0) or np.any(volume<=0):raise ValueError('Invalid physical state')
    sound=np.sqrt(gamma*pressure/density);length=np.sqrt(volume);rate=sound/length
    acoustic=mean(rate,mass);growth=q['split_weighted_growth_rate_per_time']
    cases.append({'time':mt['t'],'sites':len(sel),'sound_over_sqrt_area_mass_mean':acoustic,
                  'sound_over_sqrt_area_cell_median':float(np.median(rate)),
                  'pressure_torque_projection_rate':growth,'growth_over_mean_acoustic_rate':growth/acoustic,
                  'source_sha256':digest})
result={'initial_sha256':late['initial_sha256'],'cases':cases,'script_sha256':sha(Path(__file__)),
        'scope':'Fixed late-selected Q2 sites; instantaneous |pressure-torque projection| compared with cs/sqrt(area), not a mode measurement or convergence study. gamma matches native float parameter. Cohort is conditioned on later pairing; contrast with prior historical mechanism is suggestive only.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
