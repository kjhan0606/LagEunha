"""Compare late Q2 L modes to initial ring-cell-area modes; exploratory."""
import hashlib
import json
import socket
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ID_MASK, ordered_global

if socket.gethostname().split('.')[0].lower() in ('grammar', 'syntax'):
    raise RuntimeError('run on a compute host')
docs = Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
late = json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
spectrum = json.loads((docs/'KEPLER_Q2_RING_SPECTRUM_2026-10-29.json').read_text())
root = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1')

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''): h.update(chunk)
    return h.hexdigest()

initial = Path(late['initial'])
a, meta = ordered_global(initial)
if meta['t'] != 0 or sha(initial) != late['initial_sha256']: raise ValueError('initial hash')
ids = a['u4if'] & ID_MASK
x0 = a['x'].astype(float)-3; y0 = a['y'].astype(float)-3
r0 = np.hypot(x0,y0); angle = np.arctan2(y0,x0)
rounded = np.round(r0,10); selected = (r0>=.5)&(r0<.75)
all_r = np.unique(rounded[(r0>=.5)&(r0<2)])
rings = []
for radius in np.unique(rounded[selected]):
    ix = np.flatnonzero(selected & (rounded==radius))
    ix = ix[np.argsort(angle[ix])]
    gaps = np.diff(np.r_[angle[ix],angle[ix][0]+2*np.pi])
    if len(ix)<24 or np.max(abs(gaps/(2*np.pi/len(ix))-1))>1e-4: raise ValueError('nonuniform ring')
    k = int(np.searchsorted(all_r,radius))
    neighbors = [int(np.count_nonzero(rounded==all_r[j])) for j in (k-1,k+1) if 0<=j<len(all_r)]
    rings.append((ix,neighbors))
L0=x0*a['vy'].astype(float)-y0*a['vx'].astype(float)
paths=[Path(late['cases'][0]['endpoint']),root/'khout.011058.dat']
cases=[]
for path in paths:
    p, meta = ordered_global(path,ids); digest=sha(path)
    ref=next((v for v in spectrum['cases'] if v['time']==meta['t']),None)
    if ref is None or ref['source_sha256']!=digest or np.any(p['w2']!=0): raise ValueError('snapshot hash')
    x=p['x'].astype(float)-3;y=p['y'].astype(float)-3
    delta=x*p['vy'].astype(float)-y*p['vx'].astype(float)-L0
    rows=[]
    for ix,neighbors in rings:
        n=len(ix); area=a['volume'][ix].astype(float)
        shape=area/np.mean(area)-1
        u=delta[ix]-np.mean(delta[ix])
        A=np.fft.rfft(shape)/n; U=np.fft.rfft(u)/n
        ps=np.abs(A)**2; pu=np.abs(U)**2
        ms=int(np.argmax(ps[1:])+1);mu=int(np.argmax(pu[1:])+1)
        baseline=1/(len(ps)-1)
        frac=float(ps[mu]/np.sum(ps[1:])) if np.sum(ps[1:])>0 else None
        rank=int(1+np.count_nonzero(ps[1:]>ps[mu]))
        rows.append({'n':n,'neighbor_counts':neighbors,'neighbor_count_differences':[abs(n-v) for v in neighbors],
                     'initial_area_peak_m':ms,'late_L_peak_m':mu,'initial_area_power_fraction_at_late_peak':frac,
                     'uniform_mode_fraction':baseline,'initial_area_power_rank_at_late_peak':rank,
                     'phase_difference_rad':float(np.angle(U[mu]*np.conj(A[mu])))})
    cases.append({'time':meta['t'],'source_sha256':digest,'rings':rows,
                  'matching_peak_count':sum(q['initial_area_peak_m']==q['late_L_peak_m'] for q in rows)})
result={'initial_sha256':late['initial_sha256'],'cases':cases,'script_sha256':sha(Path(__file__)),
        'scope':'Exploratory comparison of t0 stored cell-area ring spectra against late angular-momentum residual peaks in a bounded Q2 run. The late peaks are selected after evolution; mode rank and uniform-mode fractions are descriptive, not preregistered significance tests or evidence of causation. Cell areas and ring counts do not establish the physical forcing or temporal evolution; no spatial periodic boundary is used.'}
with (docs/'KEPLER_Q2_INITIAL_GEOMETRY_MODES_2026-10-29.json').open('x') as f:
    json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
