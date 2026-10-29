"""Equal-time bounded Q2 ring-phase control using preserved native states."""
import hashlib
import json
import socket
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ID_MASK, ordered_global

if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):
    raise RuntimeError('run on compute host')
docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
triad=json.loads((docs/'KEPLER_Q2_THREE_IC_NATIVE_T65_2026-10-27.json').read_text())
prior=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
ring_ref=json.loads((docs/'KEPLER_Q2_INITIAL_PRESSURE_FORCE_RING_2026-10-27.json').read_text())

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''): h.update(chunk)
    return h.hexdigest()

sources={'aligned':(Path(prior['initial']),prior['initial_sha256']),
         'ring_coherent_5pct':(Path(ring_ref['source_snapshot']),ring_ref['source_snapshot_sha256'])}
cases={}
for name,(initial,initial_hash) in sources.items():
    endpoint=Path(triad['cases'][name]['endpoint'])
    if sha(initial)!=initial_hash or sha(endpoint)!=triad['cases'][name]['endpoint_sha256']:
        raise ValueError('hash mismatch')
    a,ma=ordered_global(initial);p,mp=ordered_global(endpoint,a['u4if']&ID_MASK)
    if ma['t']!=0 or mp['t']!=6.5 or np.any(p['w2']!=0) or not np.array_equal(a['mass'],p['mass']):
        raise ValueError('time, W0, or fixed-mass contract')
    x=a['x'].astype(float)-3;y=a['y'].astype(float)-3
    r=np.hypot(x,y);phi=np.arctan2(y,x);rounded=np.round(r,10)
    selected=(r>=.5)&(r<.75);L0=x*a['vy'].astype(float)-y*a['vx'].astype(float)
    xe=p['x'].astype(float)-3;ye=p['y'].astype(float)-3
    delta=xe*p['vy'].astype(float)-ye*p['vx'].astype(float)-L0
    rows=[]
    for radius in np.unique(rounded[selected]):
        ix=np.flatnonzero(selected & (rounded==radius));ix=ix[np.argsort(phi[ix])];n=len(ix)
        gaps=np.diff(np.r_[phi[ix],phi[ix][0]+2*np.pi])
        if n<24 or np.max(abs(gaps/(2*np.pi/n)-1))>1e-4: raise ValueError('not a uniform ring')
        area=a['volume'][ix].astype(float);area=area/np.mean(area)-1
        u=delta[ix]-np.mean(delta[ix]);A=np.fft.rfft(area)/n;U=np.fft.rfft(u)/n
        ps=np.abs(A)**2;pu=np.abs(U)**2
        m=int(np.argmax(pu[1:])+1)
        rows.append({'n':n,'area_peak':int(np.argmax(ps[1:])+1),'L_peak':m,
                     'L_rms':float(np.sqrt(np.mean(u*u))),
                     'initial_area_rms':float(np.sqrt(np.mean(area*area))),
                     'initial_area_phase_at_L_peak':float(np.angle(U[m]*np.conj(A[m]))),
                     'initial_area_rank_at_L_peak':int(1+np.count_nonzero(ps[1:]>ps[m]))})
    cases[name]={'initial_sha256':initial_hash,'endpoint_sha256':sha(endpoint),'rings':rows}
if [r['n'] for r in cases['aligned']['rings']]!=[r['n'] for r in cases['ring_coherent_5pct']['rings']]:
    raise ValueError('ring counts do not match')
result={'time':6.5,'cases':cases,'script_sha256':sha(Path(__file__)),
        'scope':'Equal-time exploratory comparison of aligned vs 5% coherent ring-phase shifted W0 native endpoints; 11 matched original inner-disk rings. Different ICs have different discrete cell masses. t6.5 is not a test of late t30/t47.5 mode growth, and selected peaks and phases do not prove causal seeding. Angular DFT is ring indexing, not a spatial periodic boundary.'}
with (docs/'KEPLER_Q2_COHERENT_SEED_MODES_2026-10-29.json').open('x') as f:
    json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
