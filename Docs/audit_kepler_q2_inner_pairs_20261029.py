"""Measure nearest inner-disk generator pairs against a circular reference."""
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
prior=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
history=json.loads((docs/'KEPLER_Q2_PRESSURE_HISTORY_2026-10-29.json').read_text())
out=docs/'KEPLER_Q2_INNER_PAIRS_2026-10-29.json'
if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('Run analysis on compute host')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()
initial=Path(prior['initial'])
if sha(initial)!=prior['initial_sha256']:raise ValueError('Initial changed')
a,ma=ordered_global(initial);ids=a['u4if']&ID_MASK
xy0=np.column_stack((a['x'],a['y']));r=np.hypot(xy0[:,0]-3.,xy0[:,1]-3.)
sel=np.flatnonzero((r>=.5)&(r<.75));mass=a['mass'][sel].astype(float)
root=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1')
paths=[Path(prior['cases'][0]['endpoint']),root/'khout.005767.dat',root/'khout.011058.dat',Path(prior['cases'][1]['endpoint'])]
rows=[]
for path in paths:
    p,meta=ordered_global(path,ids);digest=sha(path)
    expected=next((v for v in history['cases'] if v['time']==meta['t']),None)
    if ma['t']!=0 or expected is None or digest!=expected['endpoint_sha256'] or not np.array_equal(a['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('Unexpected state')
    xy=np.column_stack((p['x'],p['y']));v=np.column_stack((p['vx'],p['vy']))
    dist,ix=cKDTree(xy).query(xy[sel],k=2,workers=1);j=ix[:,1];d=dist[:,1]
    exact=orbit(xy0,meta['t']);dex=np.linalg.norm(exact[j]-exact[sel],axis=1)
    if np.any(d<=0) or np.any(dex<=0):raise ValueError('Coincident points')
    normal=(xy[j]-xy[sel])/d[:,None]
    approach=np.sum((v[j]-v[sel])*normal,axis=1)
    omega=(r*r+.01175**2)**(-.75)
    ve=np.column_stack((-(exact[:,1]-3.)*omega,(exact[:,0]-3.)*omega))
    normal_e=(exact[j]-exact[sel])/dex[:,None]
    approach_e=np.sum((ve[j]-ve[sel])*normal_e,axis=1)
    ratio=d/dex;density=np.abs(p['den'][sel].astype(float)/a['den'][sel].astype(float)-1.)
    group=ratio<.8
    contrib=mass*density
    rows.append({'time':meta['t'],'cells':len(sel),'nearest_same_peer_distance_ratio_actual_to_circular':{
        'mass_mean':mean(ratio,mass),'cell_p1':float(np.percentile(ratio,1)),
        'cell_p10':float(np.percentile(ratio,10)),'cell_median':float(np.median(ratio))},
        'sites_with_ratio_below_0p8':int(group.sum()),
        'fraction_density_abs_norm_on_those_sites':float(np.sum(contrib[group])/np.sum(contrib)),
        'approaching_nearest_neighbor_fraction':float(np.mean(approach<0)),
        'nearest_neighbor_radial_relative_speed':{'mass_mean_signed':mean(approach,mass),'cell_p1':float(np.percentile(approach,1)),
                                                 'cell_p99':float(np.percentile(approach,99))},
        'same_pair_circular_relative_speed':{'mass_mean_signed':mean(approach_e,mass),
                                             'approaching_fraction':float(np.mean(approach_e<0)),
                                             'cell_p1':float(np.percentile(approach_e,1)),
                                             'cell_p99':float(np.percentile(approach_e,99))},
        'source_sha256':digest})
result={'initial_sha256':prior['initial_sha256'],'cases':rows,'script_sha256':sha(Path(__file__)),
        'scope':'Nearest actual generator per selected inner-disk site, distances compared with the same pair at exact circular positions. Sites can count the same pair twice. Velocity convergence includes expected Kepler shear; ratios and correlations alone do not diagnose AV, shock heating, or causal pairing.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
