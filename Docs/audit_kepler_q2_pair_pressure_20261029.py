"""Reconstruct total pressure-acceleration differences of close Q2 pairs."""
import hashlib
import json
from pathlib import Path
import socket
import numpy as np
from scipy.spatial import cKDTree
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from audit_kepler_bounded_pressure_force_t10_20261027 import finite_faces_touching_active
from audit_kepler_t3_support_20261029 import orbit
from audit_kepler_native_radial_20261026 import mean

docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
prior=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
hist=json.loads((docs/'KEPLER_Q2_PRESSURE_HISTORY_2026-10-29.json').read_text())
out=docs/'KEPLER_Q2_PAIR_PRESSURE_2026-10-29.json'
if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('Run geometry on compute host')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
initial=Path(prior['initial'])
if sha(initial)!=prior['initial_sha256']:raise ValueError('Initial changed')
a,ma=ordered_global(initial);ids=a['u4if']&ID_MASK
xy0=np.column_stack((a['x'],a['y']));r0=np.hypot(xy0[:,0]-3.,xy0[:,1]-3.)
active=np.flatnonzero(r0<2.);flag=np.zeros(len(a),bool);flag[active]=True
lookup=np.full(len(a),-1,int);lookup[active]=np.arange(len(active))
sel=np.flatnonzero((r0>=.5)&(r0<.75));mass=a['mass'][sel].astype(float)
root=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1')
paths=[Path(prior['cases'][0]['endpoint']),root/'khout.011058.dat',Path(prior['cases'][1]['endpoint'])]
cases=[]
for path in paths:
    p,mt=ordered_global(path,ids);digest=sha(path)
    ref=next((r for r in hist['cases'] if r['time']==mt['t']),None)
    if ma['t']!=0 or ref is None or digest!=ref['endpoint_sha256'] or not np.array_equal(a['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('Unexpected state')
    xy=np.column_stack((p['x']-3.,p['y']-3.))
    area,g,tri,centers=power_areas(xy,np.zeros(len(a)),active,lookup,return_geometry=True)
    err=float(np.max(abs(area/p['volume'][active].astype(float)-1)))
    if g['hidden_active'] or err>=1e-5:raise ValueError('Bad geometry')
    i,j,c,l=finite_faces_touching_active(tri,centers,flag)
    delta=xy[j]-xy[i];d=np.linalg.norm(delta,axis=1)
    if np.any(d<=0):raise ValueError('Duplicate generators')
    n=delta/d[:,None];offset=c-.5*(xy[i]+xy[j]);dp=p['pressure'][i].astype(float)-p['pressure'][j].astype(float)
    bi=l[:,None]*(.5*n+offset/d[:,None]);bj=l[:,None]*(.5*n-offset/d[:,None])
    f=np.zeros((len(p),2));ii=flag[i];jj=flag[j]
    np.add.at(f,i[ii],dp[ii,None]*bi[ii]);np.add.at(f,j[jj],dp[jj,None]*bj[jj])
    accel=f/p['mass'].astype(float)[:,None]
    dist,neighbors=cKDTree(xy).query(xy[sel],k=2,workers=1);peer=neighbors[:,1];dactual=dist[:,1]
    if np.any(~flag[peer]) or np.any(dactual<=0):raise ValueError('Unreconstructed nearest neighbor')
    exact=orbit(xy0,mt['t']);dex=np.linalg.norm(exact[peer]-exact[sel],axis=1)
    normal=(xy[peer]-xy[sel])/dactual[:,None]
    arel=np.sum((accel[peer]-accel[sel])*normal,axis=1)
    speed=np.sum((np.column_stack((p['vx'],p['vy']))[peer]-np.column_stack((p['vx'],p['vy']))[sel])*normal,axis=1)
    ratio=dactual/dex;close=ratio<.8
    radial=xy[sel]/np.linalg.norm(xy[sel],axis=1)[:,None]
    alignment=abs(np.sum(normal*radial,axis=1))
    groups={}
    for label,mask in [('close',close),('other',~close)]:
        if not np.any(mask):continue
        groups[label]={'sites':int(mask.sum()),'relative_pressure_accel_mass_mean':mean(arel[mask],mass[mask]),
          'relative_pressure_accel_negative_fraction':float(np.mean(arel[mask]<0)),
          'relative_speed_mass_mean':mean(speed[mask],mass[mask]),
          'approaching_fraction':float(np.mean(speed[mask]<0)),
          'same_initial_ring_fraction':float(np.mean(abs(r0[sel[mask]]-r0[peer[mask]])<1e-10)),
          'radial_alignment_median':float(np.median(alignment[mask])),
          'radial_dominant_fraction':float(np.mean(alignment[mask]>2**-.5))}
    cases.append({'time':mt['t'],'source_sha256':digest,'max_volume_relative_error':err,'groups':groups})
result={'initial_sha256':prior['initial_sha256'],'cases':cases,'script_sha256':sha(Path(__file__)),
        'scope':'Instantaneous difference of independently reconstructed TOTAL W0 pressure accelerations along nearest-pair line, not an individual-face force. Negative means accelerations drive the separation downward at that instant; gravity/velocity curvature affect the full separation acceleration. No direct native Q2 RHS validation, time integral, AV attribution or causal proof. Sites can count a mutual pair twice.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
