"""Compare close Q2 neighbors' angular speeds with local circular shear."""
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
out=docs/'KEPLER_Q2_CLOSE_SHEAR_2026-10-29.json'
if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('Run analysis on compute host')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
initial=Path(prior['initial'])
if sha(initial)!=prior['initial_sha256']:raise ValueError('Initial changed')
a,ma=ordered_global(initial);ids=a['u4if']&ID_MASK
xy0=np.column_stack((a['x'],a['y']));r0=np.hypot(xy0[:,0]-3.,xy0[:,1]-3.)
sel=np.flatnonzero((r0>=.5)&(r0<.75));mass=a['mass'][sel].astype(float)
root=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1')
paths=[root/'khout.011058.dat',Path(prior['cases'][1]['endpoint'])]
cases=[]
for path in paths:
    p,mt=ordered_global(path,ids);digest=sha(path)
    expected=next((v for v in history['cases'] if v['time']==mt['t']),None)
    if ma['t']!=0 or expected is None or digest!=expected['endpoint_sha256'] or not np.array_equal(a['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('Unexpected state')
    xy=np.column_stack((p['x']-3.,p['y']-3.));vel=np.column_stack((p['vx'],p['vy']))
    d,idx=cKDTree(xy).query(xy[sel],k=2,workers=1);j=idx[:,1];dist=d[:,1]
    exact=orbit(xy0,mt['t']);refdist=np.linalg.norm(exact[j]-exact[sel],axis=1);close=dist/refdist<.8
    radius=np.linalg.norm(xy,axis=1);omega=(xy[:,0]*vel[:,1]-xy[:,1]*vel[:,0])/(radius*radius)
    circular=(radius*radius+.01175**2)**(-.75)
    actual=omega[j]-omega[sel];expected_shear=circular[j]-circular[sel];residual=actual-expected_shear
    dr=radius[j]-radius[sel];w=mass[close];actual=actual[close];expected_shear=expected_shear[close];residual=residual[close];dr=dr[close]
    if not len(w) or np.any(~np.isfinite(actual)):raise ValueError('Invalid pair data')
    rms=lambda x:float(np.sqrt(mean(x*x,w)))
    covariance=mean(actual*expected_shear,w)
    correlation=float(covariance/(rms(actual)*rms(expected_shear)))
    cases.append({'time':mt['t'],'sites':len(w),'same_initial_ring_fraction':float(np.mean(abs(r0[j[close]]-r0[sel[close]])<1e-10)),
        'rms_radial_offset':rms(dr),'rms_actual_angular_speed_difference':rms(actual),
        'rms_circular_shear_prediction':rms(expected_shear),'rms_residual_angular_speed_difference':rms(residual),
        'uncentered_correlation':correlation,'same_sign_fraction':float(np.mean(actual*expected_shear>0)),
        'source_sha256':digest})
result={'initial_sha256':prior['initial_sha256'],'cases':cases,'script_sha256':sha(Path(__file__)),
        'scope':'Same closest neighbors of each selected site, counting mutual pairs twice, evaluated at individual snapshots. Circular prediction uses current radii; residual includes noncircular velocities. Instantaneous correlation is not a cause, trajectory or conserved-flux test.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
