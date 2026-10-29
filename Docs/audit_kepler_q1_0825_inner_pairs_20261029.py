"""Compare Q1 nearest-neighbor and angular-momentum split at equal time."""
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
receipt=json.loads((docs/'KEPLER_Q1_GRAVITY_CAP_2026-10-29/ETA006_VS_UNCAPPED_0825.json').read_text())
initial=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/ic_v1/khout.000001.dat')
uncapped=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/chain_h200_v2/to0p825/khout.001218.dat')
capped=Path('/gpfs/kjhan/LagForce/Kepler_Q1_OrbitalCap006_Ada_20261029/to0p825/khout.010798.dat')
out=docs/'KEPLER_Q1_GRAVITY_CAP_2026-10-29/INNER_PAIRS_0825.json'
if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('Run on a compute host')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
if sha(initial)!=receipt['input_sha256']:raise ValueError('Initial changed')
a,ma=ordered_global(initial);ids=a['u4if']&ID_MASK
xy0=np.column_stack((a['x'],a['y']));r0=np.hypot(xy0[:,0]-3.,xy0[:,1]-3.)
sel=np.flatnonzero((r0>=.5)&(r0<.75));w=a['mass'][sel].astype(float)
cases={}
for name,path,key in [('uncapped',uncapped,'uncapped_endpoint_sha256'),('eta006',capped,'eta006_endpoint_sha256')]:
    if sha(path)!=receipt[key]:raise ValueError('Changed endpoint')
    p,mt=ordered_global(path,ids)
    if ma['t']!=0 or mt['t']!=.825 or not np.array_equal(a['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('Unexpected state')
    xy=np.column_stack((p['x'],p['y']));d,ix=cKDTree(xy).query(xy[sel],k=2,workers=1);peer=ix[:,1]
    exact=orbit(xy0,mt['t']);dr=np.linalg.norm(exact[peer]-exact[sel],axis=1);ratio=d[:,1]/dr
    close=ratio<.8;rel=xy-3.;specific_l=rel[:,0]*p['vy']-rel[:,1]*p['vx'];split=specific_l[peer]-specific_l[sel]
    error=np.linalg.norm(xy[sel]-exact[sel],axis=1)/np.sqrt(a['volume'][sel].astype(float))
    entropy=np.log(p['ie'][sel].astype(float)/a['ie'][sel].astype(float))+(float(np.float32(1.666667))-1.)*np.log(p['volume'][sel].astype(float)/a['volume'][sel].astype(float))
    cases[name]={'cells':len(sel),'sites_distance_ratio_lt_0p8':int(close.sum()),'distance_ratio_p1':float(np.percentile(ratio,1)),
                 'same_initial_ring_fraction':float(np.mean(abs(r0[sel]-r0[peer])<1e-10)),
                 'rms_specific_L_difference':float(np.sqrt(mean(split*split,w))),
                 'position_error_over_sqrt_initial_area_mass_mean':mean(error,w),
                 'absolute_entropy_change_mass_mean':mean(abs(entropy),w),'endpoint_sha256':sha(path)}
result={'time':.825,'cases':cases,'initial_sha256':receipt['input_sha256'],'script_sha256':sha(Path(__file__)),
        'scope':'Same-time Q1 W0 native endpoints. Each endpoint chooses its own nearest neighbors, which may differ. Sites may count mutual pairs twice. Specific-L differences are total, not isolated pressure torque. No long-time or thermal64 conclusion.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
