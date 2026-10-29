"""Compare saved weighted/W0 Q2 inner-disk nearest-pair momentum and torque."""
import hashlib
import json
from pathlib import Path
import socket
import numpy as np
from scipy.spatial import cKDTree
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from audit_kepler_bounded_pressure_force_t10_20261027 import finite_faces_touching_active
from audit_kepler_native_radial_20261026 import mean

docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
w0=json.loads(Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/to10_near32_a10_v1/audit_t10.json').read_text())
weighted=json.loads((docs/'KEPLER_Q2_WEIGHTED_T10_EQUILIBRIUM_2026-10-27/aligned_corrected_audit.json').read_text())
out=docs/'KEPLER_Q2_WEIGHTED_T10_PAIR_TORQUE_2026-10-29.json'
if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('Run geometry on compute host')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
if w0['initial_sha256']!=weighted['initial_sha256'] or sha(Path(w0['initial']))!=w0['initial_sha256']:raise ValueError('Different initial conditions')
a,ma=ordered_global(Path(w0['initial']));ids=a['u4if']&ID_MASK
r0=np.hypot(a['x']-3.,a['y']-3.);sel=np.flatnonzero((r0>=.5)&(r0<.75));m=a['mass'][sel].astype(float)
active=np.flatnonzero(r0<2.);flag=np.zeros(len(a),bool);flag[active]=True
lookup=np.full(len(a),-1,int);lookup[active]=np.arange(len(active))
baseline,mb=ordered_global(Path(w0['endpoint']),ids)
peer=cKDTree(np.column_stack((baseline['x'],baseline['y']))).query(np.column_stack((baseline['x'][sel],baseline['y'][sel])),k=2,workers=1)[1][:,1]
if np.any(~flag[peer]):raise ValueError('Unreconstructed baseline neighbor')
cases={}
for name,rec in [('w0',w0),('weighted',weighted)]:
    path=Path(rec['endpoint'])
    if sha(path)!=rec['endpoint_sha256']:raise ValueError('Changed endpoint')
    p,mt=ordered_global(path,ids)
    if ma['t']!=0 or mt['t']!=10. or not np.array_equal(a['mass'],p['mass']):raise ValueError('Unexpected state')
    if name=='w0' and np.any(p['w2']!=0):raise ValueError('Baseline not W0')
    xy=np.column_stack((p['x']-3.,p['y']-3.));weights=p['w2'].astype(float)
    area,g,tri,centers=power_areas(xy,weights,active,lookup,return_geometry=True)
    err=float(np.max(abs(area/p['volume'][active].astype(float)-1)))
    if g['hidden_active'] or err>=1e-5:raise ValueError('Power geometry mismatch')
    i,j,c,l=finite_faces_touching_active(tri,centers,flag)
    delta=xy[j]-xy[i];dist=np.linalg.norm(delta,axis=1)
    if np.any(dist<=0):raise ValueError('Coincident generators')
    normal=delta/dist[:,None];offset=c-.5*(xy[i]+xy[j]);dp=p['pressure'][i].astype(float)-p['pressure'][j].astype(float)
    bi=l[:,None]*(.5*normal+offset/dist[:,None]);bj=l[:,None]*(.5*normal-offset/dist[:,None])
    f=np.zeros((len(p),2));ii=flag[i];jj=flag[j]
    np.add.at(f,i[ii],dp[ii,None]*bi[ii]);np.add.at(f,j[jj],dp[jj,None]*bj[jj])
    acc=f/p['mass'].astype(float)[:,None]
    angular=xy[:,0]*p['vy']-xy[:,1]*p['vx'];torque=xy[:,0]*acc[:,1]-xy[:,1]*acc[:,0]
    dl=angular[peer]-angular[sel];dt=torque[peer]-torque[sel]
    rms=lambda x:float(np.sqrt(mean(x*x,m)))
    cases[name]={'cells':len(sel),'rms_specific_L_difference':rms(dl),'rms_specific_torque_difference':rms(dt),
                 'positive_torque_L_split_site_fraction':float(np.mean(dl*dt>0)),
                 'mean_L_tau_over_mean_L_squared_per_time':mean(dl*dt,m)/mean(dl*dl,m),
                 'absolute_density_change_mass_mean':mean(abs(p['den'][sel].astype(float)/a['den'][sel].astype(float)-1),m),
                 'max_relative_area_error':err,'max_abs_weight':float(np.max(abs(weights))), 'endpoint_sha256':rec['endpoint_sha256']}
result={'initial_sha256':w0['initial_sha256'],'time':10.,'cases':cases,'script_sha256':sha(Path(__file__)),
        'scope':'Same initial aligned Q2 and same W0-selected nearest pairs at t10; weighted path enabled only after t6.501. Instantaneous reconstructed variational forces, not direct weighted native RHS, pairwise first-law validation, or long-time stability. Fixed masses; global weighted-work projection does not establish local conservation.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
