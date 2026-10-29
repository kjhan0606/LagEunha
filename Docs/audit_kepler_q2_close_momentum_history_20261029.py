"""Trace fixed late-close Q2 pairs' angular momentum and pressure torques."""
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
out=docs/'KEPLER_Q2_CLOSE_MOMENTUM_HISTORY_2026-10-29.json'
if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('Run analysis on compute host')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
initial=Path(prior['initial']);a,ma=ordered_global(initial);ids=a['u4if']&ID_MASK
if sha(initial)!=prior['initial_sha256'] or ma['t']!=0:raise ValueError('Initial changed')
xy0=np.column_stack((a['x'],a['y']));r0=np.hypot(xy0[:,0]-3.,xy0[:,1]-3.)
sel0=np.flatnonzero((r0>=.5)&(r0<.75));active=np.flatnonzero(r0<2.)
flag=np.zeros(len(a),bool);flag[active]=True
lookup=np.full(len(a),-1,int);lookup[active]=np.arange(len(active))
root=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1')
last=root/'khout.011058.dat';z,mz=ordered_global(last,ids)
xy=np.column_stack((z['x'],z['y']));d,ix=cKDTree(xy).query(xy[sel0],k=2,workers=1);peer0=ix[:,1]
exact=orbit(xy0,mz['t']);reference=np.linalg.norm(exact[peer0]-exact[sel0],axis=1)
close=d[:,1]/reference<.8;sel=sel0[close];peer=peer0[close];weights=a['mass'][sel].astype(float)
if len(sel)!=244 or np.any(~flag[peer]):raise ValueError('Unexpected selected pairs')
paths=[Path(prior['cases'][0]['endpoint'])]+[root/('khout.'+n+'.dat') for n in ('001329','002762','004411','005767','007138','008700','011058')]
def rms(q):return float(np.sqrt(mean(q*q,weights)))
rows=[];all_dl=[];all_torque=[];times=[]
for path in paths:
    p,m=ordered_global(path,ids);ref=next((v for v in hist['cases'] if v['time']==m['t']),None)
    if ref is None or sha(path)!=ref['endpoint_sha256'] or not np.array_equal(a['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('Invalid state')
    xy=np.column_stack((p['x']-3.,p['y']-3.));r=np.linalg.norm(xy,axis=1)
    area,g,tri,centers=power_areas(xy,np.zeros(len(a)),active,lookup,return_geometry=True)
    areaerr=float(np.max(abs(area/p['volume'][active].astype(float)-1)))
    if g['hidden_active'] or areaerr>=1e-5:raise ValueError('Invalid reconstructed geometry')
    i,j,c,l=finite_faces_touching_active(tri,centers,flag)
    delta=xy[j]-xy[i];dist=np.linalg.norm(delta,axis=1);n=delta/dist[:,None];offset=c-.5*(xy[i]+xy[j])
    dp=p['pressure'][i].astype(float)-p['pressure'][j].astype(float)
    bi=l[:,None]*(.5*n+offset/dist[:,None]);bj=l[:,None]*(.5*n-offset/dist[:,None])
    force=np.zeros((len(p),2));ii=flag[i];jj=flag[j]
    np.add.at(force,i[ii],dp[ii,None]*bi[ii]);np.add.at(force,j[jj],dp[jj,None]*bj[jj])
    acc=force/p['mass'].astype(float)[:,None]
    specific_l=xy[:,0]*p['vy']-xy[:,1]*p['vx'];torque=xy[:,0]*acc[:,1]-xy[:,1]*acc[:,0]
    dl=specific_l[peer]-specific_l[sel];dtau=torque[peer]-torque[sel]
    dr=r[peer]-r[sel];rm=.5*(r[peer]+r[sel]);s=rm*rm+.01175**2
    derivative=2*rm*s**(-.75)-1.5*rm**3*s**(-1.75)
    pred=dl/derivative;resid=dr-pred
    rows.append({'time':m['t'],'rms_pair_radial_difference':rms(dr),'rms_specific_angular_momentum_difference':rms(dl),
                 'rms_circular_angular_momentum_predicted_radial_difference':rms(pred),'rms_radial_residual':rms(resid),
                 'rms_pair_specific_pressure_torque_difference':rms(dtau),'max_volume_relative_error':areaerr,'source_sha256':ref['endpoint_sha256']})
    all_dl.append(dl);all_torque.append(dtau);times.append(m['t'])
impulse=np.zeros(len(sel))
for k in range(1,len(times)):
    impulse+=(times[k]-times[k-1])*.5*(all_torque[k-1]+all_torque[k])
actual=all_dl[-1]-all_dl[0];residual=actual-impulse
result={'initial_sha256':prior['initial_sha256'],'fixed_late_selected_sites':len(sel),'cases':rows,
        'coarse_trapezoid_t30_to_t47p5':{'rms_actual_pair_specific_L_change':rms(actual),
             'rms_integrated_reconstructed_pressure_torque':rms(impulse),'rms_residual':rms(residual)},
        'script_sha256':sha(Path(__file__)),
        'scope':'Fixed nearest pairs selected at t47.5, each site may count a mutual pair twice. Pressure torques are independently reconstructed at snapshots 2.5 time units apart, not direct native Q2 RHS or a resolved impulse; trapezoid can alias orbital-time dynamics. L/r relation is instantaneous, not a causal proof.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
