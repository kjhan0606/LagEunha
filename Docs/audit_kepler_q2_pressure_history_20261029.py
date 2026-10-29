"""Reconstruct Q2 pressure and orbit metrics at the preserved intermediate checkpoints."""
import hashlib
import json
from pathlib import Path
import socket
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from audit_kepler_bounded_pressure_force_t10_20261027 import finite_faces_touching_active
from audit_kepler_q1_t10_pressure_timescales_20261029 import describe

docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
prior=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
out=docs/'KEPLER_Q2_PRESSURE_HISTORY_2026-10-29.json'
if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('Run geometry on compute host')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
initial=Path(prior['initial'])
if sha(initial)!=prior['initial_sha256']:raise ValueError('Initial changed')
historical=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1')
paths=[historical/('khout.'+step+'.dat') for step in ('001329','002762','004411','005767','007138','008700','011058')]
for path in paths:
    _,sample=ordered_global(path)
    prior['cases'].append({'endpoint':str(path),'endpoint_sha256':sha(path),'time':sample['t']})
prior['cases'].sort(key=lambda c:c['time'])
a,ma=ordered_global(initial);ids=a['u4if']&ID_MASK
r0=np.hypot(a['x']-3.,a['y']-3.);active=np.flatnonzero(r0<2.)
flag=np.zeros(len(a),bool);flag[active]=True
lookup=np.full(len(a),-1,dtype=int);lookup[active]=np.arange(len(active))
cases=[]
for old in prior['cases']:
    path=Path(old['endpoint'])
    if sha(path)!=old['endpoint_sha256']:raise ValueError('Changed endpoint')
    p,meta=ordered_global(path,ids)
    if ma['t']!=0 or meta['t']!=old['time'] or not np.array_equal(a['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('Invalid input')
    xy=np.column_stack((p['x']-3.,p['y']-3.))
    area,geom,tri,centers=power_areas(xy,np.zeros(len(a)),active,lookup,return_geometry=True)
    err=float(np.max(abs(area/p['volume'][active].astype(float)-1)))
    if geom['hidden_active'] or np.any(area<=0) or err>=1e-5:raise ValueError('Invalid independent geometry')
    i,j,c,l=finite_faces_touching_active(tri,centers,flag)
    delta=xy[j]-xy[i];dist=np.linalg.norm(delta,axis=1)
    if np.any(dist<=0):raise ValueError('Coincident generators')
    normal=delta/dist[:,None];offset=c-.5*(xy[i]+xy[j]);pressure=p['pressure'].astype(float)
    bi=l[:,None]*(.5*normal+offset/dist[:,None]);bj=l[:,None]*(.5*normal-offset/dist[:,None])
    force=np.zeros((len(p),2));dv=np.zeros(len(p));ii=flag[i];jj=flag[j];dp=pressure[i]-pressure[j]
    np.add.at(force,i[ii],dp[ii,None]*bi[ii]);np.add.at(force,j[jj],dp[jj,None]*bj[jj])
    vel=np.column_stack((p['vx'],p['vy']))
    speed=.5*np.einsum('ij,ij->i',vel[i]+vel[j],normal)-np.einsum('ij,ij->i',vel[j]-vel[i],offset)/dist
    flux=l*speed;np.add.at(dv,i[ii],flux[ii]);np.add.at(dv,j[jj],-flux[jj])
    radius=np.hypot(xy[:,0],xy[:,1]);gravity=radius*(radius*radius+.01175**2)**(-1.5)
    omega=(radius*radius+.01175**2)**(-.75)
    accel=np.linalg.norm(force,axis=1)/p['mass'].astype(float)
    g_ratio=np.divide(accel,gravity,out=np.zeros(len(p)),where=gravity>0)
    thermal=np.abs(pressure*dv)/(p['ie'].astype(float)*omega);mass=p['mass'].astype(float)
    regions={}
    for name,lo,hi in [('background',.1,.5),('inner_disk',.5,.75),('disk',.5,2.)]:
        mask=(r0>=lo)&(r0<hi)
        if not np.all(gravity[mask]>0):raise ValueError('Undefined gravity ratio')
        phase=meta['t']*(r0*r0+.01175**2)**(-.75)
        xx=a['x']-3.;yy=a['y']-3.
        poserr=np.hypot(p['x']-(3.+xx*np.cos(phase)-yy*np.sin(phase)),p['y']-(3.+xx*np.sin(phase)+yy*np.cos(phase)))
        dens=np.abs(p['den'].astype(float)/a['den'].astype(float)-1.)
        regions[name]={'cells':int(mask.sum()),'pressure_acceleration':describe(accel[mask],mass[mask]),
                       'pressure_accel_over_gravity':describe(g_ratio[mask],mass[mask]),
                       'absolute_thermal_rate_over_Uomega':describe(thermal[mask],mass[mask]),
                       'density_abs_change':describe(dens[mask],mass[mask]),
                       'position_error_over_sqrt_initial_area':describe((poserr/np.sqrt(a['volume'].astype(float)))[mask],mass[mask])}
    cases.append({'time':meta['t'],'endpoint_sha256':old['endpoint_sha256'],'active_cells':len(active),
                  'faces_touching_active':len(i),'max_relative_area_error':err,'regions':regions})
result={'initial_sha256':prior['initial_sha256'],'cases':cases,'script_sha256':sha(Path(__file__)),
        'scope':'Independent W0 variational force and pressure work reconstruction at saved Q2 states, including intermediate times from checkpoint footers. The method was cross-checked against eight Q1 native RHS samples at t10, not against native Q2 RHS here. Ratios are not accumulated velocity impulses, heating, or causal proof.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
