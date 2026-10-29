"""Separate Q2 pressure force into adiabatic and entropy contributions."""
import hashlib
import json
import socket
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ID_MASK, ordered_global
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from audit_kepler_bounded_pressure_force_t10_20261027 import finite_faces_touching_active

if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):
    raise RuntimeError('Run analysis on compute host')
docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
late=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
history=json.loads((docs/'KEPLER_Q2_TORQUE_FEEDBACK_2026-10-29.json').read_text())
root=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1')

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

initial=Path(late['initial']); a,meta=ordered_global(initial)
if meta['t']!=0 or sha(initial)!=late['initial_sha256']:raise ValueError('changed initial snapshot')
ids=a['u4if']&ID_MASK
xy0=np.column_stack((a['x'].astype(float)-3,a['y'].astype(float)-3))
r0=np.linalg.norm(xy0,axis=1);phi=np.arctan2(xy0[:,1],xy0[:,0])
L0=xy0[:,0]*a['vy'].astype(float)-xy0[:,1]*a['vx'].astype(float)
disk=(r0>=.5)&(r0<.75);rounded=np.round(r0,10)
rings=[]
for r in np.unique(rounded[disk]):
    ix=np.flatnonzero(disk & (rounded==r));ix=ix[np.argsort(phi[ix])]
    gaps=np.diff(np.r_[phi[ix],phi[ix][0]+2*np.pi])
    if len(ix)<24 or np.max(abs(gaps/(2*np.pi/len(ix))-1))>1e-4:raise ValueError('unexpected ring')
    rings.append(ix)
active=np.flatnonzero(r0<2);flag=np.zeros(len(a),bool);flag[active]=True
lookup=np.full(len(a),-1,int);lookup[active]=np.arange(len(active))
paths=[Path(late['cases'][0]['endpoint']),root/'khout.011058.dat']
cases=[]
for path in paths:
    p,meta=ordered_global(path,ids);digest=sha(path)
    previous=next((q for q in history['cases'] if q['time']==meta['t']),None)
    if previous is None or previous['source_sha256']!=digest or not np.array_equal(a['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('snapshot contract')
    xy=np.column_stack((p['x'].astype(float)-3,p['y'].astype(float)-3))
    area,g,tri,centers=power_areas(xy,np.zeros(len(a)),active,lookup,return_geometry=True)
    error=float(np.max(np.abs(area/p['volume'][active].astype(float)-1)))
    if g['hidden_active'] or error>=1e-5:raise ValueError('mesh mismatch')
    i,j,c,l=finite_faces_touching_active(tri,centers,flag)
    delta=xy[j]-xy[i];dist=np.linalg.norm(delta,axis=1);n=delta/dist[:,None];off=c-.5*(xy[i]+xy[j])
    bi=l[:,None]*(.5*n+off/dist[:,None]);bj=l[:,None]*(.5*n-off/dist[:,None]);ii=flag[i];jj=flag[j]
    if any(np.any(~np.isfinite(v)) or np.any(v<=0) for v in (a['den'],a['pressure'],p['den'],p['pressure'])):raise ValueError('invalid state')
    gamma=float(np.float32(1.666667))
    actual=p['pressure'].astype(float)
    adiabatic=a['pressure'].astype(float)*(p['den'].astype(float)/a['den'].astype(float))**gamma
    def accel(pressure):
        dp=pressure[i]-pressure[j];f=np.zeros((len(p),2))
        np.add.at(f,i[ii],dp[ii,None]*bi[ii]);np.add.at(f,j[jj],dp[jj,None]*bj[jj])
        return f/p['mass'].astype(float)[:,None]
    acc=accel(actual);acc_iso=accel(adiabatic)
    torque=xy[:,0]*acc[:,1]-xy[:,1]*acc[:,0]
    torque_iso=xy[:,0]*acc_iso[:,1]-xy[:,1]*acc_iso[:,0]
    mass=p['mass'].astype(float);norm=lambda v:float(np.sqrt(np.average(np.sum(v[disk]**2,axis=1),weights=mass[disk])))
    actual_norm=norm(acc);iso_norm=norm(acc_iso);diff_norm=norm(acc-acc_iso)
    entropy_log=np.log(actual/adiabatic)

    residual=xy[:,0]*p['vy'].astype(float)-xy[:,1]*p['vx'].astype(float)-L0
    # Each band holds sum(m_ring * Parseval-normalized mode power) and
    # sum(m_ring * Re(conj(L_k) tau_k)), with Nyquist counted only once.
    power=np.zeros(5);projection=np.zeros(5);projection_iso=np.zeros(5)
    for ix in rings:
        u=residual[ix]-np.mean(residual[ix]);v=torque[ix]-np.mean(torque[ix]);N=len(ix)
        U=np.fft.rfft(u)/N;V=np.fft.rfft(v)/N
        v_iso=torque_iso[ix]-np.mean(torque_iso[ix]);V_iso=np.fft.rfft(v_iso)/N
        factor=np.full(len(U),2.0);factor[0]=1.0
        if N%2==0:factor[-1]=1.0
        mass=float(np.sum(a['mass'][ix].astype(float)))
        band=np.minimum((np.fft.rfftfreq(N)*10).astype(int),4)
        for k in range(5):
            m=band==k;power[k]+=mass*np.sum(factor[m]*np.abs(U[m])**2)
            projection[k]+=mass*np.sum(factor[m]*np.real(np.conj(U[m])*V[m]))
            projection_iso[k]+=mass*np.sum(factor[m]*np.real(np.conj(U[m])*V_iso[m]))
    if np.any(power<=0):raise ValueError('zero mode power')
    cases.append({'time':meta['t'],'source_sha256':digest,'max_volume_relative_error':error,
                  'pressure_torque_projected_growth_rate_by_band':(projection/power).tolist(),
                  'adiabatic_torque_projected_rate_by_band':(projection_iso/power).tolist(),
                  'acceleration_mass_rms_actual':actual_norm,'acceleration_mass_rms_adiabatic':iso_norm,'acceleration_mass_rms_difference':diff_norm,'acceleration_difference_over_actual':diff_norm/actual_norm,
                  'inner_log_entropy_mass_mean':float(np.average(entropy_log[disk],weights=p['mass'][disk])),
                  'inner_log_entropy_mass_rms':float(np.sqrt(np.average(entropy_log[disk]**2,weights=p['mass'][disk]))),
                  'mode_power_by_band':power.tolist()})
result={'initial_sha256':late['initial_sha256'],'rings':len(rings),'sites':int(np.sum(disk)),
        'bands_m_over_n':['[0,.1)','[.1,.2)','[.2,.3)','[.3,.4)','[.4,.5]'],
        'cases':cases,'script_sha256':sha(Path(__file__)),
        'scope':'Frozen bounded Q2 density and geometry: adiabatic pressure retains each site initial P/rho^gamma, compared to actual reconstructed pressure torque. Difference measures entropy contribution to this instantaneous counterfactual force, not physical entropy rate, independent evolution or causal stability proof.'}
out=docs/'KEPLER_Q2_ADIABATIC_TORQUE_2026-10-29.json'
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
