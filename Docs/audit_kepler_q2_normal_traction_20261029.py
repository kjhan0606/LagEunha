"""Counterfactual shared normal-pressure traction on saved Q2 geometry."""
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
paths=[Path(late['cases'][0]['endpoint'])]+[root/('khout.'+n+'.dat') for n in ('001329','002762','004411','005767','007138','008700','011058')]
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
    delta=xy[j]-xy[i];dist=np.linalg.norm(delta,axis=1);n=delta/dist[:,None]
    # A single common, normal traction, equal and opposite per face.
    # Arithmetic face pressure is a diagnostic closure, not a Riemann solver.
    pface=.5*(p['pressure'][i].astype(float)+p['pressure'][j].astype(float))
    face_force=-pface[:,None]*l[:,None]*n
    force=np.zeros((len(p),2));ii=flag[i];jj=flag[j]
    np.add.at(force,i[ii],face_force[ii]);np.add.at(force,j[jj],-face_force[jj])
    acc=force/p['mass'].astype(float)[:,None]
    torque=xy[:,0]*acc[:,1]-xy[:,1]*acc[:,0]
    residual=xy[:,0]*p['vy'].astype(float)-xy[:,1]*p['vx'].astype(float)-L0
    # Each band holds sum(m_ring * Parseval-normalized mode power) and
    # sum(m_ring * Re(conj(L_k) tau_k)), with Nyquist counted only once.
    power=np.zeros(5);projection=np.zeros(5)
    for ix in rings:
        u=residual[ix]-np.mean(residual[ix]);v=torque[ix]-np.mean(torque[ix]);N=len(ix)
        U=np.fft.rfft(u)/N;V=np.fft.rfft(v)/N
        factor=np.full(len(U),2.0);factor[0]=1.0
        if N%2==0:factor[-1]=1.0
        mass=float(np.sum(a['mass'][ix].astype(float)))
        band=np.minimum((np.fft.rfftfreq(N)*10).astype(int),4)
        for k in range(5):
            m=band==k;power[k]+=mass*np.sum(factor[m]*np.abs(U[m])**2)
            projection[k]+=mass*np.sum(factor[m]*np.real(np.conj(U[m])*V[m]))
    if np.any(power<=0):raise ValueError('zero mode power')
    cases.append({'time':meta['t'],'source_sha256':digest,'max_volume_relative_error':error,
                  'normal_traction_torque_projection_rate_by_band':(projection/power).tolist(),
                  'mode_power_by_band':power.tolist()})
result={'initial_sha256':late['initial_sha256'],'rings':len(rings),'sites':int(np.sum(disk)),
        'bands_m_over_n':['[0,.1)','[.1,.2)','[.2,.3)','[.3,.4)','[.4,.5]'],
        'cases':cases,'script_sha256':sha(Path(__file__)),
        'scope':'Counterfactual instantaneous arithmetic-mean normal face pressure on frozen bounded W0 Q2 snapshots; not native evolution, not a Riemann state, energy work, spin flux, or validated method. Conservative pair momentum by construction for each interior face. Projection rates cannot predict the alternative trajectory or demonstrate stability.'}
out=docs/'KEPLER_Q2_NORMAL_TRACTION_2026-10-29.json'
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
