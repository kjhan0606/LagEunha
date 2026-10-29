"""Frozen Q2 diagnostic of a fixed-momentum affine kinetic-energy gradient.

This is a counterfactual discrete Hamiltonian, not a physical Euler closure.
The analytic face gradient is checked against independent finite differences.
"""
import hashlib
import json
import socket
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ID_MASK,ordered_global
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from audit_kepler_bounded_pressure_force_t10_20261027 import finite_faces_touching_active
from audit_kepler_polygon_moment_native_faces_20261027 import oriented_faces
from audit_kepler_affine_spin_energy_20261027 import tensor_from_faces,energies

if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('run on compute host')
docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
late=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
prior=json.loads((docs/'KEPLER_Q2_MODE_TORQUE_2026-10-29.json').read_text())
initial=Path(late['initial']);path=Path(late['cases'][0]['endpoint'])
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
a,ma=ordered_global(initial);p,meta=ordered_global(path,a['u4if']&ID_MASK)
ref=next(v for v in prior['cases'] if v['time']==meta['t'])
if sha(initial)!=late['initial_sha256'] or sha(path)!=ref['source_sha256'] or ma['t']!=0 or np.any(p['w2']!=0) or not np.array_equal(a['mass'],p['mass']):raise ValueError('snapshot mismatch')
xy0=np.column_stack((a['x']-3,a['y']-3));xy=np.column_stack((p['x']-3,p['y']-3));r0=np.linalg.norm(xy0,axis=1)
active=np.flatnonzero(r0<2.35);flag=np.zeros(len(a),bool);flag[active]=True;lookup=np.full(len(a),-1,int);lookup[active]=np.arange(len(active))
mass=p['mass'].astype(float);zero=np.zeros(len(p))
def moments(pos):
    area,g,tri,centers=power_areas(pos,zero,active,lookup,return_geometry=True)
    if g['hidden_active'] or np.any(area<=0):raise ValueError('incomplete geometry')
    owner,u,v,offset,counts=oriented_faces(tri,centers,pos,active,lookup)
    area2,cl,qxx,qxy,qyy=tensor_from_faces(u,v,owner,offset,counts,mass[active])
    if np.max(abs(area2/area-1))>1e-8 or np.any(qxx*qyy-qxy*qxy<=0):raise ValueError('invalid moments')
    return area,cl,qxx,qxy,qyy,tri,centers
area,cl,qxx,qxy,qyy,tri,centers=moments(xy)
inertia=qxx+qyy;grad=[p['stress.'+k][active].astype(float) for k in ('gUxx','gUxy','gUyx','gUyy')]
spin,rot_energy,affine_energy,_=energies(*grad,qxx,qxy,qyy)
G=np.stack((np.stack((grad[0],grad[1]),axis=1),np.stack((grad[2],grad[3]),axis=1)),axis=1)
S=np.empty((len(active),2,2));S[:,0,0]=qxx;S[:,0,1]=S[:,1,0]=qxy;S[:,1,1]=qyy
# A = G S is the transpose of the centered momentum-cross moment. Hold A
# fixed during the hypothetical shape derivative, not the velocity gradient.
A=G@S
center=np.zeros_like(xy);center[active]=xy[active]+cl
rho=np.zeros(len(p));rho[active]=mass[active]/area
Gfull=np.zeros((len(p),2,2));Gfull[active]=G
kin_per_mass=np.zeros(len(p));kin_per_mass[active]=2*affine_energy/mass[active]
i,j,c,length=finite_faces_touching_active(tri,centers,flag)
inner=(r0>=.5)&(r0<.75)
if np.any((inner[i]&~flag[j]) | (inner[j]&~flag[i])):raise ValueError('missing inner neighbor moments')
delta=xy[j]-xy[i];distance=np.linalg.norm(delta,axis=1);normal=delta/distance[:,None]
tangent=np.column_stack((-normal[:,1],normal[:,0]));mid=.5*(xy[i]+xy[j]);force=np.zeros_like(xy)
for sign in (-1,1):
    q=c+sign*(length/(2*np.sqrt(3)))[:,None]*tangent
    ui=np.einsum('nij,nj->ni',Gfull[i],q-center[i]);uj=np.einsum('nij,nj->ni',Gfull[j],q-center[j])
    hi=.5*rho[i]*(np.sum(ui*ui,axis=1)-kin_per_mass[i])
    hj=.5*rho[j]*(np.sum(uj*uj,axis=1)-kin_per_mass[j])
    delta_h=hi-hj;off=(q-mid)/distance[:,None]
    fi=.5*length[:,None]*delta_h[:,None]*(.5*normal+off)
    fj=.5*length[:,None]*delta_h[:,None]*(.5*normal-off)
    np.add.at(force,i,fi);np.add.at(force,j,fj)

# Check -dK/dq with the momentum-cross tensor A held fixed. This fixed-A
# derivative alone is not a complete Hamiltonian or spin evolution law.
rng=np.random.default_rng(20261029);direction=rng.normal(0,.01,size=xy.shape)
pred=-float(np.sum(force*direction,dtype=np.longdouble))
fd=[]
for h in (1e-4,5e-5):
    ip=moments(xy+h*direction);im=moments(xy-h*direction)
    def kinetic(v):
        s=np.empty((len(active),2,2));s[:,0,0]=v[2];s[:,0,1]=s[:,1,0]=v[3];s[:,1,1]=v[4]
        return float(np.sum(.5*np.einsum('nij,nij->n',A@np.linalg.inv(s),A),dtype=np.longdouble))
    ep=kinetic(ip);em=kinetic(im)
    fd.append({'step':h,'energy_derivative':(ep-em)/(2*h),'absolute_difference':abs((ep-em)/(2*h)-pred)})
residual=(xy[:,0]*p['vy']-xy[:,1]*p['vx'])-(xy0[:,0]*a['vy']-xy0[:,1]*a['vx'])
torque=(xy[:,0]*force[:,1]-xy[:,1]*force[:,0])/mass
phi=np.arctan2(xy0[:,1],xy0[:,0]);rad=np.round(r0,10);power=0.;projection=0.
for radius in np.unique(rad[inner]):
    idx=np.flatnonzero(inner&(rad==radius));idx=idx[np.argsort(phi[idx])];n=len(idx)
    u=residual[idx]-np.mean(residual[idx]);v=torque[idx]-np.mean(torque[idx]);U=np.fft.rfft(u)/n;V=np.fft.rfft(v)/n
    factor=np.full(len(U),2.);factor[0]=1.
    if n%2==0:factor[-1]=1.
    frequency=np.fft.rfftfreq(n);band=(frequency>=.2)&(frequency<.3);weight=float(np.sum(mass[idx]))
    power+=weight*np.sum(factor[band]*abs(U[band])**2)
    projection+=weight*np.sum(factor[band]*np.real(np.conj(U[band])*V[band]))
if power<=0:raise ValueError('zero modal power')
result={'time':meta['t'],'snapshot_sha256':ref['source_sha256'],'sites_for_rotor_energy':len(active),
        'spin_energy_total':float(np.sum(rot_energy,dtype=np.longdouble)),
        'affine_kinetic_energy_total':float(np.sum(affine_energy,dtype=np.longdouble)),
        'thermal_energy_same_cells':float(np.sum(p['ie'][active],dtype=np.longdouble)),
        'main_band_fixed_affine_moment_torque_projection':float(projection/power),
        'main_band_native_pressure_torque_projection':ref['pressure_torque_projected_growth_rate_per_time'] if 'pressure_torque_projected_growth_rate_per_time' in ref else ref['pressure_torque_projected_growth_rate_by_band'][2],
        'finite_difference_predicted_energy_rate':pred,'finite_difference':fd,'script_sha256':sha(Path(__file__)),
        'scope':'Counterfactual derivative of affine kinetic energy with fixed centered momentum-cross tensor and uniform-density geometric second moments, evaluated on saved Q2 state. The tensor is reconstructed from stored gradients, not evolved. This isolated derivative is not a complete Hamiltonian, Euler closure, proof of angular momentum/energy conservation, stability or a native solution. Gradient checked only in one direction.'}
out=docs/'KEPLER_Q2_AFFINE_VARIATIONAL_T30_2026-10-29.json'
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
