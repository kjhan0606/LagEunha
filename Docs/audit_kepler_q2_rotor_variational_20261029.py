"""Frozen Q2 diagnostic of a conservative geometric-inertia rotor candidate.

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
initial=Path(late['initial']);path=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1/khout.011058.dat')
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
spin,rot_energy,_,_=energies(*grad,qxx,qxy,qyy)
center=np.zeros_like(xy);center[active]=xy[active]+cl
coef=np.zeros(len(p));coef[active]=(mass[active]/area)*(spin*spin/(2*inertia*inertia))
i,j,c,length=finite_faces_touching_active(tri,centers,flag)
inner=(r0>=.5)&(r0<.75)
if np.any((inner[i]&~flag[j]) | (inner[j]&~flag[i])):raise ValueError('missing inner neighbor moments')
delta=xy[j]-xy[i];distance=np.linalg.norm(delta,axis=1);normal=delta/distance[:,None]
tangent=np.column_stack((-normal[:,1],normal[:,0]));mid=.5*(xy[i]+xy[j]);force=np.zeros_like(xy)
radius2=np.zeros(len(p));radius2[active]=inertia/mass[active]
for sign in (-1,1):
    q=c+sign*(length/(2*np.sqrt(3)))[:,None]*tangent
    hi=np.sum((q-center[i])**2,axis=1)-radius2[i];hj=np.sum((q-center[j])**2,axis=1)-radius2[j]
    delta_h=coef[i]*hi-coef[j]*hj;off=(q-mid)/distance[:,None]
    fi=.5*length[:,None]*delta_h[:,None]*(.5*normal+off)
    fj=.5*length[:,None]*delta_h[:,None]*(.5*normal-off)
    np.add.at(force,i,fi);np.add.at(force,j,fj)

# Check -dH/dq dot direction with spins held fixed; all cells touching the
# selected inertia set receive force, including neighboring nonactive cells.
rng=np.random.default_rng(20261029);direction=rng.normal(0,.01,size=xy.shape)
pred=-float(np.sum(force*direction,dtype=np.longdouble))
fd=[]
for h in (1e-4,5e-5):
    ip=moments(xy+h*direction);im=moments(xy-h*direction)
    ep=float(np.sum(spin*spin/(2*(ip[2]+ip[4])),dtype=np.longdouble))
    em=float(np.sum(spin*spin/(2*(im[2]+im[4])),dtype=np.longdouble))
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
        'thermal_energy_same_cells':float(np.sum(p['ie'][active],dtype=np.longdouble)),
        'main_band_rotor_torque_projection':float(projection/power),
        'main_band_native_pressure_torque_projection':ref['pressure_torque_projected_growth_rate_per_time'] if 'pressure_torque_projected_growth_rate_per_time' in ref else ref['pressure_torque_projected_growth_rate_by_band'][2],
        'finite_difference_predicted_energy_rate':pred,'finite_difference':fd,'script_sha256':sha(Path(__file__)),
        'scope':'Counterfactual Hamiltonian sum s^2/(2I) with fixed spin s and uniform-density geometric inertia, evaluated on saved Q2 states. Spin is reconstructed from stored velocity gradients and not an evolved native state. Force is a proposed artificial rotor term, not proven Euler hydrodynamics, stability, physical energy or a native solution; boundary of the active inertia set is held outside the sum. Gradient checked only for the tested direction.'}
out=docs/'KEPLER_Q2_ROTOR_VARIATIONAL_2026-10-29.json'
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
