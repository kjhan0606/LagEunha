"""Frozen-state pressure torque on uniform-density Q2 polygon centroids."""
import hashlib
import json
import socket
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ID_MASK, ordered_global
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from audit_kepler_bounded_pressure_force_t10_20261027 import finite_faces_touching_active
from audit_kepler_polygon_moment_native_faces_20261027 import oriented_faces
from audit_kepler_affine_spin_energy_20261027 import tensor_from_faces

if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):
    raise RuntimeError('Run on a compute host')
docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
late=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
prior=json.loads((docs/'KEPLER_Q2_NORMAL_TRACTION_2026-10-29.json').read_text())
root=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1')
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
initial=Path(late['initial']);a,meta=ordered_global(initial)
if meta['t']!=0 or sha(initial)!=late['initial_sha256']:raise ValueError('initial changed')
ids=a['u4if']&ID_MASK
xy0=np.column_stack((a['x'].astype(float)-3,a['y'].astype(float)-3))
r0=np.linalg.norm(xy0,axis=1)
inner=np.flatnonzero((r0>=.5)&(r0<.75))
flag=np.zeros(len(a),bool);flag[inner]=True
lookup=np.full(len(a),-1,int);lookup[inner]=np.arange(len(inner))
paths=[Path(late['cases'][0]['endpoint']),root/'khout.011058.dat']
cases=[]
for path in paths:
    p,meta=ordered_global(path,ids);digest=sha(path)
    ref=next((v for v in prior['cases'] if v['time']==meta['t']),None)
    if ref is None or ref['source_sha256']!=digest or not np.array_equal(p['mass'],a['mass']) or np.any(p['w2']!=0):raise ValueError('changed state')
    xy=np.column_stack((p['x'].astype(float)-3,p['y'].astype(float)-3))
    areas,g,tri,centers=power_areas(xy,np.zeros(len(a)),inner,lookup,return_geometry=True)
    if g['hidden_active']:raise ValueError('incomplete geometry')
    owner,start,end,offset,counts=oriented_faces(tri,centers,xy,inner,lookup)
    mass=p['mass'][inner].astype(float)
    area,centroid_local,qxx,qxy,qyy=tensor_from_faces(start,end,owner,offset,counts,mass)
    inertia=qxx+qyy
    area_error=float(np.max(np.abs(area/p['volume'][inner].astype(float)-1)))
    if area_error>1e-5 or np.any(qxx*qyy-qxy*qxy<=0):raise ValueError('bad cell geometry')
    centroids=np.zeros_like(xy);centroids[inner]=xy[inner]+centroid_local
    i,j,c,l=finite_faces_touching_active(tri,centers,flag)
    d=xy[j]-xy[i];norm=np.linalg.norm(d,axis=1);n=d/norm[:,None]
    pressure=.5*(p['pressure'][i].astype(float)+p['pressure'][j].astype(float))
    force=-pressure[:,None]*l[:,None]*n
    cross=lambda u,v:u[:,0]*v[:,1]-u[:,1]*v[:,0]
    ia=flag[i];ja=flag[j];ti=cross(c[ia]-centroids[i[ia]],force[ia]);tj=cross(c[ja]-centroids[j[ja]],-force[ja])
    net=np.zeros(len(inner));abs_faces=np.zeros(len(inner))
    np.add.at(net,lookup[i[ia]],ti);np.add.at(net,lookup[j[ja]],tj)
    np.add.at(abs_faces,lookup[i[ia]],abs(ti));np.add.at(abs_faces,lookup[j[ja]],abs(tj))
    both=ia&ja
    orbital=cross(centroids[i[both]]-centroids[j[both]],force[both])
    spin_i=cross(c[both]-centroids[i[both]],force[both]);spin_j=cross(c[both]-centroids[j[both]],-force[both])
    closure=float(np.max(abs(orbital+spin_i+spin_j)))
    scale=float(np.max(abs(orbital)+abs(spin_i)+abs(spin_j)))
    omega=(np.sum(xy[inner]**2,axis=1)+.01175**2)**(-.75)
    spin_accel=abs(net/inertia)/(omega*omega)
    cases.append({'time':meta['t'],'source_sha256':digest,'cells':len(inner),
                  'max_volume_relative_error':area_error,
                  'net_abs_spin_torque_over_sum_abs_face_torque':float(np.sum(abs(net))/np.sum(abs_faces)),
                  'spin_accel_over_omega2_mass_mean':float(np.average(spin_accel,weights=mass)),
                  'spin_accel_over_omega2_p99':float(np.quantile(spin_accel,.99)),
                  'pair_angular_momentum_closure_relative':closure/max(scale,1e-300)})
result={'initial_sha256':late['initial_sha256'],'cases':cases,'script_sha256':sha(Path(__file__)),
        'scope':'Instantaneous arithmetic-mean normal pressure traction on frozen Q2 W0 cells. Uniform density assumed to identify polygon centroid and inertia with mass moments; no evolved spin, actual mass-centroid state, Riemann traction, energy work, I-dot or evolution. Tau/I divided by Omega squared is a local scale, not measured spin evolution or relaxation.'}
out=docs/'KEPLER_Q2_PRESSURE_SPIN_2026-10-29.json'
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
