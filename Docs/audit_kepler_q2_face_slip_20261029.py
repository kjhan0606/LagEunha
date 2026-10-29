"""Compare exact Voronoi face speed with saved affine velocity reconstruction."""
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
prior=json.loads((docs/'KEPLER_Q2_MODE_TORQUE_2026-10-29.json').read_text())
root=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1')
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
initial=Path(late['initial']);a,meta=ordered_global(initial)
if meta['t']!=0 or sha(initial)!=late['initial_sha256']:raise ValueError('initial changed')
ids=a['u4if']&ID_MASK
xy0=np.column_stack((a['x'].astype(float)-3,a['y'].astype(float)-3));r0=np.linalg.norm(xy0,axis=1)
inner=np.flatnonzero((r0>=.5)&(r0<.75));flag=np.zeros(len(a),bool);flag[inner]=True
lookup=np.full(len(a),-1,int);lookup[inner]=np.arange(len(inner))
cases=[]
for path in [Path(late['cases'][0]['endpoint']),root/'khout.011058.dat']:
    p,meta=ordered_global(path,ids);digest=sha(path)
    ref=next((v for v in prior['cases'] if v['time']==meta['t']),None)
    if ref is None or ref['source_sha256']!=digest or not np.array_equal(p['mass'],a['mass']) or np.any(p['w2']!=0):raise ValueError('changed state')
    xy=np.column_stack((p['x'].astype(float)-3,p['y'].astype(float)-3));vel=np.column_stack((p['vx'],p['vy'])).astype(float)
    area,g,tri,centers=power_areas(xy,np.zeros(len(a)),inner,lookup,return_geometry=True)
    if g['hidden_active'] or np.max(abs(area/p['volume'][inner]-1))>1e-5:raise ValueError('invalid geometry')
    i,j,c,l=finite_faces_touching_active(tri,centers,flag)
    both=flag[i]&flag[j];i=i[both];j=j[both];c=c[both];l=l[both]
    dr=xy[j]-xy[i];dist=np.linalg.norm(dr,axis=1);n=dr/dist[:,None];offset=c-.5*(xy[i]+xy[j])
    # Derivative of (x_face - midpoint).(x_j-x_i)=0.
    w=np.sum(.5*(vel[i]+vel[j])*n,axis=1)-np.sum((vel[j]-vel[i])*offset,axis=1)/dist
    def face_v(k):
        rel=c-xy[k]
        return np.column_stack((vel[k,0]+p['stress.gUxx'][k]*rel[:,0]+p['stress.gUxy'][k]*rel[:,1],
                                vel[k,1]+p['stress.gUyx'][k]*rel[:,0]+p['stress.gUyy'][k]*rel[:,1]))
    vi=face_v(i);vj=face_v(j)
    candidate=np.sum(.5*(vi+vj)*n,axis=1)
    cs_all=np.sqrt(float(np.float32(1.666667))*p['pressure'].astype(float)/p['den'].astype(float))
    cs=.5*(cs_all[i]+cs_all[j])
    if np.any(cs<=0) or not np.all(np.isfinite(candidate)):raise ValueError('invalid fields')
    slip=np.abs(candidate-w);ratio=slip/cs
    cases.append({'time':meta['t'],'source_sha256':digest,'interior_faces':len(i),
                  'length_weighted_mean_abs_slip':float(np.average(slip,weights=l)),
                  'length_weighted_mean_slip_over_cs':float(np.average(ratio,weights=l)),
                  'face_p95_slip_over_cs':float(np.quantile(ratio,.95)),
                  'length_fraction_slip_over_cs_gt_1':float(np.sum(l[ratio>1])/np.sum(l))})
result={'initial_sha256':late['initial_sha256'],'cases':cases,'script_sha256':sha(Path(__file__)),
        'scope':'Frozen W0 Voronoi face-normal geometry speed vs arithmetic mean of saved-gradient reconstructed velocities. Both sites initially in .5<=r<.75. A kinematic indicator only, not the native Riemann contact velocity, actual mass flux, proof of fixed-mass incompatibility, or a spin evolution. Stored gradients may belong to a nearby RK stage. Shared polygon domain is bounded.'}
out=docs/'KEPLER_Q2_FACE_SLIP_2026-10-29.json'
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
