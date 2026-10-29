"""Fit one scalar weight rate per generator to Kepler material face speeds."""
import hashlib
import json
import socket
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import lsmr
from kepler_native_global_ids_20261029 import ID_MASK, ordered_global
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from audit_kepler_bounded_pressure_force_t10_20261027 import finite_faces_touching_active

if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('run on compute host')
docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
late=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
history=json.loads((docs/'KEPLER_Q2_TORQUE_FEEDBACK_2026-10-29.json').read_text())
root=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return h.hexdigest()
initial=Path(late['initial']);a,meta=ordered_global(initial)
if meta['t']!=0 or sha(initial)!=late['initial_sha256']:raise ValueError('initial hash')
ids=a['u4if']&ID_MASK;r0=np.hypot(a['x']-3,a['y']-3)
active=np.flatnonzero((r0>=.5)&(r0<.75));flag=np.zeros(len(a),bool);flag[active]=True
lookup=np.full(len(a),-1,int);lookup[active]=np.arange(len(active))
cases=[]
for path in (Path(late['cases'][0]['endpoint']),root/'khout.011058.dat'):
    p,meta=ordered_global(path,ids);ref=next(v for v in history['cases'] if v['time']==meta['t'])
    if sha(path)!=ref['source_sha256'] or np.any(p['w2']!=0) or not np.array_equal(a['mass'],p['mass']):raise ValueError('snapshot changed')
    xy=np.column_stack((p['x'].astype(float)-3,p['y'].astype(float)-3))
    area,info,tri,centers=power_areas(xy,np.zeros(len(p)),active,lookup,return_geometry=True)
    area_error=float(np.max(abs(area/p['volume'][active].astype(float)-1)))
    if info['hidden_active'] or area_error>=1e-5:raise ValueError('geometry mismatch')
    i,j,c,length=finite_faces_touching_active(tri,centers,flag)
    delta=xy[j]-xy[i];distance=np.linalg.norm(delta,axis=1);n=delta/distance[:,None]
    tangent=np.column_stack((-n[:,1],n[:,0]));mid=.5*(xy[i]+xy[j])
    def velocity(z):
        om=(np.sum(z*z,axis=1)+.01175**2)**(-.75)
        return np.column_stack((-z[:,1]*om,z[:,0]*om))
    vel=velocity(xy);slips=[]
    for sign in (-1,1):
        point=c+sign*(length/(2*np.sqrt(3)))[:,None]*tangent
        geom=np.sum(.5*(vel[i]+vel[j])*n,axis=1)-np.sum((vel[j]-vel[i])*(point-mid),axis=1)/distance
        slips.append(np.sum(velocity(point)*n,axis=1)-geom)
    slip=np.stack(slips,axis=1)
    # One row per Gauss point. Interior pair rates enter with opposite signs.
    cols=np.unique(np.r_[i,j]);colmap=np.full(len(p),-1,int);colmap[cols]=np.arange(len(cols))
    M=len(i);rows=np.repeat(np.arange(2*M),2)
    columns=np.column_stack((np.repeat(colmap[i],2),np.repeat(colmap[j],2))).ravel()
    coeff=np.repeat(np.sqrt(length/2)/(2*distance),2)
    values=np.column_stack((coeff,-coeff)).ravel()
    J=coo_matrix((values,(rows,columns)),shape=(2*M,len(cols))).tocsr()
    b=(np.sqrt(length/2)[:,None]*slip).ravel()
    rigid=np.column_stack((-xy[:,1],xy[:,0]));rigid_max=0.
    for sign in (-1,1):
        point=c+sign*(length/(2*np.sqrt(3)))[:,None]*tangent
        speed=np.sum(.5*(rigid[i]+rigid[j])*n,axis=1)-np.sum((rigid[j]-rigid[i])*(point-mid),axis=1)/distance
        fluid=np.sum(np.column_stack((-point[:,1],point[:,0]))*n,axis=1)
        rigid_max=max(rigid_max,float(np.max(abs(fluid-speed))))
    gauge_error=float(np.max(abs(J@np.ones(len(cols)))))
    probe=np.random.default_rng(20261029).normal(size=len(cols))
    expected=np.repeat(np.sqrt(length/2)*(probe[colmap[i]]-probe[colmap[j]])/(2*distance),2)
    jacobian_error=float(np.max(abs(J@probe-expected)))
    if rigid_max>1e-12 or gauge_error>1e-12 or jacobian_error>1e-12:raise ValueError('rigid/gauge/Jacobian control failed')
    norms=np.sqrt(np.asarray(J.power(2).sum(axis=0)).ravel())
    if np.any(norms<=0) or not np.linalg.norm(b)>0:raise ValueError('degenerate fit')
    B=J@diags(1/norms);fit=lsmr(B,b,atol=1e-10,btol=1e-10,maxiter=10000)
    residual=B@fit[0]-b
    normal=float(np.linalg.norm(B.T@residual)/(np.linalg.norm(B.T@b)+1e-300))
    rates=fit[0]/norms;mean_slip=np.mean(slip,axis=1)
    correction=(rates[colmap[i]]-rates[colmap[j]])/(2*distance)
    mean_residual=correction-mean_slip
    face_mean_relative=float(np.linalg.norm(np.sqrt(length)*mean_residual)/np.linalg.norm(np.sqrt(length)*mean_slip))
    face_mean_power=float(np.sum(length*mean_slip**2)/np.dot(b,b))
    cell_before=np.zeros(len(active));cell_after=np.zeros(len(active))
    for sites,side in ((i,1.),(j,-1.)):
        mask=flag[sites];loc=lookup[sites[mask]]
        np.add.at(cell_before,loc,side*length[mask]*mean_slip[mask])
        np.add.at(cell_after,loc,side*length[mask]*mean_residual[mask])
    cell_net_ratio=float(np.linalg.norm(cell_after/area)/np.linalg.norm(cell_before/area))
    # A separate free constant per face is a lower bound, since a single
    # power weight per cell imposes graph-cycle compatibility.
    no_graph=(np.sqrt(length/2)[:,None]*(slip-np.mean(slip,axis=1)[:,None])).ravel()
    cases.append({'time':meta['t'],'snapshot_sha256':ref['source_sha256'],'inner_cells':len(active),'faces':M,
        'weight_variables':len(cols),'area_relative_error_max':area_error,'rigid_rotation_max_abs_slip':rigid_max,
        'constant_weight_gauge_max_abs':gauge_error,'jacobian_probe_max_abs':jacobian_error,
        'least_squares_istop':int(fit[1]),'least_squares_iterations':int(fit[2]),
        'relative_face_slip_residual':float(np.linalg.norm(residual)/np.linalg.norm(b)),
        'independent_face_constant_lower_bound':float(np.linalg.norm(no_graph)/np.linalg.norm(b)),
        'face_mean_slip_relative_residual':face_mean_relative,
        'face_mean_fraction_of_total_slip_squared':face_mean_power,
        'cell_net_area_rate_slip_relative_rms':cell_net_ratio,
        'normal_equation_relative_residual':normal,
        'rms_baseline_slip_per_sqrt_face_length':float(np.linalg.norm(b)/np.sqrt(np.sum(length))),
        'abs_weight_rate_over_median_area_omega_p95':float(np.quantile(abs(fit[0]/norms)/(np.median(area)*np.median((np.sum(xy[active]**2,axis=1)+.01175**2)**(-.75))),.95))})
result={'initial_sha256':sha(initial),'cases':cases,'script_sha256':sha(Path(__file__)),
    'scope':'At saved W0 Q2 positions, use exact Kepler velocities for generators and fluid; fit unconstrained instantaneous signed scalar weight rates at all inner and touching neighbor sites to normal-speed equality at two Gauss points per face. W=0 is not maintained if rates are nonzero. This is a frozen linearized no-mass-flux material-face test, not variable-mass/pressure/energy closure, finite-weight evolution, a rigorous singular-value proof, or long-time stability.'}
with (docs/'KEPLER_Q2_MATERIAL_FACE_WEIGHT_FIT_2026-10-29.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
