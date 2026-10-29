"""Fit weight rates to the material/geometric second-moment mismatch."""
import hashlib
import json
import socket
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import lsmr
from kepler_native_global_ids_20261029 import ID_MASK,ordered_global
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from audit_kepler_bounded_pressure_force_t10_20261027 import finite_faces_touching_active
from audit_kepler_polygon_moment_native_faces_20261027 import oriented_faces
from audit_kepler_affine_spin_energy_20261027 import tensor_from_faces

if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('run on compute host')
docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
late=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
prior=json.loads((docs/'KEPLER_Q2_MODE_TORQUE_2026-10-29.json').read_text())
root=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1')
initial=Path(late['initial'])
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
a,ma=ordered_global(initial);ids=a['u4if']&ID_MASK
if sha(initial)!=late['initial_sha256'] or ma['t']!=0:raise ValueError('initial mismatch')
r0=np.hypot(a['x']-3,a['y']-3);active=np.flatnonzero((r0>=.5)&(r0<.75))
flag=np.zeros(len(a),bool);flag[active]=True;lookup=np.full(len(a),-1,int);lookup[active]=np.arange(len(active))
mass=a['mass'][active].astype(float);zero=np.zeros(len(a))
def geom(xy,weights=None):
    area,g,tri,centers=power_areas(xy,zero if weights is None else weights,active,lookup,return_geometry=True)
    if g['hidden_active'] or np.any(area<=0):raise ValueError('invalid geometry')
    owner,u,v,offset,counts=oriented_faces(tri,centers,xy,active,lookup)
    area2,cl,qxx,qxy,qyy=tensor_from_faces(u,v,owner,offset,counts,mass)
    if np.max(abs(area2/area-1))>1e-8 or np.any(qxx*qyy-qxy*qxy<=0):raise ValueError('invalid moments')
    S=np.empty((len(active),2,2));S[:,0,0]=qxx;S[:,0,1]=S[:,1,0]=qxy;S[:,1,1]=qyy
    return area,xy[active]+cl,S,tri,centers
def norm(v):return np.sqrt(np.sum(v*v,axis=(1,2)))
cases=[]
for path in (Path(late['cases'][0]['endpoint']),root/'khout.011058.dat'):
    p,meta=ordered_global(path,ids);ref=next(v for v in prior['cases'] if v['time']==meta['t'])
    if sha(path)!=ref['source_sha256'] or np.any(p['w2']!=0) or not np.array_equal(a['mass'],p['mass']):raise ValueError('state mismatch')
    xy=np.column_stack((p['x']-3,p['y']-3));vel=np.column_stack((p['vx'],p['vy'])).astype(float)
    area,C,S,tri,centers=geom(xy)
    G=np.stack((np.stack((p['stress.gUxx'][active],p['stress.gUxy'][active]),axis=1),np.stack((p['stress.gUyx'][active],p['stress.gUyy'][active]),axis=1)),axis=1)
    material=G@S+S@np.swapaxes(G,1,2)
    i,j,c,length=finite_faces_touching_active(tri,centers,flag)
    delta=xy[j]-xy[i];d=np.linalg.norm(delta,axis=1);n=delta/d[:,None];tan=np.column_stack((-n[:,1],n[:,0]));mid=.5*(xy[i]+xy[j]);dS=np.zeros_like(S)
    for sign in (-1,1):
        q=c+sign*(length/(2*np.sqrt(3)))[:,None]*tan
        wn=np.sum(.5*(vel[i]+vel[j])*n,axis=1)-np.sum((vel[j]-vel[i])*(q-mid),axis=1)/d
        for sites,side in ((i,1.),(j,-1.)):
            selected=flag[sites];local=lookup[sites[selected]];z=q[selected]-C[local]
            shape=z[:,:,None]*z[:,None,:]-S[local]/mass[local,None,None]
            values=side*(.5*length[selected]*wn[selected]*mass[local]/area[local])[:,None,None]*shape
            np.add.at(dS,local,values)
    # Directly recompute geometric moments at moved generators, same masses.
    checks=[]
    for h in (2e-5,1e-5):
        plus=geom(xy+h*vel)[2];minus=geom(xy-h*vel)[2];fd=(plus-minus)/(2*h)
        checks.append({'step':h,'relative_frobenius_rms':float(np.linalg.norm(fd-dS)/np.linalg.norm(fd))})
    omega=(np.sum(xy[active]*xy[active],axis=1)+.01175**2)**(-.75)
    scale=omega*np.trace(S,axis1=1,axis2=2)
    difference=norm(dS-material)/scale;geom_norm=norm(dS)/scale;material_norm=norm(material)/scale
    # For weights W_i, w_n,weight=(Wdot_i-Wdot_j)/(2d). For a
    # constant face speed the exact segment average of zz^T is
    # z_mid z_mid^T + length^2 tangent tangent^T/12.
    cols=np.unique(np.r_[i,j]);colmap=np.full(len(a),-1,int);colmap[cols]=np.arange(len(cols))
    rows=[];columns=[];values=[];normalizer=np.sqrt(mass/np.mean(mass))/scale
    for sites,other in ((i,j),(j,i)):
        selected=flag[sites];local=lookup[sites[selected]];z=c[selected]-C[local];t=tan[selected];L=length[selected]
        mat=z[:,:,None]*z[:,None,:]+(L*L/12)[:,None,None]*t[:,:,None]*t[:,None,:]-S[local]/mass[local,None,None]
        coeff=(mass[local]/area[local]*L/(2*d[selected])*normalizer[local])[:,None]
        comps=np.column_stack((mat[:,0,0],np.sqrt(2)*mat[:,0,1],mat[:,1,1]))*coeff
        for k in range(3):
            rr=3*local+k;v=comps[:,k]
            rows.extend((rr,rr));columns.extend((colmap[sites[selected]],colmap[other[selected]]));values.extend((v,-v))
    J=coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(columns))),shape=(3*len(active),len(cols))).tocsr()
    target=(material-dS)*normalizer[:,None,None]
    b=np.column_stack((target[:,0,0],np.sqrt(2)*target[:,0,1],target[:,1,1])).reshape(-1)
    column_norm=np.sqrt(np.asarray(J.power(2).sum(axis=0)).ravel())
    if np.any(column_norm<=0):raise ValueError('zero Jacobian column')
    preconditioned=J@diags(1/column_norm)
    probe=np.random.default_rng(20261029).normal(0,1e-4,len(cols));full_probe=np.zeros(len(a));full_probe[cols]=probe
    hweight=1e-4;fdweight=(geom(xy,hweight*full_probe)[2]-geom(xy,-hweight*full_probe)[2])/(2*hweight)
    weighted_fd=fdweight*normalizer[:,None,None]
    vfd=np.column_stack((weighted_fd[:,0,0],np.sqrt(2)*weighted_fd[:,0,1],weighted_fd[:,1,1])).reshape(-1)
    jacobian_error=float(np.linalg.norm(J@probe-vfd)/np.linalg.norm(vfd))
    fit=lsmr(preconditioned,b,atol=1e-10,btol=1e-10,maxiter=1500)
    rates=fit[0]/column_norm;residual=J@rates-b
    scaled_rate=rates/(np.median(area)*np.median(omega))
    cases.append({'time':meta['t'],'snapshot_sha256':ref['source_sha256'],'cells':len(active),
                  'geometry_area_relative_max':float(np.max(abs(area/p['volume'][active].astype(float)-1))),
                  'geometric_rate_over_omega_trace_mass_mean':float(np.average(geom_norm,weights=mass)),
                  'material_rate_over_omega_trace_mass_mean':float(np.average(material_norm,weights=mass)),
                  'difference_over_omega_trace_mass_mean':float(np.average(difference,weights=mass)),
                  'difference_over_omega_trace_p95':float(np.quantile(difference,.95)),
                  'weight_variables':len(cols),'least_squares_istop':int(fit[1]),'least_squares_iterations':int(fit[2]),
                  'weight_jacobian_probe_relative_rms':jacobian_error,
                  'weight_fit_relative_residual':float(np.linalg.norm(residual)/np.linalg.norm(b)),
                  'weighted_normal_equation_relative_residual':float(np.linalg.norm(preconditioned.T@residual)/(np.linalg.norm(preconditioned.T@b)+1e-300)),
                  'abs_weight_rate_over_median_area_omega_median':float(np.median(abs(scaled_rate))),
                  'abs_weight_rate_over_median_area_omega_p95':float(np.quantile(abs(scaled_rate),.95)),
                  'finite_difference_checks':checks})
result={'initial_sha256':sha(initial),'cases':cases,'script_sha256':sha(Path(__file__)),
        'scope':'Unconstrained least-squares instantaneous signed weight-rate fit to fixed-mass uniform-density Q2 inner-disk second-moment mismatch; individual exterior-neighbor weight rates are free. No weight positivity, mass flux, geometry feasibility, temporal evolution, pressure work, native solver, or performance is represented. Iterative residual is not a proof of the global minimum unless convergence is established.'}
out=docs/'KEPLER_Q2_WEIGHT_MOMENT_FIT_2026-10-29.json'
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
