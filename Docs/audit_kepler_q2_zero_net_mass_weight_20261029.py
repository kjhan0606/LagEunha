"""Fit weight rates to zero net ALE mass flux in frozen uniform-density Q2."""
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

if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('run on compute host')
docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
late=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
history=json.loads((docs/'KEPLER_Q2_TORQUE_FEEDBACK_2026-10-29.json').read_text())
root=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
initial=Path(late['initial']);a,meta=ordered_global(initial)
if meta['t']!=0 or sha(initial)!=late['initial_sha256']:raise ValueError('initial hash')
ids=a['u4if']&ID_MASK;r0=np.hypot(a['x']-3,a['y']-3)
active=np.flatnonzero((r0>=.5)&(r0<.75));flag=np.zeros(len(a),bool);flag[active]=True
lookup=np.full(len(a),-1,int);lookup[active]=np.arange(len(active));cases=[]
for path in (Path(late['cases'][0]['endpoint']),root/'khout.011058.dat'):
    p,meta=ordered_global(path,ids);ref=next(v for v in history['cases'] if v['time']==meta['t'])
    if sha(path)!=ref['source_sha256'] or np.any(p['w2']!=0) or not np.array_equal(a['mass'],p['mass']):raise ValueError('state changed')
    xy=np.column_stack((p['x'].astype(float)-3,p['y'].astype(float)-3));area,g,tri,centers=power_areas(xy,np.zeros(len(p)),active,lookup,return_geometry=True)
    err=float(np.max(abs(area/p['volume'][active].astype(float)-1)))
    if g['hidden_active'] or err>=1e-5:raise ValueError('geometry mismatch')
    i,j,c,length=finite_faces_touching_active(tri,centers,flag);dxy=xy[j]-xy[i];dist=np.linalg.norm(dxy,axis=1)
    normal=dxy/dist[:,None];tan=np.column_stack((-normal[:,1],normal[:,0]));mid=.5*(xy[i]+xy[j])
    def velocity(z):
        om=(np.sum(z*z,axis=1)+.01175**2)**(-.75)
        return np.column_stack((-z[:,1]*om,z[:,0]*om))
    v=velocity(xy);nodes,weights=np.polynomial.legendre.leggauss(4);slip=np.empty((len(i),4))
    for k,node in enumerate(nodes):
        q=c+(.5*node*length)[:,None]*tan
        w=np.sum(.5*(v[i]+v[j])*normal,axis=1)-np.sum((v[j]-v[i])*(q-mid),axis=1)/dist
        slip[:,k]=np.sum(velocity(q)*normal,axis=1)-w
    mean=slip@(.5*weights)
    cols=np.unique(np.r_[i,j]);colmap=np.full(len(p),-1,int);colmap[cols]=np.arange(len(cols))
    norm=area*(np.sum(xy[active]**2,axis=1)+.01175**2)**(-.75)
    mass=p['mass'][active].astype(float);factor=np.sqrt(mass/np.mean(mass))/norm
    rows=[];columns=[];values=[];target=np.zeros(len(active))
    for sites,side in ((i,1.),(j,-1.)):
        mask=flag[sites];local=lookup[sites[mask]]
        coef=side*length[mask]/(2*dist[mask])*factor[local]
        rows.extend((local,local));columns.extend((colmap[i[mask]],colmap[j[mask]]));values.extend((coef,-coef))
        np.add.at(target,local,side*length[mask]*mean[mask]*factor[local])
    J=coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(columns))),shape=(len(active),len(cols))).tocsr()
    cnorm=np.sqrt(np.asarray(J.power(2).sum(axis=0)).ravel())
    if np.any(cnorm<=0) or not np.linalg.norm(target)>0:raise ValueError('degenerate fit')
    B=J@diags(1/cnorm);fit=lsmr(B,target,atol=1e-11,btol=1e-11,maxiter=10000)
    rates=fit[0]/cnorm;residual=B@fit[0]-target
    normalres=float(np.linalg.norm(B.T@residual)/(np.linalg.norm(B.T@target)+1e-300))
    correction=(rates[colmap[i]]-rates[colmap[j]])/(2*dist)
    post=slip-correction[:,None];sqrtw=np.sqrt(length[:,None]*.5*weights[None,:])
    raw=float(np.linalg.norm(sqrtw*slip));after=float(np.linalg.norm(sqrtw*post))
    mean_post=post@(.5*weights);raw_net=np.zeros(len(active));net=np.zeros(len(active))
    for sites,side in ((i,1.),(j,-1.)):
        mask=flag[sites];local=lookup[sites[mask]]
        np.add.at(raw_net,local,side*length[mask]*mean[mask]);np.add.at(net,local,side*length[mask]*mean_post[mask])
    identity=float(np.linalg.norm(net*factor+residual)/(np.linalg.norm(target)+1e-300))
    if identity>1e-10:raise ValueError('flux/Jacobian mismatch')
    cases.append({'time':meta['t'],'source_sha256':ref['source_sha256'],'cells':len(active),'faces':len(i),'weight_variables':len(cols),
        'area_relative_error_max':err,'lsmr_stop':int(fit[1]),'lsmr_iterations':int(fit[2]),'normal_equation_relative_residual':normalres,
        'net_mass_rate_relative_weighted_rms':float(np.linalg.norm(residual)/np.linalg.norm(target)),
        'flux_jacobian_identity_relative_rms':identity,'pointwise_slip_relative_rms':after/raw,
        'face_gross_exchange_before':float(np.sum(length*(abs(slip)@(.5*weights)))),
        'face_gross_exchange_after':float(np.sum(length*(abs(post)@(.5*weights)))),
        'face_mean_slip_rms_after_over_before':float(np.linalg.norm(np.sqrt(length)*mean_post)/np.linalg.norm(np.sqrt(length)*mean)),
        'rate_abs_over_median_area_omega_p95':float(np.quantile(abs(rates)/(np.median(area)*np.median((np.sum(xy[active]**2,axis=1)+.01175**2)**(-.75))),.95))})
result={'initial_sha256':sha(initial),'cases':cases,'script_sha256':sha(Path(__file__)),
        'scope':'Frozen Q2 W0 geometry, exact analytic Kepler velocities, uniform unit density; unconstrained Wdot fit minimizes per-cell net ALE mass flux over inner cells with touching-neighbor rates free. Nonzero Wdot leaves W0. Numerical face fluxes remain nonzero and require momentum/energy/moment transport. No Wdot law, nonuniform density, Riemann solver, finite-time W0 reduction, feasibility, positivity, native GPU or stability claimed.'}
with (docs/'KEPLER_Q2_ZERO_NET_MASS_WEIGHT_2026-10-29.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
