"""Least-squares material-face weight rate on a saved W0 Kepler state."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import LinearOperator,cg

from audit_kepler_bounded_pressure_force_t10_20261027 import finite_faces_touching_active
from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest


N=55351
EPS=.01175
DT=.00193983


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('case',choices=('aligned','jitter'))
    ap.add_argument('receipt',type=Path)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    receipt=json.loads(args.receipt.read_text())
    if abs(float(receipt['time'])-6.5)>1e-10 or not all(receipt['checks'].values()):
        raise ValueError('Expected accepted native t=6.5 state')
    for key in ('initial','endpoint'):
        if digest(Path(receipt[key]))!=receipt[key+'_sha256']:
            raise ValueError(f'Changed {key} snapshot')
    initial,_=ordered(Path(receipt['initial']),N)
    p,_=ordered(Path(receipt['endpoint']),N)
    if np.any(p['w2']!=0) or not np.array_equal(initial['mass'],p['mass']):
        raise ValueError('Expected W0 fixed-mass state')
    xy0=np.column_stack((initial['x'].astype(float)-3.,initial['y'].astype(float)-3.))
    xy=np.column_stack((p['x'].astype(float)-3.,p['y'].astype(float)-3.))
    vel=np.column_stack((p['vx'].astype(float),p['vy'].astype(float)))
    r0=np.linalg.norm(xy0,axis=1)
    active=np.flatnonzero((r0>=.5)&(r0<2.))
    flag=np.zeros(N,dtype=bool);flag[active]=True
    lookup=np.full(N,-1,dtype=int);lookup[active]=np.arange(len(active))
    area,meta,triangles,centers=power_areas(xy,np.zeros(N),active,lookup,return_geometry=True)
    if meta['hidden_active'] or np.any(area<=0):
        raise ValueError('Invalid active Voronoi geometry')
    area_error=float(np.max(np.abs(area/p['volume'][active].astype(float)-1.)))
    if area_error>1e-5:
        raise ValueError('Independent area does not match native snapshot')
    i,j,c,length=finite_faces_touching_active(triangles,centers,flag)
    inside=flag[i]&flag[j]
    i,j,c,length=i[inside],j[inside],c[inside],length[inside]
    delta=xy[j]-xy[i]
    distance=np.linalg.norm(delta,axis=1)
    if np.any(distance<=0):
        raise ValueError('Duplicate site pair')
    normal=delta/distance[:,None]
    w0=(np.einsum('ij,ij->i',xy[j],vel[j])-np.einsum('ij,ij->i',xy[i],vel[i])
        -np.einsum('ij,ij->i',vel[j]-vel[i],c))/distance
    radius2=np.einsum('ij,ij->i',c,c)
    fluid=np.column_stack((-c[:,1],c[:,0]))*(radius2+EPS**2)[:,None]**(-.75)
    slip=np.einsum('ij,ij->i',fluid,normal)-w0
    ni=lookup[i];nj=lookup[j]
    e=len(i);n=len(active)
    rows=np.repeat(np.arange(e),2)
    cols=np.column_stack((ni,nj)).ravel()
    coef=1./(2*distance)
    vals=np.column_stack((coef,-coef)).ravel()
    B=coo_matrix((vals,(rows,cols)),shape=(e,n)).tocsr()
    sqrt_length=np.sqrt(length)
    C=B.multiply(sqrt_length[:,None]).tocsr()
    laplacian=(C.T@C).tocsr()
    components,_=connected_components(laplacian,directed=False)
    if components!=1:
        raise ValueError(f'Disk face graph has {components} components')
    rhs=np.asarray(C.T@(sqrt_length*slip)).ravel()
    reduced=laplacian[1:,1:].tocsr()
    diag=reduced.diagonal()
    if np.any(diag<=0):
        raise ValueError('Invalid pinned graph Laplacian')
    preconditioner=LinearOperator(reduced.shape,matvec=lambda x:x/diag,dtype=float)
    count=[0]
    def tick(_): count[0]+=1
    sol,info=cg(reduced,rhs[1:],M=preconditioner,rtol=1e-8,atol=0.,
                maxiter=2500,callback=tick)
    rate=np.r_[0.,sol]
    corrected=slip-np.asarray(B@rate).ravel()
    grad=np.asarray(C.T@(sqrt_length*corrected)).ravel()
    gradient_ratio=float(np.linalg.norm(grad)/max(np.linalg.norm(rhs),1e-300))
    def rms(v): return float(np.sqrt(np.average(v*v,weights=length)))
    def meanabs(v): return float(np.average(np.abs(v),weights=length))
    density=p['den'].astype(float)
    mass=p['mass'].astype(float)
    def flux_stats(v):
        edge=.5*(density[i]+density[j])*length*v
        net=np.zeros(n)
        np.add.at(net,ni,edge);np.add.at(net,nj,-edge)
        frac=DT*np.abs(net)/mass[active]
        return {'gross_mass_crossing_over_disk_mass':float(DT*np.sum(np.abs(edge),dtype=np.longdouble)
                    /np.sum(mass[active],dtype=np.longdouble)),
                'net_cell_mass_fraction_p99':float(np.quantile(frac,.99)),
                'net_cell_mass_fraction_mean':float(np.average(frac,weights=mass[active]))}
    delta_rate=rate[ni]-rate[nj]
    pressure=p['pressure'].astype(float)
    face_work=(pressure[i]-pressure[j])*length*delta_rate/(2*distance)
    thermal=float(np.sum(p['ie'][active].astype(float),dtype=np.longdouble))
    checks={'graph_connected':components==1,'cg_converged':info==0,
            'least_squares_gradient':gradient_ratio<1e-6,
            'objective_not_increased':rms(corrected)<=rms(slip)*(1+1e-12),
            'finite_rate':bool(np.all(np.isfinite(rate)))}
    result={'case':args.case,'time':6.5,'receipt_sha256':digest(args.receipt),
            'endpoint_sha256':receipt['endpoint_sha256'],'disk_sites':n,'disk_disk_faces':e,
            'reference_power_area_error_max':area_error,'cg_iterations':count[0],
            'cg_info':int(info),'least_squares_gradient_relative':gradient_ratio,
            'length_weighted_rms_slip_w0':rms(slip),'length_weighted_rms_slip_optimal':rms(corrected),
            'length_weighted_mean_abs_slip_w0':meanabs(slip),
            'length_weighted_mean_abs_slip_optimal':meanabs(corrected),
            'rms_ratio_optimal_over_w0':rms(corrected)/rms(slip),
            'flux_proxy_w0':flux_stats(slip),'flux_proxy_optimal':flux_stats(corrected),
            'max_one_step_abs_delta_weight_over_d2':float(np.max(DT*np.abs(delta_rate)/distance**2)),
            'disk_disk_signed_weight_pressure_work_dt_over_internal_energy':float(
                DT*np.sum(face_work,dtype=np.longdouble)/thermal),
            'disk_disk_absolute_weight_pressure_work_dt_over_internal_energy':float(
                DT*np.sum(np.abs(face_work),dtype=np.longdouble)/thermal),
            'checks':{k:bool(v) for k,v in checks.items()},
            'scope':'One-state linear least-squares weight-rate lower bound on disk-disk '
                    'analytic Kepler face-normal slip, pinned gauge and fixed generator velocity. '
                    'Disk boundary faces, actual Riemann states, pressure-work constraints, '
                    'finite-step power geometry and mass/momentum/energy/spin flux evolution '
                    'are excluded; this is not a native hydrodynamic solution.'}
    with args.output.open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:result[k] for k in ('case','cg_iterations','cg_info',
        'rms_ratio_optimal_over_w0','flux_proxy_w0','flux_proxy_optimal','checks')},indent=2))
    if not all(checks.values()):
        raise SystemExit('Optimal material-face rate audit failed')


if __name__=='__main__':
    main()
