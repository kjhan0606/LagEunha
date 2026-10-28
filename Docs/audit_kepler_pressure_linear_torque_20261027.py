"""Compare exact linear-pressure and midpoint face moments on native polygons."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_affine_spin_energy_20261027 import tensor_from_faces
from audit_kepler_polygon_moment_native_faces_20261027 import oriented_faces
from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest

N=55351


def local_pressure_gradients(xy,p,triangles,ids,lookup):
    edges=np.concatenate((triangles[:,(0,1)],triangles[:,(1,2)],triangles[:,(2,0)]))
    edges.sort(axis=1);edges=np.unique(edges,axis=0)
    active=lookup>=0
    left=edges[active[edges[:,0]]];right=edges[active[edges[:,1]]]
    owner=np.concatenate((left[:,0],right[:,1]));neighbor=np.concatenate((left[:,1],right[:,0]))
    local=lookup[owner];delta=xy[neighbor]-xy[owner];dp=p[neighbor]-p[owner]
    d2=np.einsum('ij,ij->i',delta,delta)
    if np.any(d2<=0):raise ValueError('Duplicate neighbor')
    w=1./d2;n=len(ids)
    accum=lambda x:np.bincount(local,weights=x,minlength=n)
    xx=accum(w*delta[:,0]**2);xy_m=accum(w*delta[:,0]*delta[:,1]);yy=accum(w*delta[:,1]**2)
    bx=accum(w*delta[:,0]*dp);by=accum(w*delta[:,1]*dp);det=xx*yy-xy_m*xy_m
    if np.any(det<=1e-12):raise ValueError('Degenerate pressure fit')
    return np.column_stack(((yy*bx-xy_m*by)/det,(xx*by-xy_m*bx)/det))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('receipt',type=Path)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    receipt=json.loads(args.receipt.read_text())
    if not all(receipt['checks'].values()):raise ValueError('Native gate failed')
    for key in ('initial','endpoint'):
        if digest(Path(receipt[key]))!=receipt[key+'_sha256']:raise ValueError(f'Changed {key}')
    p0,m0=ordered(Path(receipt['initial']),N);p,m=ordered(Path(receipt['endpoint']),N)
    if abs(m0['t'])>1e-12 or abs(m['t']-receipt['time'])>1e-10:raise ValueError('Time mismatch')
    xy0=np.column_stack((p0['x'].astype(float)-3.,p0['y'].astype(float)-3.))
    xy=np.column_stack((p['x'].astype(float)-3.,p['y'].astype(float)-3.))
    r0=np.linalg.norm(xy0,axis=1);ids=np.flatnonzero((r0>=.5)&(r0<2.));lookup=np.full(N,-1,dtype=int);lookup[ids]=np.arange(len(ids))
    area,meta,triangles,centers=power_areas(xy,p['w2'].astype(float),ids,lookup,return_geometry=True)
    if meta['hidden_active'] or np.any(area<=0):raise ValueError('Invalid power cells')
    owner,a,b,offset,counts=oriented_faces(triangles,centers,xy,ids,lookup)
    mass=p['mass'][ids].astype(float);area2,centroid,qxx,qxy,qyy=tensor_from_faces(a,b,owner,offset,counts,mass)
    gradient=local_pressure_gradients(xy,p['pressure'].astype(float),triangles,ids,lookup)
    e=b-a;normal=np.column_stack((e[:,1],-e[:,0]));mid=.5*(a+b);r=mid-centroid[owner]
    pressure=p['pressure'][ids].astype(float)[owner]+np.sum(gradient[owner]*mid,axis=1)
    cross=r[:,0]*normal[:,1]-r[:,1]*normal[:,0]
    midpoint_face=-pressure*cross
    corr_face=-(np.sum(gradient[owner]*e,axis=1)*(e[:,0]*normal[:,1]-e[:,1]*normal[:,0]))/12.
    starts=offset[:-1];midpoint=np.add.reduceat(midpoint_face,starts)
    correction=np.add.reduceat(corr_face,starts);exact=midpoint+correction
    face_abs=np.add.reduceat(np.abs(midpoint_face)+np.abs(corr_face),starts)
    rel_exact=np.abs(exact)/np.maximum(face_abs,1e-300)
    radius=np.linalg.norm(xy[ids],axis=1);omega2=(radius*radius+.01175**2)**(-1.5)
    scale=(qxx+qyy)*omega2;normalized=np.abs(midpoint)/scale
    checks={'tensor_positive':bool(np.all(qxx*qyy-qxy*qxy>0)),
            'area_consistent':bool(np.max(np.abs(area2/area-1.))<1e-9),
            'linear_pressure_exact_torque_zero':bool(np.max(rel_exact)<1e-9),
            'finite':bool(np.all(np.isfinite(normalized)))}
    bands=[]
    for lo,hi in ((.5,.75),(.75,1.),(1.,1.5),(1.5,2.)):
        mask=(r0[ids]>=lo)&(r0[ids]<hi)
        bands.append({'r_min':lo,'r_max':hi,'cells':int(mask.sum()),
                      'normalized_midpoint_torque_rms':float(np.sqrt(np.average(normalized[mask]**2,weights=mass[mask]))),
                      'normalized_midpoint_torque_p95':float(np.quantile(normalized[mask],.95))})
    result={'receipt':str(args.receipt.resolve()),'receipt_sha256':digest(args.receipt),
            'endpoint_sha256':receipt['endpoint_sha256'],'time':receipt['time'],
            'weighted':bool(np.any(p['w2']!=0)),'disk_sites':len(ids),
            'pressure_gradient':'independent inverse-distance weighted least-squares fit to neighboring saved cell pressures',
            'exact_torque_relative_to_face_torque_scale_max':float(np.max(rel_exact)),
            'normalized_midpoint_torque_rms':float(np.sqrt(np.average(normalized**2,weights=mass))),
            'normalized_midpoint_torque_p95':float(np.quantile(normalized,.95)),
            'radial_bands':bands,'checks':checks,
            'scope':'Hypothetical independent linear-pressure normal traction on reconstructed '
                    'uniform-density polygons. Midpoint versus exact face integration is a '
                    'geometric diagnostic, not the native variational force, physical spin '
                    'source, accumulated torque, or a heating measurement.'}
    with args.output.open('x') as out:json.dump(result,out,indent=2,allow_nan=False);out.write('\n')
    print(json.dumps({'time':result['time'],'weighted':result['weighted'],
                      'normalized_rms':result['normalized_midpoint_torque_rms'],'checks':checks},indent=2))
    if not all(checks.values()):raise SystemExit('Linear-pressure torque audit failed')


if __name__=='__main__':main()
