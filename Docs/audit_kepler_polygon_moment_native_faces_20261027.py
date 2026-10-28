"""Compare polygon-moment kernel with independent triangle integration."""

import argparse
import ctypes
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest


N=55351


def oriented_faces(triangles,centers,xy,active,lookup):
    edges=np.vstack((triangles[:,[0,1]],triangles[:,[1,2]],triangles[:,[2,0]]))
    tids=np.tile(np.arange(len(triangles)),3)
    edges.sort(axis=1)
    order=np.lexsort((edges[:,1],edges[:,0]));edges=edges[order];tids=tids[order]
    starts=np.r_[0,1+np.flatnonzero(np.any(edges[1:]!=edges[:-1],axis=1))]
    count=np.diff(np.r_[starts,len(edges)])
    flag=np.zeros(N,dtype=bool);flag[active]=True
    touch=flag[edges[starts,0]]|flag[edges[starts,1]]
    if np.any(count[touch]!=2):raise ValueError('Nonfinite face touches active disk')
    starts=starts[touch]
    i,j=edges[starts,0],edges[starts,1]
    a=centers[tids[starts]].copy();b=centers[tids[starts+1]].copy()
    d=xy[j]-xy[i]
    outward=np.column_stack((b[:,1]-a[:,1],-(b[:,0]-a[:,0])))
    swap=np.einsum('ij,ij->i',outward,d)<0
    temp=a[swap].copy();a[swap]=b[swap];b[swap]=temp
    ia=flag[i];ja=flag[j]
    owner=np.r_[lookup[i[ia]],lookup[j[ja]]]
    first=np.vstack((a[ia]-xy[i[ia]],b[ja]-xy[j[ja]]))
    second=np.vstack((b[ia]-xy[i[ia]],a[ja]-xy[j[ja]]))
    order=np.argsort(owner,kind='stable')
    owner=owner[order];first=first[order];second=second[order]
    counts=np.bincount(owner,minlength=len(active))
    if np.any(counts<3):raise ValueError('Disk polygon has fewer than three finite faces')
    offset=np.r_[0,np.cumsum(counts)].astype(np.int32)
    return owner,first,second,offset,counts


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('receipt',type=Path)
    ap.add_argument('--expected-failed-check',action='append',default=[])
    ap.add_argument('--library',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    receipt=json.loads(args.receipt.read_text())
    failed={k for k,v in receipt['checks'].items() if not v}
    if failed!=set(args.expected_failed_check):raise ValueError(f'Unexpected gate failures: {sorted(failed)}')
    for key in ('initial','endpoint'):
        if digest(Path(receipt[key]))!=receipt[key+'_sha256']:
            raise ValueError(f'Changed {key} snapshot')
    initial,_=ordered(Path(receipt['initial']),N)
    p,meta_snapshot=ordered(Path(receipt['endpoint']),N)
    if abs(float(meta_snapshot['t'])-float(receipt['time']))>1e-10:
        raise ValueError('Checkpoint time mismatch')
    xy0=np.column_stack((initial['x'].astype(float)-3.,initial['y'].astype(float)-3.))
    xy=np.column_stack((p['x'].astype(float)-3.,p['y'].astype(float)-3.))
    r0=np.linalg.norm(xy0,axis=1)
    active=np.flatnonzero((r0>=.5)&(r0<2.))
    lookup=np.full(N,-1,dtype=int);lookup[active]=np.arange(len(active))
    weights=p['w2'].astype(float)
    area,meta,triangles,centers=power_areas(xy,weights,active,lookup,return_geometry=True)
    if meta['hidden_active'] or np.any(area<=0):raise ValueError('Invalid active power cells')
    owner,a,b,offset,counts=oriented_faces(triangles,centers,xy,active,lookup)
    n=len(active)
    x1=np.ascontiguousarray(a[:,0]);y1=np.ascontiguousarray(a[:,1])
    x2=np.ascontiguousarray(b[:,0]);y2=np.ascontiguousarray(b[:,1])
    mass=np.ascontiguousarray(p['mass'][active].astype(float))
    area_c=np.empty(n);cx_c=np.empty(n);cy_c=np.empty(n);inertia_c=np.empty(n)
    status=np.empty(n,dtype=np.int32)
    pd=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    pi=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
    lib=ctypes.CDLL(str(args.library.resolve()))
    fn=lib.gfs_polygon_moment_batch
    fn.argtypes=[ctypes.c_int,pi,pd,pd,pd,pd,pd,pd,pd,pd,pd,pi]
    fn.restype=ctypes.c_int
    code=fn(n,offset,x1,y1,x2,y2,mass,area_c,cx_c,cy_c,inertia_c,status)
    if code!=0 or np.any(status!=0):raise ValueError(f'Moment routine failed: code={code}, bad={np.count_nonzero(status)}')
    starts=offset[:-1]
    anchor_x=np.add.reduceat(x1,starts)/counts
    anchor_y=np.add.reduceat(y1,starts)/counts
    ax=anchor_x[owner];ay=anchor_y[owner]
    triangle_area=.5*((x1-ax)*(y2-ay)-(y1-ay)*(x2-ax))
    if np.min(triangle_area)<-1e-14:raise ValueError('Inconsistent face orientation')
    ref_area=np.add.reduceat(triangle_area,starts)
    ref_cx=np.add.reduceat(triangle_area*(ax+x1+x2)/3.,starts)/ref_area
    ref_cy=np.add.reduceat(triangle_area*(ay+y1+y2)/3.,starts)/ref_area
    ux=ax-ref_cx[owner];uy=ay-ref_cy[owner]
    vx=x1-ref_cx[owner];vy=y1-ref_cy[owner]
    wx=x2-ref_cx[owner];wy=y2-ref_cy[owner]
    q=(ux*ux+uy*uy+vx*vx+vy*vy+wx*wx+wy*wy+
       ux*vx+uy*vy+vx*wx+vy*wy+wx*ux+wy*uy)
    ref_inertia=mass*np.add.reduceat(triangle_area*q/6.,starts)/ref_area
    area_rel=np.abs(area_c/ref_area-1.)
    centroid_abs=np.hypot(cx_c-ref_cx,cy_c-ref_cy)
    inertia_rel=np.abs(inertia_c/ref_inertia-1.)
    native_area_rel=np.abs(area_c/p['volume'][active].astype(float)-1.)
    independent_area_rel=np.abs(area_c/area-1.)
    checks={'all_moments_positive':bool(np.all(area_c>0) and np.all(inertia_c>0)),
            'triangle_area_agreement':bool(np.max(area_rel)<1e-10),
            'triangle_centroid_agreement':bool(np.max(centroid_abs)<1e-10),
            'triangle_inertia_agreement':bool(np.max(inertia_rel)<1e-10),
            'independent_power_area_agreement':bool(np.max(independent_area_rel)<1e-10),
            'native_saved_area_agreement':bool(np.max(native_area_rel)<1e-5)}
    result={'receipt':str(args.receipt.resolve()),'receipt_sha256':digest(args.receipt),
            'endpoint_sha256':receipt['endpoint_sha256'],'time':float(receipt['time']),
            'native_failed_checks':sorted(failed),'weighted':bool(np.any(weights!=0)),
            'disk_sites':n,'directed_disk_faces':len(owner),'min_faces':int(counts.min()),
            'max_faces':int(counts.max()),'area_relative_error_max':float(area_rel.max()),
            'centroid_absolute_error_max':float(centroid_abs.max()),
            'inertia_relative_error_max':float(inertia_rel.max()),
            'inertia_relative_error_p99':float(np.quantile(inertia_rel,.99)),
            'independent_power_area_error_max':float(independent_area_rel.max()),
            'native_saved_area_error_max':float(native_area_rel.max()),
            'checks':checks,'scope':'CPU execution of CUDA-compatible polygon formula on '
            'independently reconstructed native checkpoint faces; triangle integration is '
            'a separate moment reference. Native GPU face extraction and GPU cost remain separate.'}
    with args.output.open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'time':result['time'],'weighted':result['weighted'],
                      'faces':len(owner),'checks':checks,
                      'inertia_relative_error_max':result['inertia_relative_error_max']},indent=2))
    if not all(checks.values()):raise SystemExit('Native polygon moment comparison failed')


if __name__=='__main__':main()
