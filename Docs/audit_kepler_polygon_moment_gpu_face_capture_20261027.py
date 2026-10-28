"""Check polygon moments on captured native CUDA face endpoints."""

import argparse
import csv
import ctypes
import json
from pathlib import Path

import numpy as np

from prepare_kepler_sparse_area_restart_20261027 import digest


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('case',choices=('aligned','jitter'))
    ap.add_argument('--library',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    docs=Path(__file__).resolve().parent
    reference=docs/'KEPLER_Q2_BOUNDED_FACEFORCE_PAIR_T6501_2026-10-27.json'
    row=json.loads(reference.read_text())['cases'][args.case]
    root=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/faceforce_pair_t6501_a10_v1')/args.case
    faces=root/'face_dump_r0_call1.tsv';inputs=root/'face_input_r0_call1.tsv'
    if digest(faces)!=row['face_dump_sha256'] or digest(inputs)!=row['face_input_sha256']:
        raise ValueError('Native GPU capture changed')
    with inputs.open() as f:
        inputs_by_id={int(v['id']):v for v in csv.DictReader(f,delimiter='\t')}
    grouped={}
    with faces.open() as f:
        for v in csv.DictReader(f,delimiter='\t'):
            if int(v['wall'])!=0:raise ValueError('Unexpected wall face in disk capture')
            grouped.setdefault(int(v['id_i']),[]).append(v)
    ids=sorted(grouped)
    if len(ids)!=64 or any(i not in inputs_by_id for i in ids):
        raise ValueError('Expected 64 selected GPU cells and matching inputs')
    counts=np.array([len(grouped[i]) for i in ids],dtype=np.int32)
    if np.any(counts<3):raise ValueError('Incomplete captured polygon')
    offset=np.r_[0,np.cumsum(counts)].astype(np.int32)
    a=[];b=[];mass=[];native_volume=[]
    for i in ids:
        source=inputs_by_id[i];x=float(source['x']);y=float(source['y'])
        mass.append(float(source['mass']));native_volume.append(float(source['V']))
        for v in grouped[i]:
            a.append((float(v['x1'])-x,float(v['y1'])-y))
            b.append((float(v['x2'])-x,float(v['y2'])-y))
    a=np.asarray(a);b=np.asarray(b);n=len(ids)
    x1=np.ascontiguousarray(a[:,0]);y1=np.ascontiguousarray(a[:,1])
    x2=np.ascontiguousarray(b[:,0]);y2=np.ascontiguousarray(b[:,1])
    mass=np.ascontiguousarray(mass,dtype=float)
    area=np.empty(n);cx=np.empty(n);cy=np.empty(n);inertia=np.empty(n)
    status=np.empty(n,dtype=np.int32)
    pd=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    pi=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
    lib=ctypes.CDLL(str(args.library.resolve()));fn=lib.gfs_polygon_moment_batch
    fn.argtypes=[ctypes.c_int,pi,pd,pd,pd,pd,pd,pd,pd,pd,pd,pi];fn.restype=ctypes.c_int
    code=fn(n,offset,x1,y1,x2,y2,mass,area,cx,cy,inertia,status)
    if code!=0 or np.any(status!=0):raise ValueError(f'Invalid captured cells: {np.count_nonzero(status)}')
    starts=offset[:-1];owner=np.repeat(np.arange(n),counts)
    ax=np.add.reduceat(x1,starts)/counts;ay=np.add.reduceat(y1,starts)/counts
    tx=ax[owner];ty=ay[owner]
    tri=.5*((x1-tx)*(y2-ty)-(y1-ty)*(x2-tx))
    if np.min(tri)<-1e-14:raise ValueError('Captured face orientation is inconsistent')
    ref_area=np.add.reduceat(tri,starts)
    ref_cx=np.add.reduceat(tri*(tx+x1+x2)/3.,starts)/ref_area
    ref_cy=np.add.reduceat(tri*(ty+y1+y2)/3.,starts)/ref_area
    ux=tx-ref_cx[owner];uy=ty-ref_cy[owner]
    vx=x1-ref_cx[owner];vy=y1-ref_cy[owner]
    wx=x2-ref_cx[owner];wy=y2-ref_cy[owner]
    q=(ux*ux+uy*uy+vx*vx+vy*vy+wx*wx+wy*wy+
       ux*vx+uy*vy+vx*wx+vy*wy+wx*ux+wy*uy)
    ref_inertia=mass*np.add.reduceat(tri*q/6.,starts)/ref_area
    area_error=float(np.max(np.abs(area/ref_area-1.)))
    centroid_error=float(np.max(np.hypot(cx-ref_cx,cy-ref_cy)))
    inertia_error=float(np.max(np.abs(inertia/ref_inertia-1.)))
    native_area_error=float(np.max(np.abs(area/np.asarray(native_volume)-1.)))
    checks={'area_vs_triangle':area_error<1e-10,'centroid_vs_triangle':centroid_error<1e-10,
            'inertia_vs_triangle':inertia_error<1e-10,'area_vs_native':native_area_error<1e-5}
    result={'case':args.case,'time':6.5,'cells':n,'directed_faces':len(a),
            'face_dump_sha256':digest(faces),'face_input_sha256':digest(inputs),
            'reference_sha256':digest(reference),'area_relative_error_max':area_error,
            'centroid_absolute_error_max':centroid_error,'inertia_relative_error_max':inertia_error,
            'area_vs_native_relative_error_max':native_area_error,'checks':checks,
            'scope':'CPU execution of CUDA-compatible moment routine on 64 captured native '
                    'GPU cell-face polygons; not device execution of the moment kernel.'}
    with args.output.open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(result,indent=2,allow_nan=False))
    if not all(checks.values()):raise SystemExit('Native GPU face-capture moment audit failed')


if __name__=='__main__':main()
