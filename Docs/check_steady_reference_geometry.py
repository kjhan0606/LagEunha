"""Check original face endpoints used by the exact-volume reference."""
import json,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');sys.path.insert(0,str(base))
from run_canonical_av_branches_20261018 import read
h=np.load(base/'Codex_continuous_initial_rhs_quadratic_split_20261018/rhs_r0.55_f8.npz')
ic,t=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');assert t==0
allx=np.column_stack([ic['x']-3,ic['y']-3])
faces=np.loadtxt(base/'Codex_pair_hg_modes_t0_20261018/raw_faces.dat.rank0000')
order=np.argsort(faces[:,0],kind='stable');owners=faces[order,0]
vertex_error=0.;segment_error=0.;normal_error=0.;relative_area_error=0.
for i in h['ids']:
    rows=faces[order[np.searchsorted(owners,i):np.searchsorted(owners,i,side='right')]]
    v=rows[:,6:8];poly=v[np.argsort(np.arctan2(v[:,1],v[:,0]))]
    ends=rows[:,8:10]
    vertex_error=max(vertex_error,float(np.linalg.norm(ends[:,None,:]-poly[None,:,:],axis=2).min(axis=1).max()))
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        mismatch=np.minimum(np.linalg.norm(rows[:,6:8]-a,axis=1)+np.linalg.norm(rows[:,8:10]-b,axis=1),
                            np.linalg.norm(rows[:,6:8]-b,axis=1)+np.linalg.norm(rows[:,8:10]-a,axis=1))
        segment_error=max(segment_error,float(mismatch.min()))
    lengths=np.linalg.norm(rows[:,8:10]-rows[:,6:8],axis=1)
    delta=allx[rows[:,1].astype(int)]-allx[i]
    normal=delta/np.linalg.norm(delta,axis=1)[:,None]
    normal_error=max(normal_error,float(np.linalg.norm(lengths@normal)))
    area=.5*np.sum(poly[:,0]*np.roll(poly[:,1],-1)-poly[:,1]*np.roll(poly[:,0],-1))
    assert area>0
    relative_area_error=max(relative_area_error,abs(area/ic['volume'][i]-1))
assert vertex_error<1e-10 and segment_error<1e-10 and normal_error<1e-9 and relative_area_error<1e-9
print(json.dumps(dict(passed=True,cells=len(h['ids']),maximum_second_endpoint_to_polygon_distance=vertex_error,
 maximum_polygon_edge_to_face_distance=segment_error,maximum_closed_boundary_normal_sum=normal_error,
 maximum_relative_area_difference=relative_area_error,scope='Original cached t0 geometry only; not moving-geometry backend validation.'),indent=2))
