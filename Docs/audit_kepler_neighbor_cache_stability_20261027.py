"""Measure how often real-cell neighbors persist across saved native steps."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest

N=55351


def directed_pairs(snapshot,ids):
    xy=np.column_stack((snapshot['x'].astype(float)-3.,snapshot['y'].astype(float)-3.))
    lookup=np.full(N,-1,dtype=np.int64);lookup[ids]=np.arange(len(ids))
    area,meta,triangles,_=power_areas(xy,snapshot['w2'].astype(float),ids,lookup,return_geometry=True)
    if meta['hidden_active'] or np.any(area<=0):raise ValueError('Invalid power cells')
    edges=np.concatenate((triangles[:,(0,1)],triangles[:,(1,2)],triangles[:,(2,0)]))
    edges.sort(axis=1);edges=np.unique(edges,axis=0)
    active=np.zeros(N,dtype=bool);active[ids]=True
    left=edges[active[edges[:,0]]];right=edges[active[edges[:,1]]]
    owners=np.concatenate((left[:,0],right[:,1]))
    neigh=np.concatenate((left[:,1],right[:,0]))
    keys=owners.astype(np.int64)*N+neigh
    return np.unique(keys)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('run',type=Path)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    receipt=json.loads((args.run/'audit.json').read_text())
    if not all(receipt['checks'].values()):raise ValueError('Native replay gate failed')
    before=Path(receipt.get('restart',receipt.get('parent')));after=Path(receipt['endpoint'])
    expected=receipt.get('restart_sha256',receipt.get('parent_sha256'))
    if digest(before)!=expected or digest(after)!=receipt['endpoint_sha256']:
        raise ValueError('Changed checkpoint')
    p0,m0=ordered(args.run/'initial.dat',N)
    a,ma=ordered(before,N);b,mb=ordered(after,N)
    initial_xy=np.column_stack((p0['x'].astype(float)-3.,p0['y'].astype(float)-3.))
    r0=np.linalg.norm(initial_xy,axis=1);ids=np.flatnonzero((r0>=.5)&(r0<2.))
    old=directed_pairs(a,ids);new=directed_pairs(b,ids)
    common=np.intersect1d(old,new,assume_unique=True)
    removed=np.setdiff1d(old,new,assume_unique=True);added=np.setdiff1d(new,old,assume_unique=True)
    changed_ids=np.unique(np.concatenate((removed//N,added//N)))
    new_count=np.bincount(new//N,minlength=N);common_count=np.bincount(common//N,minlength=N)
    if np.any(new_count[ids]==0):raise ValueError('Missing disk neighbor sets')
    bands=[]
    for lo,hi in ((.5,.75),(.75,1.),(1.,1.5),(1.5,2.)):
        members=ids[(r0[ids]>=lo)&(r0[ids]<hi)]
        changed=np.isin(members,changed_ids,assume_unique=True)
        bands.append({'r_min':lo,'r_max':hi,'cells':len(members),
                      'unchanged_neighbor_set_fraction':float(1.-np.mean(changed)),
                      'new_neighbor_coverage':float(common_count[members].sum()/new_count[members].sum())})
    result={'run':str(args.run.resolve()),'receipt_sha256':digest(args.run/'audit.json'),
            'initial_sha256':digest(args.run/'initial.dat'),'before_sha256':digest(before),
            'after_sha256':digest(after),'before_time':float(ma['t']),'after_time':float(mb['t']),
            'weighted':bool(np.any(a['w2']!=0) or np.any(b['w2']!=0)),
            'disk_sites':len(ids),'before_directed_neighbors':len(old),
            'after_directed_neighbors':len(new),'retained_directed_neighbors':len(common),
            'added_directed_neighbors':len(added),'removed_directed_neighbors':len(removed),
            'unchanged_neighbor_set_fraction':float(1.-len(changed_ids)/len(ids)),
            'new_neighbor_coverage_by_old':float(len(common)/len(new)),
            'radial_bands':bands,'checks':{'native_gate':True,'valid_finite_power_diagrams':True},
            'scope':'Independently reconstructed regular-triangulation adjacency at two '
                    'native snapshots; not a CUDA cache implementation or proof that cached '
                    'indices, MPI halo images or early preclipping are safe.'}
    with args.output.open('x') as out:json.dump(result,out,indent=2,allow_nan=False);out.write('\n')
    print(json.dumps({'weighted':result['weighted'],'unchanged':result['unchanged_neighbor_set_fraction'],
                      'coverage':result['new_neighbor_coverage_by_old'],'added':len(added)},indent=2))


if __name__=='__main__':main()
