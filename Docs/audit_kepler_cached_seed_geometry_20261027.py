"""Compare old-neighbor and directional seed geometry on native snapshots."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest

N=55351
CELL=.0625
XMIN=-.0625
MX=98


def clip(poly,dx,dy,threshold):
    result=[]
    for i in range(len(poly)):
        a=poly[i];b=poly[(i+1)%len(poly)]
        fa=dx*a[0]+dy*a[1]-threshold;fb=dx*b[0]+dy*b[1]-threshold
        ain=fa<=0.;bin_=fb<=0.
        if ain:result.append(a)
        if ain!=bin_:
            t=fa/(fa-fb);result.append(a+t*(b-a))
    return np.asarray(result,dtype=float).reshape((-1,2))


def initial_polygon():
    h=CELL
    return np.array([[-h,-h],[h,-h],[h,h],[-h,h]],dtype=float)


def clip_indices(poly,indices,xy,weights,i):
    for j in indices:
        delta=xy[j]-xy[i];d2=float(delta@delta)
        if d2<=0:continue
        frac=np.clip(.5+.5*(weights[i]-weights[j])/d2,.05,.95)
        poly=clip(poly,delta[0],delta[1],frac*d2)
        if len(poly)<3:raise ValueError('Seed clipping produced an empty cell')
    return poly


def main():
    ap=argparse.ArgumentParser();ap.add_argument('run',type=Path)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--samples',type=int,default=1024)
    args=ap.parse_args()
    # Analytic right-half-plane intersection guards the prototype clipper.
    test=clip(initial_polygon(),1.,0.,0.)
    if len(test)!=4 or abs(np.max(test[:,0]))>1e-15 or np.min(test[:,0])!=-CELL:
        raise ValueError('Analytic clipping check failed')
    receipt=json.loads((args.run/'audit.json').read_text())
    if not all(receipt['checks'].values()):raise ValueError('Native parent gate failed')
    old_path=Path(receipt.get('restart',receipt.get('parent')));new_path=Path(receipt['endpoint'])
    old_sha=receipt.get('restart_sha256',receipt.get('parent_sha256'))
    if digest(old_path)!=old_sha or digest(new_path)!=receipt['endpoint_sha256']:
        raise ValueError('Checkpoint changed')
    p0,_=ordered(args.run/'initial.dat',N);old,oldmeta=ordered(old_path,N);new,newmeta=ordered(new_path,N)
    xy0=np.column_stack((p0['x'].astype(float)-3.,p0['y'].astype(float)-3.))
    r0=np.linalg.norm(xy0,axis=1);disk=np.flatnonzero((r0>=.5)&(r0<2.))
    lookup=np.full(N,-1,dtype=int);lookup[disk]=np.arange(len(disk))
    oldxy=np.column_stack((old['x'].astype(float)-3.,old['y'].astype(float)-3.))
    _,meta,triangles,_=power_areas(oldxy,old['w2'].astype(float),disk,lookup,return_geometry=True)
    if meta['hidden_active']:raise ValueError('Hidden old cell')
    edges=np.concatenate((triangles[:,(0,1)],triangles[:,(1,2)],triangles[:,(2,0)]))
    edges.sort(axis=1);edges=np.unique(edges,axis=0)
    left=edges[lookup[edges[:,0]]>=0];right=edges[lookup[edges[:,1]]>=0]
    owners=np.concatenate((left[:,0],right[:,1]));neigh=np.concatenate((left[:,1],right[:,0]))
    order=np.argsort(owners,kind='stable');owners=owners[order];neigh=neigh[order]
    offsets=np.r_[0,np.cumsum(np.bincount(owners,minlength=N))]
    xy=np.column_stack((new['x'].astype(float),new['y'].astype(float)))
    weights=new['w2'].astype(float)
    ij=((xy-XMIN)/CELL).astype(int)
    if np.any(ij<0) or np.any(ij>=MX):raise ValueError('Out-of-range linked-cell index')
    bins={}
    for k,(ix,iy) in enumerate(ij):bins.setdefault((int(ix),int(iy)),[]).append(k)
    rng=np.random.default_rng(20261027)
    sample=np.sort(rng.choice(disk,size=min(args.samples,len(disk)),replace=False))
    fallback=0;near_counts=[];cache_counts=[];old_radius=[];cache_radius=[];old_pass=[];cache_pass=[]
    for i in sample:
        ix,iy=ij[i]
        nearby=np.asarray([j for cy in range(iy-2,iy+3) for cx in range(ix-2,ix+3)
                           for j in bins.get((cx,cy),()) if j!=i],dtype=int)
        delta=xy[nearby]-xy[i];d2=np.sum(delta*delta,axis=1);valid=d2>0
        nearby,delta,d2=nearby[valid],delta[valid],d2[valid]
        frac=np.clip(.5+.5*(weights[i]-weights[nearby])/d2,.05,.95)
        dr2=frac*frac*d2
        sectors=np.where(np.abs(delta[:,0])>=np.abs(delta[:,1]),
             np.where(delta[:,0]>=0,np.where(delta[:,1]>=0,0,7),np.where(delta[:,1]>=0,3,4)),
             np.where(delta[:,1]>=0,np.where(delta[:,0]>=0,1,2),np.where(delta[:,0]>=0,6,5)))
        anchors=[int(nearby[np.flatnonzero(sectors==s)[np.argmin(dr2[sectors==s])]])
                 for s in range(8) if np.any(sectors==s)]
        baseline=clip_indices(initial_polygon(),anchors,xy,weights,i)
        previous=neigh[offsets[i]:offsets[i+1]]
        cached=clip_indices(initial_polygon(),previous,xy,weights,i)
        closed=not np.any(np.abs(cached)>=CELL*(1-32*np.finfo(float).eps))
        if not closed:fallback+=1;cached=baseline
        rb=float(np.max(np.sum(baseline*baseline,axis=1)))
        rc=float(np.max(np.sum(cached*cached,axis=1)))
        near_counts.append(len(anchors));cache_counts.append(len(previous));old_radius.append(rb);cache_radius.append(rc)
        old_pass.append(int(np.sum(dr2<=rb)));cache_pass.append(int(np.sum(dr2<=rc)))
    old_radius=np.asarray(old_radius);cache_radius=np.asarray(cache_radius)
    result={'run':str(args.run.resolve()),'receipt_sha256':digest(args.run/'audit.json'),
            'before_sha256':old_sha,'after_sha256':receipt['endpoint_sha256'],
            'before_time':float(oldmeta['t']),'after_time':float(newmeta['t']),
            'weighted':bool(np.any(old['w2']!=0) or np.any(new['w2']!=0)),
            'sample_cells':len(sample),'sampling_seed':20261027,'fallback_cells':fallback,
            'directional_seed_planes_mean':float(np.mean(near_counts)),
            'cached_seed_planes_mean':float(np.mean(cache_counts)),
            'cached_vs_directional_radius2_median':float(np.median(cache_radius/old_radius)),
            'cached_vs_directional_radius2_p95':float(np.quantile(cache_radius/old_radius,.95)),
            'seed_candidates_within_initial_radius_mean_directional':float(np.mean(old_pass)),
            'seed_candidates_within_initial_radius_mean_cached':float(np.mean(cache_pass)),
            'checks':{'native_gate':True,'analytic_clip':True,'finite_power_diagram':True},
            'scope':'Sampled CPU geometric prototype only. Candidate counts are upper bounds '
                    'before further streaming clips; the real GPU scans and capacity fallback '
                    'are not implemented or timed. Not proof of bitwise or MPI equivalence.'}
    with args.output.open('x') as out:json.dump(result,out,indent=2,allow_nan=False);out.write('\n')
    print(json.dumps({k:result[k] for k in ('weighted','sample_cells','fallback_cells',
          'cached_vs_directional_radius2_median','seed_candidates_within_initial_radius_mean_directional',
          'seed_candidates_within_initial_radius_mean_cached')},indent=2))


if __name__=='__main__':main()
