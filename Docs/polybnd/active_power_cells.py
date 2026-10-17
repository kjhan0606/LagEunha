"""Experimental all-constraint clipping with vectorized violation selection.

Every site remains a constraint. Select a violated plane, clip, and repeat;
discard no neighbor on a geometric heuristic. Applied planes are checked again
at the end to bound accumulated floating-point clipping error.
"""
import numpy as np
from shapely.geometry import Polygon
from weighted_mass_patch import clip_halfplane


def power_cells_active(sites,weights,domain,seed_count=0):
    sites=np.asarray(sites,dtype=float);weights=np.asarray(weights,dtype=float);domain=np.asarray(domain,dtype=float)
    if sites.ndim!=2 or sites.shape[1]!=2 or weights.shape!=(len(sites),):raise ValueError('invalid sites/weights')
    if not all(np.all(np.isfinite(v)) for v in (sites,weights,domain)):raise ValueError('finite geometry required')
    cells=[];scale=max(1.,float(np.ptp(domain,axis=0).max()));seeds=None
    if seed_count:
        if seed_count<0:raise ValueError('nonnegative seed count required')
        from scipy.spatial import cKDTree
        _,seeds=cKDTree(sites).query(sites,k=min(len(sites),int(seed_count)+1))
        seeds=np.asarray(seeds).reshape(len(sites),-1)
    for i,site in enumerate(sites):
        delta=sites-site;norm=2*np.linalg.norm(delta,axis=1)
        if np.count_nonzero(norm==0)>1:raise ValueError('distinct sites required')
        normals=2*delta;bounds=np.einsum('ij,ij->i',delta,delta)+weights[i]-weights
        active=np.zeros(len(sites),dtype=bool);active[i]=True
        p=domain-site
        # Seed only the initial intersection. Every unseeded plane remains
        # in the global violation search and final all-plane certificate.
        if seeds is not None:
            for j in seeds[i]:
                if active[j]:continue
                if np.max(p@normals[j]-bounds[j])>0:p=clip_halfplane(p,normals[j],bounds[j])
                active[j]=True
                if len(p)<3:break
        if len(p)<3:cells.append(Polygon());continue
        for _ in range(len(sites)-1):
            residual=(normals@p.T-bounds[:,None]).max(axis=1)
            distance=np.divide(residual,norm,out=np.full(len(sites),-np.inf),where=norm>0)
            distance[active]=-np.inf;j=int(np.argmax(distance))
            if distance[j]<=0:break
            p=clip_halfplane(p,normals[j],bounds[j]);active[j]=True
            if len(p)<3:break
        if len(p)>=3:
            residual=(normals@p.T-bounds[:,None]).max(axis=1)
            tolerance=128*np.finfo(float).eps*(abs(normals)@np.max(abs(p),axis=0)+abs(bounds)+scale**2)
            if np.any(residual>tolerance):raise RuntimeError('active clipping failed full constraint verification')
            cells.append(Polygon(p+site))
        else:cells.append(Polygon())
    return cells
