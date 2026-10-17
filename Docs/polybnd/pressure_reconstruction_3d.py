"""Cell-mean preserving quadratic pressure on 3D convex cells.

A finite-volume least-squares fit uses geometric neighbor moments. Tetrahedral
quadratic Bernstein coefficients certify positivity after a mean-preserving
scaling. This is an explicit spatial model, never a fallback.
"""
import numpy as np
from scipy.spatial import ConvexHull
from moment_density_3d import PAIRS


def features(points,center,scale):
    z=(np.atleast_2d(points)-center)/scale
    return np.column_stack([z]+[z[:,i]*z[:,j] for i,j in PAIRS])


def value(metadata,points):
    return metadata['mean']+(features(points,metadata['center'],metadata['scale'])-metadata['feature_mean'])@metadata['coefficient']


def prepare(states,polys):
    centers=np.array([s['density'].origin for s in states]);volumes=np.array([s['density'].volume for s in states]);P=np.array([s['P'] for s in states]);cov=[]
    for s,c,V in zip(states,centers,volumes):
        d=s['density'];dx=d.points-c;cov.append(dx.T@(d.weights[:,None]*dx)/V)
    cov=np.array(cov);result=[];limited=0;min_theta=1.
    for i,(s,poly) in enumerate(zip(states,polys)):
        h=volumes[i]**(1/3);z=(centers-centers[i])/h;radius=np.linalg.norm(z,axis=1)/3
        ids=np.flatnonzero((radius<1)&(np.arange(len(states))!=i));second=cov/h**2+z[:,:,None]*z[:,None,:]
        expected=np.column_stack([z]+[second[:,j,k] for j,k in PAIRS]);mean=expected[i];design=expected[ids]-mean
        r=radius[ids];w=(1-r)**4*(1+4*r)
        coefficient,_,rank,_=np.linalg.lstsq(design*np.sqrt(w)[:,None],(P[ids]-P[i])*np.sqrt(w),rcond=1e-12)
        if rank<9:raise RuntimeError('3D quadratic pressure stencil deficient')
        metadata=dict(mean=P[i],center=centers[i],scale=h,feature_mean=mean,coefficient=coefficient)
        vertex=np.vstack([centers[i],poly]);pv=value(metadata,vertex);lower=float(pv.min());hull=ConvexHull(poly)
        for face in hull.simplices:
            tetra=np.r_[0,face+1]
            for a in range(4):
                for b in range(a+1,4):
                    u,v=tetra[a],tetra[b];bernstein=2*value(metadata,(vertex[u]+vertex[v])/2)[0]-.5*(pv[u]+pv[v]);lower=min(lower,float(bernstein))
        floor=1e-12*P[i];theta=1.
        if lower<floor:
            theta=(P[i]-floor)/(P[i]-lower);metadata['coefficient']=coefficient*theta;limited+=1
        metadata['positivity_scale']=theta;min_theta=min(min_theta,theta);result.append(metadata)
    return result,dict(pressure_limited_cells=limited,pressure_min_positivity_scale=min_theta)
