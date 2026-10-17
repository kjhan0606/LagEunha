"""Vectorized forms of the existing triangle integrals and planar cut areas."""
import numpy as np


def triangle_moments_many(parts):
    v=np.asarray(parts);a=v[:,1]-v[:,0];b=v[:,2]-v[:,0]
    area=abs(a[:,0]*b[:,1]-a[:,1]*b[:,0])/2;s=v.sum(1)
    Q=area[:,None,None]*(np.einsum('nvi,nvj->nij',v,v)+s[:,:,None]*s[:,None,:])/12
    return np.column_stack([area,area[:,None]*s/3,Q.reshape(-1,4)])


def normalized_many(raw,center,scale):
    m=raw[:,0];q=raw[:,1:3];Q=raw[:,3:7].reshape(-1,2,2)
    S=(Q-center[None,:,None]*q[:,None,:]-q[:,:,None]*center[None,None,:]+m[:,None,None]*np.outer(center,center))/scale**2
    return np.column_stack([m,(q-m[:,None]*center)/scale,S[:,0,0],S[:,0,1],S[:,1,1]])


def phase_averages(parts,metadata,low,high):
    d=np.asarray(parts)@np.asarray(metadata['normal'])-metadata['offset'];fraction=np.zeros(len(d))
    nonnegative=d.min(1)>=0;fraction[nonnegative]=1.
    mixed=(d.max(1)>0)&(d.min(1)<0);positive=(d>0).sum(1)
    one=mixed&(positive==1)
    p=d[one].max(1);denominator=np.prod(np.where(d[one]>0,1.,p[:,None]-d[one]),axis=1)
    fraction[one]=p*p/denominator
    two=mixed&(positive==2)
    n=d[two].min(1);denominator=np.prod(np.where(d[two]<0,1.,n[:,None]-d[two]),axis=1)
    fraction[two]=1-n*n/denominator
    return low+(high-low)*np.clip(fraction,0.,1.)
