"""Explicit pressure fit alternatives; both retain the original cell mean.

The tensor/free-intercept fit includes the target cell and uses separable
compact weights. After fitting, its intercept is discarded and the existing
zero-mean polynomial basis retains the prescribed thermal-energy mean.
"""
import os
import numpy as np


def mode():
    value=os.environ.get('LAGFORCE_PRESSURE_FIT','radial_fixed_mean')
    if value not in ('radial_fixed_mean','tensor_free_mean'):raise ValueError('unknown pressure fit')
    return value


def tensor_stencil(centers,center,scale,support):
    distance=abs((centers-center)/(support*scale));ids=np.flatnonzero(np.max(distance,axis=1)<1.)
    r=distance[ids];weight=np.prod((1-r)**4*(1+4*r),axis=1)
    return ids,weight


def free_mean_fit(A,target,weight,ridge=1e-14):
    matrix=np.sqrt(weight)[:,None]*np.column_stack([np.ones(len(A)),A])
    augmented=np.vstack([matrix,np.sqrt(ridge)*np.eye(A.shape[1]+1)])
    result=np.linalg.lstsq(augmented,np.r_[np.sqrt(weight)*target,np.zeros(A.shape[1]+1)],rcond=None)[0]
    return result[1:]
