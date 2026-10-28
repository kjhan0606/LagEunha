"""Outward intervals for encoded polygon coordinates and cell mean density.

Default coordinate errors are a working affine-roundoff envelope, not a proof
of all errors in the upstream clipping geometry. Supplied errors may include
known geometric uncertainty.
"""
import numpy as np

def outward(lo,hi):
    return np.nextafter(lo,-np.inf),np.nextafter(hi,np.inf)

def add(a,b):return outward(a[0]+b[0],a[1]+b[1])
def subtract(a,b):return outward(a[0]-b[1],a[1]-b[0])
def multiply(a,b):
    corners=np.array([a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]])
    return outward(corners.min(axis=0),corners.max(axis=0))
def divide_positive(a,b):
    if np.any(np.asarray(b[0])<=0):raise ValueError('strictly positive interval denominator required')
    reciprocal=outward(1/b[1],1/b[0])
    return multiply(a,reciprocal)

def working_vertex_error(poly,area):
    u=np.finfo(float).eps/2;gamma4=4*u/(1-4*u)
    return gamma4*(abs(np.asarray(poly))+np.sqrt(area))

def affine_vertex_error(poly,error,rotation,shift):
    u=np.finfo(float).eps/2;gamma4=4*u/(1-4*u)
    return np.asarray(error)@abs(np.asarray(rotation)).T+gamma4*(abs(np.asarray(poly))@abs(np.asarray(rotation)).T+abs(np.asarray(shift)))

def area_interval(poly,error):
    p=np.asarray(poly);e=np.asarray(error)
    if e.shape!=p.shape or not np.isfinite(e).all() or np.any(e<0):raise ValueError('invalid vertex error bounds')
    lo,hi=outward(p-e,p+e)
    r=subtract((lo,hi),(lo[0],hi[0]));total=(0.,0.)
    for i in range(len(p)):
        j=(i+1)%len(p)
        positive=multiply((r[0][i,0],r[1][i,0]),(r[0][j,1],r[1][j,1]))
        negative=multiply((r[0][i,1],r[1][i,1]),(r[0][j,0],r[1][j,0]))
        total=add(total,subtract(positive,negative))
    result=multiply(total,(.5,.5))
    if not np.isfinite(result).all() or result[0]<=0:
        raise ValueError('polygon orientation/area unresolved within coordinate uncertainty')
    return tuple(float(v) for v in result)

def mean_density_interval(poly,area,mass,vertex_error=None,mass_error=0.):
    if not np.isfinite(mass_error) or mass_error<0:raise ValueError('invalid mass error bound')
    e=working_vertex_error(poly,area) if vertex_error is None else np.asarray(vertex_error)
    area_bounds=area_interval(poly,e);mass_bounds=outward(mass-mass_error,mass+mass_error)
    if mass_bounds[0]<=0:raise ValueError('positive mass not resolved')
    rho=divide_positive(mass_bounds,area_bounds)
    return area_bounds,tuple(float(v) for v in rho),e
