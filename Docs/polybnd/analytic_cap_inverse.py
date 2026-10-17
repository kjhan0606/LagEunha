"""Experimental convex-polygon cap inversion; not enabled in evolution.

Triangulate the polygon once. Each triangle's cap area is a piecewise quadratic
function of offset; solve the polynomial on the appropriate vertex interval.
The returned polygon still uses the existing clipping function.
"""
import numpy as np
from shapely.geometry import Polygon
from discontinuous_density_audit import clip


def cap_at_fraction(vertices,normal,fraction):
    vertices=np.asarray(vertices,dtype=float);normal=np.asarray(normal,dtype=float)
    if not np.isfinite(fraction) or not 0<=fraction<=1 or not np.all(np.isfinite(normal)) or np.linalg.norm(normal)==0:
        raise ValueError('invalid cap parameters')
    area=Polygon(vertices).area
    if area<=0:raise ValueError('positive convex polygon required')
    heights=vertices@normal
    if fraction==0:return float(heights.max()),Polygon()
    if fraction==1:return float(heights.min()),Polygon(vertices)
    if fraction>.5:
        offset,_=cap_at_fraction(vertices,-normal,1-fraction)
        return -offset,clip(vertices,normal,-offset)
    triangles=np.array([[vertices[0],vertices[i],vertices[i+1]] for i in range(1,len(vertices)-1)])
    u=triangles[:,1]-triangles[:,0];v=triangles[:,2]-triangles[:,0]
    weights=abs(u[:,0]*v[:,1]-u[:,1]*v[:,0])/2
    keep=weights>0;weights=weights[keep];triangles=triangles[keep]
    weights*=area/weights.sum()
    z=np.sort(triangles@normal,axis=1);z0,z1,z2=z.T

    def cap(offset):
        values=np.zeros(len(z));values[offset<=z0]=weights[offset<=z0]
        lower=(offset>z0)&(offset<z1)
        values[lower]=weights[lower]*(1-(offset-z0[lower])**2/((z1[lower]-z0[lower])*(z2[lower]-z0[lower])))
        upper=(offset>=z1)&(offset<z2)
        values[upper]=weights[upper]*(z2[upper]-offset)**2/((z2[upper]-z0[upper])*(z2[upper]-z1[upper]))
        return float(values.sum())

    knots=np.unique(heights);target=fraction*area;lo=hi=None
    for left,right in zip(knots[:-1],knots[1:]):
        if cap(left)>=target and cap(right)<=target:lo,hi=float(left),float(right);break
    if lo is None:raise RuntimeError('cap interval unavailable')
    anchor=hi;middle=lo+(hi-lo)/2;aa=bb=cc=0.
    for weight,(low,mid,high) in zip(weights,z):
        if middle<=low:cc+=weight
        elif middle>=high:continue
        elif middle<mid:
            scale=weight/((mid-low)*(high-low));delta=anchor-low
            aa-=scale;bb-=2*scale*delta;cc+=weight-scale*delta**2
        else:
            scale=weight/((high-low)*(high-mid));delta=high-anchor
            aa+=scale;bb-=2*scale*delta;cc+=scale*delta**2
    cc-=target;width=hi-lo;roots=[]
    if abs(aa)*width<=64*np.finfo(float).eps*abs(bb):
        if bb!=0:roots=[-cc/bb]
    else:
        discriminant=bb*bb-4*aa*cc
        if discriminant>=0:
            q=-.5*(bb+np.copysign(np.sqrt(discriminant),bb))
            if q!=0:roots=[q/aa,cc/q]
    tolerance=64*np.finfo(float).eps*max(1.,abs(lo),abs(hi))
    candidates=[float(np.clip(anchor+root,lo,hi)) for root in roots if -width-tolerance<=root<=tolerance]
    if candidates:
        offset=min(candidates,key=lambda c:abs(cap(c)-target))
    else:
        # Retain a safe monotone fallback for nearly degenerate coefficients.
        for _ in range(52):
            midpoint=lo+(hi-lo)/2
            if cap(midpoint)>target:lo=midpoint
            else:hi=midpoint
        offset=lo+(hi-lo)/2
    return offset,clip(vertices,normal,offset)
