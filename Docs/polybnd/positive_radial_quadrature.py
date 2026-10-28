"""Positive polygon quadrature split explicitly at radial spline interfaces."""
from functools import lru_cache
import numpy as np

@lru_cache(None)
def gauss(order):
    x,w=np.polynomial.legendre.leggauss(int(order));return x,w

def rule(poly,order=24,radii=()):
    p=np.asarray(poly,dtype=float);d=np.roll(p,-1,axis=0)-p
    if p.ndim!=2 or p.shape[1]!=2 or len(p)<3:raise ValueError('polygon required')
    knots=list(p[:,0]);xmin=float(p[:,0].min());xmax=float(p[:,0].max())
    for radius in radii:
        knots.extend(x for x in (-radius,radius) if xmin<x<xmax)
        for a,v in zip(p,d):
            A=float(v@v)
            if A==0:continue
            B=2*float(a@v);C=float(a@a)-radius**2;disc=B*B-4*A*C
            if disc<0:continue
            q=-.5*(B+np.copysign(np.sqrt(disc),B))
            roots=[-B/(2*A)] if q==0 else [q/A,C/q]
            knots.extend(float((a+t*v)[0]) for t in roots if 0<t<1)
    knots=sorted(set(knots));g,w=gauss(order);points=[];weights=[]
    for xl,xr in zip(knots[:-1],knots[1:]):
        if xr<=xl:continue
        xx=(xl+xr)/2+(xr-xl)*g/2
        for x,wx in zip(xx,(xr-xl)*w/2):
            nonvertical=d[:,0]!=0;t=(x-p[nonvertical,0])/d[nonvertical,0]
            hit=(t>=0)&(t<=1);ys=p[nonvertical,1][hit]+t[hit]*d[nonvertical,1][hit]
            if len(ys)<2:raise ValueError('unresolved positive vertical polygon slice')
            yl=float(ys.min());yh=float(ys.max());cuts=[yl,yh]
            for radius in radii:
                if abs(x)<radius:
                    height=np.sqrt(radius**2-x*x)
                    cuts.extend(y for y in (-height,height) if yl<y<yh)
            cuts=sorted(set(cuts))
            for lo,hi in zip(cuts[:-1],cuts[1:]):
                if hi<=lo:continue
                yy=(lo+hi)/2+(hi-lo)*g/2
                points.append(np.column_stack([np.full(len(g),x),yy]));weights.append(wx*(hi-lo)*w/2)
    points=np.concatenate(points);weights=np.concatenate(weights)
    if not np.isfinite(points).all() or not np.all(weights>0):raise ValueError('invalid positive quadrature')
    return points,weights
