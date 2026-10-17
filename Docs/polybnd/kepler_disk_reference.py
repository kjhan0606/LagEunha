"""Canonical non-vacuum Kepler disk fields and positive split polar quadrature.

Definition: box[0,6]^2, center(3,3), rho1 on .5<=r<=2, rotating background
rho.01, uniform pressure, spline rc=.5, outer frozen radius2.5. Reference fields
are for initialization and scoring, never for repairing evolved states.
"""
import numpy as np
from scipy.spatial import cKDTree,ConvexHull
from kepler_spline_gravity import enclosed_mass,potential

GAMMA=5/3;CENTER=np.array([3.,3.]);RC=.5
PRESSURE=(GAMMA-1)*(5*GAMMA/2)*1e-5


def frozen(t,points):return np.linalg.norm(np.asarray(points)-CENTER,axis=1)>=2.5


def fields(points):
    x=np.asarray(points)-CENTER;r=np.linalg.norm(x,axis=1);rho=np.where((r>=.5)&(r<=2.),1.,.01)
    omega=np.full_like(r,np.sqrt(32/(3*RC**3)));positive=r>0
    omega[positive]=np.sqrt(enclosed_mass(r[positive]/RC)/r[positive]**3)
    velocity=omega[:,None]*np.column_stack([-x[:,1],x[:,0]])
    return rho,velocity,np.full_like(r,PRESSURE)


def polar_rule(poly,order=16,radii=(.25,.5,2.)):
    p=np.asarray(poly,dtype=float)-CENTER
    # Multiple exact bisectors can repeat a vertex; remove only identical
    # consecutive coordinates, preserving the polygonal region exactly.
    p=p[np.any(p!=np.roll(p,1,axis=0),axis=1)]
    if len(p)<3:raise ValueError('degenerate polygon')
    local=p-p[0]
    area=.5*np.sum(local[:,0]*np.roll(local[:,1],-1)-local[:,1]*np.roll(local[:,0],-1))
    if area==0:raise ValueError('zero-area polygon')
    if area<0:p=p[::-1]
    # Power cells are convex. Roundoff-sized inward vertex loops can create
    # spurious halfspaces that exclude a macroscopic wedge at the polar origin.
    # A convex hull removes those loops only if its area agrees to roundoff.
    hull=ConvexHull(p);convex=p[hull.vertices]
    diameter=float(np.max(np.ptp(p,axis=0)))
    if abs(hull.volume-abs(area))>128*np.finfo(float).eps*diameter**2:
        raise ValueError('nonconvex polygon exceeds roundoff allowance')
    p=convex
    edge=np.roll(p,-1,axis=0)-p;n=np.column_stack([edge[:,1],-edge[:,0]]);b=np.sum(n*p,axis=1)
    lengths=np.linalg.norm(n,axis=1)
    if np.any(lengths<=0):raise ValueError('degenerate polygon edge')
    angles=[-np.pi,np.pi,*np.arctan2(p[:,1],p[:,0]).tolist()]
    radius_extent=float(np.max(np.linalg.norm(p,axis=1)))
    for normal,bound,length in zip(n,b,lengths):
        direction=np.arctan2(normal[1],normal[0])
        # A face close to the polar origin creates a narrow angular layer in
        # r=b/(n.direction). Resolve it geometrically around both denominator
        # zeros; increasing one global Gauss order misses this layer.
        angular_scale=abs(bound)/(length*radius_extent)
        if 0.<angular_scale<1e-4:
            offsets=angular_scale*2.**np.arange(0,int(np.ceil(np.log2(1/angular_scale)))+1)
            for pole in (direction-np.pi/2,direction+np.pi/2):
                for offset in offsets:
                    angles.extend([((pole-offset+np.pi)%(2*np.pi))-np.pi,((pole+offset+np.pi)%(2*np.pi))-np.pi])
        for radius in radii:
            ratio=bound/(radius*length)
            if abs(ratio)<1:
                delta=np.arccos(ratio)
                angles.extend([((direction-delta+np.pi)%(2*np.pi))-np.pi,((direction+delta+np.pi)%(2*np.pi))-np.pi])
    angles=np.unique(angles);z,w=np.polynomial.legendre.leggauss(order);points=[];weights=[]
    for left,right in zip(angles[:-1],angles[1:]):
        if right-left<1e-14:continue
        # Every vertex angle is a breakpoint, so ray visibility is constant
        # inside this interval for a convex polygon. Reject empty intervals
        # once instead of evaluating every Gauss node outside the cell.
        middle=(left+right)/2;mid_a=n@np.array([np.cos(middle),np.sin(middle)])
        mid_positive=mid_a>1e-15*lengths;mid_negative=mid_a<-1e-15*lengths
        if np.any((~mid_positive&~mid_negative)&(b<0)):continue
        mid_upper=np.min(b[mid_positive]/mid_a[mid_positive],initial=np.inf)
        mid_lower=max(0.,float(np.max(b[mid_negative]/mid_a[mid_negative],initial=-np.inf)))
        if mid_upper<=mid_lower:continue
        theta=(left+right)/2+(right-left)/2*z;tw=(right-left)/2*w
        for angle,angle_weight in zip(theta,tw):
            direction=np.array([np.cos(angle),np.sin(angle)]);a=n@direction;positive=a>1e-15*lengths;negative=a<-1e-15*lengths
            if np.any((~positive&~negative)&(b<0)):continue
            upper=np.min(b[positive]/a[positive],initial=np.inf);lower=max(0.,float(np.max(b[negative]/a[negative],initial=-np.inf)))
            if upper<=lower:continue
            if not np.isfinite(upper):raise ValueError('unbounded polar ray')
            cuts=[lower]+[r for r in radii if lower<r<upper]+[upper]
            for lo,hi in zip(cuts[:-1],cuts[1:]):
                r=(lo+hi)/2+(hi-lo)/2*z;rw=(hi-lo)/2*w
                points.append(CENTER+r[:,None]*direction);weights.append(angle_weight*rw*r)
    if not points:raise ValueError('empty polar quadrature')
    weights=np.concatenate(weights)
    if abs(weights.sum()-hull.volume)>2e-10*hull.volume:
        raise ValueError('polar quadrature failed independent polygon-area check')
    return np.vstack(points),weights


def physical_raw(poly,order=16):
    x,w=polar_rule(poly,order);rho,v,P=fields(x);dm=w*rho;m=dm.sum();q=dm@x;Q=x.T@(dm[:,None]*x);p=dm@v;D=x.T@(dm[:,None]*v);Pi=v.T@(dm[:,None]*v)
    c=q/m;S=Q-np.outer(q,q)/m;C=D-np.outer(q,p)/m;B=np.linalg.solve(S,C).T;residual=v-p/m-(x-c)@B.T
    _,factor=np.linalg.qr(np.sqrt(dm)[:,None]*residual,mode='reduced')
    E=w@P/(GAMMA-1)+.5*np.trace(Pi)
    raw=np.r_[m,q,Q.ravel(),p,D.ravel(),Pi.ravel(),E]
    return raw,factor.T,dict(volume=float(w.sum()),potential=float(dm@potential(x,center=CENTER,rc=RC)))


def ring_sites(n,require_count=True):
    """Exact layout algorithm from ICgraded/mkrings_uniform.py; no file writes."""
    L=6.;spacing=L/n;xs=[];ys=[];k=1
    while k*spacing<3*np.sqrt(2)+spacing:
        r=k*spacing;dd=spacing if r<2.7 else .985*spacing;count=int(round(2*np.pi*r/dd));theta=-np.pi+(np.arange(count)+(.6180339887*k)%1.)*2*np.pi/count
        xs.append(3+r*np.cos(theta));ys.append(3+r*np.sin(theta));k+=1
    x=np.concatenate(xs);y=np.concatenate(ys);inside=(x>0)&(x<L)&(y>0)&(y<L);x=x[inside];y=y[inside]
    points=np.column_stack([x,y]);pairs=cKDTree(points,boxsize=[L,L]).query_pairs(.7*spacing);drop=set(b for a,b in pairs);keep=np.ones(len(points),dtype=bool);keep[list(drop)]=False;points=points[keep]
    need=n*n-len(points)
    if need>0 and require_count:raise ValueError('canonical layout requires padding at this resolution: '+str(need))
    if need<0:
        radius=np.hypot(points[:,0]-3.,points[:,1]-3.);order=np.argsort(-radius);keep=np.ones(len(points),dtype=bool);keep[order[:-need]]=False;points=points[keep]
    return points
