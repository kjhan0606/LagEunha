"""Small all-constraint power-cell/disk reference with true circular arcs.

Static geometry contract, not the production embedded polygon builder.
No generator projection, radius cutoff, or arc-to-chord replacement.
Closed fixed-mass solves must remove the constant-weight gauge explicitly.
"""
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class Line:
    neighbor:int
    normal:np.ndarray
    offset:float
    lo:float
    hi:float
    def quadrature(self,order=16):
        g,w=np.polynomial.legendre.leggauss(order)
        t=(self.lo+self.hi)/2+(self.hi-self.lo)/2*g
        tangent=np.array([-self.normal[1],self.normal[0]])
        x=self.offset*self.normal+t[:,None]*tangent
        return x,np.broadcast_to(self.normal,x.shape).copy(),w*(self.hi-self.lo)/2

@dataclass(frozen=True)
class Arc:
    radius:float
    lo:float
    hi:float
    def quadrature(self,order=48):
        g,w=np.polynomial.legendre.leggauss(order)
        theta=(self.lo+self.hi)/2+(self.hi-self.lo)/2*g
        normal=np.column_stack([np.cos(theta),np.sin(theta)])
        return self.radius*normal,normal,w*self.radius*(self.hi-self.lo)/2

@dataclass
class CircularCell:
    lines:list
    arcs:list
    radius:float
    def quadrature(self):
        parts=[piece.quadrature() for piece in self.lines+self.arcs]
        if not parts:raise ValueError('empty or unresolved disk-clipped cell')
        return tuple(np.concatenate([p[k] for p in parts]) for k in range(3))
    def moments(self):
        """Green formula for homogeneous volume moments through degree2."""
        x,n,w=self.quadrature();radial=np.sum(x*n,axis=1)*w
        area=np.sum(radial)/2
        q=radial@x/3
        Q=np.einsum('n,ni,nj->ij',radial,x,x)/4
        if not np.isfinite(area) or area<=0:raise ValueError('empty or unresolved disk-clipped cell')
        return area,q,Q


def cells(sites,weights,radius):
    sites=np.asarray(sites,dtype=float);weights=np.asarray(weights,dtype=float)
    if sites.ndim!=2 or sites.shape[1]!=2 or weights.shape!=(len(sites),):raise ValueError('invalid sites/weights shape')
    if not np.isfinite(sites).all() or not np.isfinite(weights).all() or not np.isfinite(radius) or radius<=0:raise ValueError('finite sites/weights and positive radius required')
    out=[]
    for i,site in enumerate(sites):
        delta=sites-site;distance=np.linalg.norm(delta,axis=1);ids=np.flatnonzero(np.arange(len(sites))!=i)
        if np.any(distance[ids]==0):raise ValueError('coincident generators unresolved')
        normals=delta[ids]/distance[ids,None]
        offsets=(np.sum(sites[ids]**2,axis=1)-site@site+weights[i]-weights[ids])/(2*distance[ids])
        lines=[];knots=[0.,2*np.pi]
        for k,(normal,b) in enumerate(zip(normals,offsets)):
            if b<=-radius:raise ValueError('empty or tangent-only disk-clipped cell')
            if b>=radius:continue
            angle=np.arctan2(normal[1],normal[0]);opening=np.arccos(b/radius)
            knots.extend([(angle-opening)%(2*np.pi),(angle+opening)%(2*np.pi)])
            tangent=np.array([-normal[1],normal[0]])
            length=np.sqrt((radius-b)*(radius+b));lo=-length;hi=length
            for l,(other,c) in enumerate(zip(normals,offsets)):
                if l==k:continue
                slope=other@tangent;remaining=c-b*(other@normal)
                if slope==0:
                    if remaining<0:hi=lo;break
                elif slope>0:hi=min(hi,remaining/slope)
                else:lo=max(lo,remaining/slope)
            if hi>lo:lines.append(Line(int(ids[k]),normal,float(b),float(lo),float(hi)))
        arcs=[];knots=sorted(set(knots))
        for lo,hi in zip(knots[:-1],knots[1:]):
            mid=(lo+hi)/2;x=radius*np.array([np.cos(mid),np.sin(mid)])
            if np.all(normals@x<=offsets):arcs.append(Arc(radius,lo,hi))
        cell=CircularCell(lines,arcs,radius);cell.moments();out.append(cell)
    return out


def mean_zero_gauge(matrix,residual):
    """Contract for future closed massroot: solve J z=-r with sum(z)=0.

    Requires a compatible residual. No frozen anchor and no ridge/pinning.
    This dense reference is for small static fixtures only.
    """
    matrix=np.asarray(matrix,dtype=float);residual=np.asarray(residual,dtype=float);n=len(residual)
    scale=max(float(np.linalg.norm(matrix,np.inf)),1.)
    tol=64*n*np.finfo(float).eps*scale
    if np.max(abs(matrix@np.ones(n)))>tol:raise ValueError('missing constant-weight nullspace')
    if abs(residual.sum())>64*n*np.finfo(float).eps*max(float(abs(residual).sum()),1.):raise ValueError('incompatible closed mass residual')
    augmented=np.zeros((n+1,n+1));augmented[:n,:n]=matrix;augmented[:n,n]=1;augmented[n,:n]=1
    return np.linalg.solve(augmented,np.r_[-residual,0.])[:n]
