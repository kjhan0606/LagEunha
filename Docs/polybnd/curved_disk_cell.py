"""Convex rational power polygon intersected with a true disk.

Supporting lines/linear clips are rational. Circle intersections carry rational
outward enclosures of their algebraic roots; tiny positive pieces are retained.
Unresolved angular ordering raises instead of replacing arcs by chords.
"""
from dataclasses import dataclass
from decimal import Decimal,localcontext
from fractions import Fraction as F
from functools import lru_cache
from math import atan2,isqrt,pi
import numpy as np

class UnresolvedCurvedTopology(RuntimeError):pass

@lru_cache(None)
def gauss(order):return np.polynomial.legendre.leggauss(order)
def fraction(x):return x if isinstance(x,F) else (F(x) if isinstance(x,str) else F(float(x)))
def ld(x):return np.longdouble(x.numerator)/np.longdouble(x.denominator)
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def sub(a,b):return tuple(x-y for x,y in zip(a,b))

def sqrt_enclosure(q):
    if q<0:raise ValueError('negative radicand')
    a=isqrt(q.numerator);b=isqrt(q.denominator)
    if a*a==q.numerator and b*b==q.denominator:return F(a,b),F(a,b)
    with localcontext() as ctx:
        ctx.prec=90
        s=(Decimal(q.numerator)/Decimal(q.denominator)).sqrt()
        lo=F(s);hi=F(s)
        while lo*lo>q:lo=F(ctx.next_minus(Decimal(lo.numerator)/Decimal(lo.denominator)))
        while hi*hi<q:hi=F(ctx.next_plus(Decimal(hi.numerator)/Decimal(hi.denominator)))
    assert lo*lo<=q<=hi*hi
    return lo,hi

def circle_roots(a,b,radius):
    """Enclosures of segment-parameter roots; exact rational discriminant."""
    d=sub(b,a);A=dot(d,d)
    if A==0:raise UnresolvedCurvedTopology('zero exact edge')
    B=2*dot(a,d);C=dot(a,a)-radius*radius;D=B*B-4*A*C
    if D<0:return []
    low,high=sqrt_enclosure(D)
    return [((-B-high)/(2*A),(-B-low)/(2*A)),((-B+low)/(2*A),(-B+high)/(2*A))]

def midpoint(interval):return (interval[0]+interval[1])/2

@dataclass
class StraightBoundary:
    neighbor:int
    A:tuple
    B:tuple
    root_enclosures:tuple=()
    def quadrature(self,order=16):
        d=sub(self.B,self.A);length2=dot(d,d)
        if length2<=0:raise UnresolvedCurvedTopology('nonpositive straight measure')
        lo,hi=sqrt_enclosure(length2);length=ld((lo+hi)/2)
        delta=np.array([ld(v) for v in d]);start=np.array([ld(v) for v in self.A]);g,w=gauss(order)
        p=start+(g.astype(np.longdouble)+1)[:,None]*delta/2
        normal=np.array([delta[1],-delta[0]])/length
        return p,np.broadcast_to(normal,p.shape).copy(),length*w.astype(np.longdouble)/2

@dataclass
class CircularBoundary:
    radius:float
    lo:float
    hi:float
    def quadrature(self,order=48):
        g,w=gauss(order);theta=np.longdouble(self.lo)+(g.astype(np.longdouble)+1)*np.longdouble(self.hi-self.lo)/2
        n=np.column_stack([np.cos(theta),np.sin(theta)])
        return np.longdouble(self.radius)*n,n,np.longdouble(self.radius)*np.longdouble(self.hi-self.lo)*w.astype(np.longdouble)/2

class DiskCell:
    def __init__(self,vertices,tags,radius,owner=None):
        self.vertices=[tuple(map(fraction,p)) for p in vertices];self.tags=list(map(int,tags));self.radius=float(radius);self.owner=owner
        if len(self.vertices)<3 or len(self.tags)!=len(self.vertices) or radius<=0:raise UnresolvedCurvedTopology('invalid enclosing polygon')
        self.R=fraction(radius);self.R2=self.R*self.R
        self.planes=[]
        for a,b in zip(self.vertices,self.vertices[1:]+self.vertices[:1]):
            d=sub(b,a)
            if dot(d,d)==0:raise UnresolvedCurvedTopology('zero exact edge')
            normal=(d[1],-d[0]);self.planes.append((normal,dot(normal,a)))
        self.lines=[];knots=[0.,2*pi];crossings=[]
        inside=[dot(p,p)<=self.R2 for p in self.vertices]
        knots.extend(atan2(float(p[1]),float(p[0]))%(2*pi) for p in self.vertices if dot(p,p)==self.R2)
        for k,(a,b) in enumerate(zip(self.vertices,self.vertices[1:]+self.vertices[:1])):
            if inside[k] and inside[(k+1)%len(inside)]:
                self.lines.append(StraightBoundary(self.tags[k],a,b));continue
            roots=circle_roots(a,b,self.R)
            if not roots:continue
            lower=(max(F(0),roots[0][0]),max(F(0),roots[0][1]));upper=(min(F(1),roots[1][0]),min(F(1),roots[1][1]))
            if upper[1]<=lower[0]:continue
            if upper[0]<=lower[1]:raise UnresolvedCurvedTopology('circle interval ordering unresolved')
            d=sub(b,a);lo=midpoint(lower);hi=midpoint(upper)
            aa=tuple(a[j]+lo*d[j] for j in (0,1));bb=tuple(a[j]+hi*d[j] for j in (0,1))
            self.lines.append(StraightBoundary(self.tags[k],aa,bb,(lower,upper)))
            for root in roots:
                if root[1]<0 or root[0]>1:continue
                if root[0]<0 or root[1]>1:
                    # An exact endpoint root is rational and has zero width.
                    raise UnresolvedCurvedTopology('circle root at segment endpoint unresolved')
                t=midpoint(root);point=tuple(a[j]+t*d[j] for j in (0,1));angle=atan2(float(point[1]),float(point[0]))%(2*pi)
                crossings.append((angle,point));knots.append(angle)
        ordered=sorted(crossings,key=lambda a:a[0])
        for (theta,a),(theta2,b) in zip(ordered,ordered[1:]):
            if theta==theta2 and a!=b:
                # Adjacent edges meeting at the same algebraic endpoint can
                # have tiny enclosure differences. Exact polygon endpoint
                # coincidences are allowed; distinct unresolved angles fail.
                if not any(dot(v,v)==self.R2 and all(abs(float(v[j]-a[j]))==0 and abs(float(v[j]-b[j]))==0 for j in (0,1)) for v in self.vertices):
                    raise UnresolvedCurvedTopology('distinct circular crossings collapse in angle encoding')
        self.arcs=[];knots=sorted(set(knots))
        relevant=[(n,c) for n,c in self.planes if not(c>=0 and c*c>=self.R2*dot(n,n))]
        for lo,hi in zip(knots[:-1],knots[1:]):
            theta=np.longdouble(lo)+np.longdouble(hi-lo)/2
            p=np.longdouble(self.radius)*np.array([np.cos(theta),np.sin(theta)])
            allowed=True;unresolved=False
            for n,c in relevant:
                normal=np.array([ld(v) for v in n]);value=normal@p-ld(c)
                bound=128*np.finfo(np.longdouble).eps*(abs(normal)@abs(p)+abs(ld(c)))+8*np.finfo(float).eps*abs(normal).sum()*self.radius
                if abs(value)<=bound:unresolved=True
                elif value>0:allowed=False;break
            if allowed and unresolved:raise UnresolvedCurvedTopology('arc halfplane classification unresolved')
            if allowed:self.arcs.append(CircularBoundary(self.radius,lo,hi))
        if not self.lines and not self.arcs:raise UnresolvedCurvedTopology('empty disk intersection')
        self._moments=None

    @classmethod
    def from_native(cls,cell,radius):
        return cls([edge['A'] for edge in cell['edges']],[edge['neighbor'] for edge in cell['edges']],radius,cell['owner'])

    def boundary_quadrature(self,line_order=16,arc_order=48):
        values=[line.quadrature(line_order) for line in self.lines]+[arc.quadrature(arc_order) for arc in self.arcs]
        return tuple(np.concatenate([v[k] for v in values]) for k in range(3))

    def moments(self,origin=None):
        p,n,w=self.boundary_quadrature();o=np.zeros(2,dtype=np.longdouble) if origin is None else np.asarray(origin,dtype=np.longdouble)
        x=p-o;radial=w*np.sum(x*n,axis=1);area=radial.sum()/2;q=radial@x/3;Q=np.einsum('n,ni,nj->ij',radial,x,x)/4
        if area<=0 or not np.isfinite(area):raise UnresolvedCurvedTopology('nonpositive disk area')
        return area,q,Q

    @property
    def area(self):
        if self._moments is None:self._moments=self.moments()
        return float(self._moments[0])
    @property
    def origin(self):
        if self._moments is None:self._moments=self.moments()
        return np.asarray(self._moments[1]/self._moments[0],dtype=float)

    def clip_halfplane(self,normal,offset,tag=-1000):
        normal=tuple(map(fraction,normal));offset=fraction(offset);out=[];tags=[]
        for k,(a,b) in enumerate(zip(self.vertices,self.vertices[1:]+self.vertices[:1])):
            sa=dot(normal,a)-offset;sb=dot(normal,b)-offset;ai=sa<=0;bi=sb<=0
            if ai:out.append(a);tags.append(self.tags[k])
            if ai!=bi:
                t=sa/(sa-sb);d=sub(b,a);out.append(tuple(a[j]+t*d[j] for j in (0,1)));tags.append(tag if ai else self.tags[k])
        clean=[];ct=[]
        for p,t in zip(out,tags):
            if clean and p==clean[-1]:ct[-1]=t
            else:clean.append(p);ct.append(t)
        if len(clean)>1 and clean[-1]==clean[0]:clean.pop();ct.pop()
        if len(clean)<3:return None
        return DiskCell(clean,ct,self.radius,self.owner)

    def positive_quadrature(self,order=24,radii=()):
        """Positive polar slices; all radial density/potential breaks explicit.

        Angular envelope changes occur at polygon vertices or intersections
        of a supporting segment and a radial interface. No wall chord enters.
        """
        knots=[0.,2*pi]
        for v in self.vertices:
            if v!=(F(0),F(0)):knots.append(atan2(float(v[1]),float(v[0]))%(2*pi))
        radii=tuple(sorted(set(float(r) for r in radii if 0<r<self.radius)))
        for a,b in zip(self.vertices,self.vertices[1:]+self.vertices[:1]):
            d=sub(b,a)
            for radius in radii+(self.radius,):
                for interval in circle_roots(a,b,fraction(radius)):
                    if interval[0]<0 or interval[1]>1:continue
                    t=midpoint(interval);point=tuple(a[j]+t*d[j] for j in (0,1))
                    knots.append(atan2(float(point[1]),float(point[0]))%(2*pi))
        for normal,c in self.planes:
            if c==0:
                angle=atan2(float(normal[1]),float(normal[0]))
                knots.extend([(angle+pi/2)%(2*pi),(angle-pi/2)%(2*pi)])
        knots=sorted(set(knots));g,w=gauss(order);points=[];weights=[]
        normals=np.array([[ld(v) for v in n] for n,c in self.planes]);offsets=np.array([ld(c) for n,c in self.planes])
        # A line close to the origin can make b/(n dot er) vary sharply
        # even between all exact envelope/radial intersection breakpoints.
        # Bound its denominator variation geometrically (at most factor2 per
        # panel) before applying Gauss. This is independent of the integrand
        # or fit residual and introduces no omitted-area cutoff.
        def conditioned(lo,hi):
            mid=np.longdouble(lo)+np.longdouble(hi-lo)/2
            unit=np.array([np.cos(mid),np.sin(mid)]);den=normals@unit
            low=np.longdouble(0);high=np.longdouble(self.radius);active=[];ilow=ihigh=None
            for k,(d,b) in enumerate(zip(den,offsets)):
                if d>0 and b/d<high:high=b/d;ihigh=k
                elif d<0 and b/d>low:low=b/d;ilow=k
                elif d==0 and b<0:return True
            if high<=low:return True
            active=[k for k in (ilow,ihigh) if k is not None and offsets[k]!=0]
            for k in active:
                n=normals[k]
                d0=n@np.array([np.cos(np.longdouble(lo)),np.sin(np.longdouble(lo))])
                d1=n@np.array([np.cos(np.longdouble(hi)),np.sin(np.longdouble(hi))])
                if d0*d1<=0 or max(abs(d0),abs(d1))>2*min(abs(d0),abs(d1)):return False
            return True
        panels=[];pending=[(lo,hi,0) for lo,hi in zip(knots[:-1],knots[1:])]
        while pending:
            lo,hi,depth=pending.pop()
            if conditioned(lo,hi):panels.append((lo,hi));continue
            mid=lo+(hi-lo)/2
            if depth>=64 or mid==lo or mid==hi:raise UnresolvedCurvedTopology('angular quadrature conditioning unresolved')
            pending.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
        for lo,hi in sorted(panels):
            theta=np.longdouble(lo)+(g.astype(np.longdouble)+1)*np.longdouble(hi-lo)/2
            direction=np.column_stack([np.cos(theta),np.sin(theta)]);den=direction@normals.T
            for direction0,den0,wtheta in zip(direction,den,w*np.longdouble(hi-lo)/2):
                lower=np.longdouble(0);upper=np.longdouble(self.radius);valid=True
                for value,offset in zip(den0,offsets):
                    if value>0:upper=min(upper,offset/value)
                    elif value<0:lower=max(lower,offset/value)
                    elif offset<0:valid=False;break
                if not valid or upper<=lower:continue
                cuts=[lower]+[np.longdouble(r) for r in radii if lower<r<upper]+[upper]
                for a,b in zip(cuts[:-1],cuts[1:]):
                    rr=a+(g.astype(np.longdouble)+1)*(b-a)/2
                    points.append(rr[:,None]*direction0);weights.append(wtheta*w.astype(np.longdouble)*(b-a)*rr/2)
        if not points:raise UnresolvedCurvedTopology('empty positive polar quadrature')
        points=np.asarray(np.concatenate(points),dtype=float);weights=np.asarray(np.concatenate(weights),dtype=float)
        if not np.isfinite(points).all() or not np.isfinite(weights).all() or not np.all(weights>0):raise UnresolvedCurvedTopology('unresolved positive quadrature')
        return points,weights
