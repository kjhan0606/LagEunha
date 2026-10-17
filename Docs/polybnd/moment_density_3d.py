"""Positive exponential-quadratic density on convex 3D cells.

Three first and six independent second moments plus mass are fitted using
positive tetrahedral quadrature. This is a component of the moment model;
no 3D fluid evolution is supplied by this module.
"""
import numpy as np
from scipy.spatial import ConvexHull

PAIRS=((0,0),(0,1),(0,2),(1,1),(1,2),(2,2))


def quadrature(vertices,order=6):
    vertices=np.asarray(vertices,dtype=float)
    if vertices.ndim!=2 or vertices.shape[1]!=3 or len(vertices)<4 or not np.all(np.isfinite(vertices)):
        raise ValueError('finite 3D convex cell vertices required')
    if order<2:raise ValueError('quadrature order must be at least two')
    hull=ConvexHull(vertices);center=vertices[hull.vertices].mean(axis=0)
    z,w=np.polynomial.legendre.leggauss(order);z=(z+1)/2;w=w/2
    u,v,t=np.meshgrid(z,z,z,indexing='ij');wu,wv,wt=np.meshgrid(w,w,w,indexing='ij')
    u=u.ravel();v=v.ravel();t=t.ravel()
    bary=np.column_stack([1-u,u*(1-v),u*v*(1-t),u*v*t])
    weights=(wu*wv*wt).ravel()*u*u*v
    points=[];masses=[]
    for face in hull.simplices:
        tetra=np.vstack([center,vertices[face]])
        jac=abs(np.linalg.det(tetra[1:]-tetra[0]))
        if jac<=0:raise ValueError('degenerate tetrahedron')
        points.append(bary@tetra);masses.append(jac*weights)
    return np.vstack(points),np.concatenate(masses)


def features(points,origin,scale):
    r=(np.atleast_2d(points)-origin)/scale
    return np.column_stack([r]+[r[:,i]*r[:,j] for i,j in PAIRS])


class Density3D:
    def __init__(self,vertices,m,q,Q,order=6,tolerance=3e-13):
        q=np.asarray(q,dtype=float);Q=np.asarray(Q,dtype=float)
        if q.shape!=(3,) or Q.shape!=(3,3) or not np.isfinite(m) or m<=0 or not np.all(np.isfinite(q)) or not np.all(np.isfinite(Q)):
            raise ValueError('finite positive mass and 3D moments required')
        if not np.allclose(Q,Q.T,rtol=0,atol=1e-13*max(1.,np.max(abs(Q)))):
            raise ValueError('symmetric second moment required')
        points,weights=quadrature(vertices,order)
        self.mass=float(m);self.volume=float(weights.sum());self.origin=weights@points/self.volume
        self.scale=self.volume**(1/3);self.points=points;self.weights=weights
        c=q/m
        if np.linalg.eigvalsh((Q+Q.T)/(2*m)-np.outer(c,c)).min()<=0:
            raise ValueError('nonpositive position covariance')
        second=Q/m-np.outer(self.origin,c)-np.outer(c,self.origin)+np.outer(self.origin,self.origin)
        target=np.r_[(c-self.origin)/self.scale,[second[i,j]/self.scale**2 for i,j in PAIRS]]
        F=features(points,self.origin,self.scale);theta=np.zeros(9)
        def evaluate(theta):
            exponent=F@theta;shift=exponent.max();a=weights*np.exp(exponent-shift);Z=a.sum();prob=a/Z
            mean=prob@F;centered=F-mean
            return mean-target,(centered*prob[:,None]).T@centered,shift+np.log(Z),prob
        for iteration in range(60):
            residual,H,logZ,prob=evaluate(theta);norm=float(max(abs(residual)))
            if norm<tolerance:break
            try:step=np.linalg.solve(H,residual)
            except np.linalg.LinAlgError as e:raise RuntimeError('3D density moments are unresolved or infeasible') from e
            for bt in range(30):
                candidate=theta-step*2.**(-bt)
                if max(abs(evaluate(candidate)[0]))<norm:theta=candidate;break
            else:raise RuntimeError('3D positive density fit stalled; no moment clipping applied')
        else:raise RuntimeError('3D positive density fit did not converge')
        self.theta=theta;self.logZ=logZ;self.residual=norm;self.iterations=iteration+1
        self.weighted_mass=m*prob

    def value(self,points):
        return self.mass*np.exp(features(points,self.origin,self.scale)@self.theta-self.logZ)
