"""Positive density with spatial second moments and one potential moment.
Orthogonalize the potential against quadratics to avoid a redundant feature.
Research fit only; transport and realizability over time require separate tests.
"""
import numpy as np

class PotentialDensity:
    def __init__(self,poly,m,q,Q,potential_moment,potential,quadrature,order=24):
        pts,w=quadrature(poly,order);self.area=w.sum();self.origin=w@pts/self.area
        self.scale=np.sqrt(self.area);self.mass=float(m);self.potential=potential
        self.points=pts;self.weights=w
        F=self.polynomial(pts);phi=potential(pts)
        self.beta=np.linalg.lstsq(np.sqrt(w)[:,None]*F,np.sqrt(w)*phi,rcond=None)[0]
        residual=phi-F@self.beta
        self.potential_scale=float(np.sqrt(w@(residual*residual)/self.area))
        if self.potential_scale<1e-13:raise ValueError('potential feature numerically redundant')
        A=np.column_stack([F[:,1:],residual/self.potential_scale])
        ld=np.longdouble;o=self.origin.astype(ld);c=np.asarray(q,dtype=ld)/ld(m)
        second=np.asarray(Q,dtype=ld)/ld(m)-np.outer(o,c)-np.outer(c,o)+np.outer(o,o)
        target_poly=np.r_[ld(1),(c-o)/ld(self.scale),second[0,0]/ld(self.area),second[0,1]/ld(self.area),second[1,1]/ld(self.area)]
        target_extra=(ld(potential_moment)/ld(m)-self.beta.astype(ld)@target_poly)/ld(self.potential_scale)
        target=np.asarray(np.r_[target_poly[1:],target_extra],dtype=float);theta=np.zeros(6)
        def evaluate(t):
            ex=A@t;shift=ex.max();weighted=w*np.exp(ex-shift);Z=weighted.sum();pr=weighted/Z
            mean=pr@A;center=A-mean;H=center.T@(pr[:,None]*center)
            return mean-target,H,shift+np.log(Z),pr
        for it in range(60):
            r,H,logZ,pr=evaluate(theta);norm=float(np.max(abs(r)))
            if norm<1e-12:break
            step=np.linalg.solve(H,r)
            for bt in range(30):
                trial=theta-step*2.**(-bt)
                if np.max(abs(evaluate(trial)[0]))<norm:theta=trial;break
            else:raise RuntimeError('potential density fit stalled')
        else:raise RuntimeError('potential density fit did not converge')
        self.theta=theta;self.logZ=logZ;self.residual=norm;self.iterations=it+1
        self.condition=float(np.linalg.cond(H));self.weighted_mass=self.mass*pr

    def polynomial(self,pts):
        z=(np.atleast_2d(pts)-self.origin)/self.scale
        return np.column_stack([np.ones(len(z)),z,z[:,0]**2,z[:,0]*z[:,1],z[:,1]**2])

    def value(self,pts):
        pts=np.atleast_2d(pts);F=self.polynomial(pts)
        A=np.column_stack([F[:,1:],(self.potential(pts)-F@self.beta)/self.potential_scale])
        return self.mass*np.exp(A@self.theta-self.logZ)
