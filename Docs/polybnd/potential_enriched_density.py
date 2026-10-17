"""Positive density matching original m,q,Q and a transported potential moment.

This representation is opt-in. The caller supplies physical potential and
positive volume quadrature; no energy or moment target is changed on failure.
"""
import numpy as np
from moment_ale_evolve import Density,features
from positive_exponential_tilt import fit


class PotentialDensity:
    def __init__(self,poly,m,q,Q,potential_moment,potential,quadrature):
        self.base=Density(poly,m,q,Q);self.mass=float(m);self.potential=potential
        self.origin=self.base.origin;self.scale=self.base.scale
        self.points,self.weights=quadrature(poly);self.area=float(self.weights.sum())
        if not np.isfinite(potential_moment):raise ValueError('finite potential moment required')
        b=self._basis(self.points);dm=self.weights*self.base.value(self.points)
        values=np.asarray(potential(self.points))
        if values.shape!=(len(dm),) or not np.isfinite(values).all():raise ValueError('invalid potential')
        A=b.T@((dm/self.mass)[:,None]*b)
        self.projection=np.linalg.solve(A,b.T@(dm*values/self.mass))
        remainder=values-b@self.projection
        self.sigma=float(np.sqrt(dm@(remainder**2)/self.mass))
        if self.sigma<=64*np.finfo(float).eps*max(1.,float(np.max(abs(values)))):
            raise ValueError('potential basis is numerically dependent on retained moments')
        # Center using extended arithmetic, preserving supplied moment targets.
        o=np.asarray(self.origin,dtype=np.longdouble);q=np.asarray(q,dtype=np.longdouble)
        Q=np.asarray(Q,dtype=np.longdouble);c=q/self.mass
        second=Q/self.mass-np.outer(o,c)-np.outer(c,o)+np.outer(o,o)
        target=np.asarray(np.r_[(c-o)/self.scale,second[0,0]/self.scale**2,
            second[0,1]/self.scale**2,second[1,1]/self.scale**2],dtype=float)
        target_b=np.r_[1.,target]
        target=np.r_[target,(potential_moment/self.mass-self.projection@target_b)/self.sigma]
        self.theta,self.logZ,p,diag=fit(dm/self.mass,self._features(self.points),target)
        self.weighted_mass=self.mass*p;self.residual=diag['residual'];self.iterations=diag['iterations']
        self.potential_moment=float(potential_moment)

    def _basis(self,x):
        x=np.atleast_2d(x)
        return np.column_stack([np.ones(len(x)),features(x,self.origin,self.scale)])

    def _features(self,x):
        b=self._basis(x)
        return np.column_stack([b[:,1:],(self.potential(np.atleast_2d(x))-b@self.projection)/self.sigma])

    def value(self,x):
        x=np.atleast_2d(x)
        return self.base.value(x)*np.exp(self._features(x)@self.theta-self.logZ)
