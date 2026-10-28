"""Five-feature curved-domain baseline; used only to precondition Psi fit.

Shared exponential/trace API is inherited from the tested polygon baseline.
The actual physical density is the six-feature ContactPotentialDensity.
"""
import numpy as np
from contact_prior_density import ContactPriorDensity,features

class CurvedContactPriorDensity(ContactPriorDensity):
    def __init__(self,domain,m,q,Q,prior_metadata,order=12,radii=(.25,.5)):
        self.mass=float(m);self.area=float(domain.area);self.origin=domain.origin.copy();self.scale=np.sqrt(self.area)
        self.cut_metadata=dict(prior_metadata);self.quadrature_order=order;self.radii=tuple(radii)
        if not np.isfinite(self.mass) or self.mass<=0 or self.area<=0:raise ValueError('positive finite cell mass and area required')
        if prior_metadata['flat']:pieces=[(domain,float(prior_metadata['value']))]
        else:
            normal=np.array(prior_metadata['normal']);offset=prior_metadata['offset']
            pieces=[(domain.clip_halfplane(normal,offset),prior_metadata['low']),
                    (domain.clip_halfplane(-normal,-offset),prior_metadata['high'])]
        self.parts=[];self.part_priors=[];all_points=[];all_weights=[];priors=[]
        for part,rho in pieces:
            if not np.isfinite(rho) or rho<=0:raise ValueError('positive finite prior required')
            if part is None:continue
            points,weights=part.positive_quadrature(order,self.radii)
            self.parts.append(part);self.part_priors.append(float(rho));all_points.append(points);all_weights.append(weights);priors.append(np.full(len(weights),rho))
        self.points=np.concatenate(all_points);self.weights=np.concatenate(all_weights);prior=np.concatenate(priors)
        F=features(self.points,self.origin,self.scale);ld=np.longdouble;o=self.origin.astype(ld);center=np.asarray(q,dtype=ld)/ld(m)
        second=np.asarray(Q,dtype=ld)/ld(m)-np.outer(o,center)-np.outer(center,o)+np.outer(o,o)
        target=np.asarray(np.r_[(center-o)/ld(self.scale),second[0,0]/ld(self.area),second[0,1]/ld(self.area),second[1,1]/ld(self.area)],dtype=float)
        theta=np.zeros(5)
        def evaluate(t):
            exponent=F@t+np.log(prior);shift=float(exponent.max());wm=self.weights*np.exp(exponent-shift);Z=float(wm.sum());p=wm/Z
            residual=p@F-target;delta=F-p@F;H=delta.T@(p[:,None]*delta)
            return residual,H,shift+np.log(Z),p
        for iteration in range(60):
            residual,H,logZ,prob=evaluate(theta);norm=float(abs(residual).max())
            if norm<3e-13:break
            step=np.linalg.solve(H,residual)
            for damping in range(30):
                trial=theta-step*2.**(-damping)
                if abs(evaluate(trial)[0]).max()<norm:theta=trial;break
            else:raise RuntimeError('curved baseline moment fit stalled')
        else:raise RuntimeError('curved baseline moment fit did not converge')
        self.theta=theta;self.logZ=logZ;self.residual=norm;self.iterations=iteration+1;self.condition=float(np.linalg.cond(H));self.weighted_mass=self.mass*prob
        if not np.isfinite(prob).all() or np.any(prob<=0):raise RuntimeError('unresolved curved baseline density')
