"""Positive m/q/Q/Psi relative entropy reconstruction candidate.

Psi is an explicit independent sidecar. Does not modify raw18 or dynamics.
Newton failure is unresolved, never an assertion of physical infeasibility.
"""
import numpy as np
from contact_prior_density import clip_halfplane,polygon_area
from positive_radial_quadrature import rule

class ReconstructionFailure(RuntimeError):
    def __init__(self,status,**details):
        self.status=status;self.details=details;super().__init__(status)

def features(points,origin,scale):
    z=(np.asarray(points,dtype=np.longdouble)-np.asarray(origin,dtype=np.longdouble))/np.longdouble(scale)
    return np.column_stack([np.ones(len(z)),z,z[:,0]**2,z[:,0]*z[:,1],z[:,1]**2])

class ContactPotentialDensity:
    def __init__(self,poly,m,q,Q,prior_metadata,baseline,potential,Psi_target,
                 reference_coefficients,order=24,radii=(),polynomial=False):
        self.mass=float(m);self.origin=baseline.origin.copy();self.scale=float(baseline.scale)
        self.area=float(baseline.area)
        self.cut_metadata=dict(prior_metadata);self.potential=potential;self.radii=tuple(radii)
        self.reference=np.asarray(reference_coefficients,dtype=np.longdouble)
        if self.reference.shape!=(6,):raise ValueError('six coefficients in [1,zx,zy,zx²,zxzy,zy²] required')
        self.Psi_target=np.longdouble(Psi_target);self.parts=[];self.part_priors=[]
        if not np.isfinite(self.mass) or self.mass<=0 or not np.isfinite(self.Psi_target) or not np.isfinite(q).all() or not np.isfinite(Q).all():
            raise ValueError('finite positive mass, moments and potential target required')
        curved=hasattr(poly,'positive_quadrature')
        if prior_metadata['flat']:pieces=[(poly if curved else np.asarray(poly),prior_metadata['value'])]
        else:
            n=np.asarray(prior_metadata['normal']);o=prior_metadata['offset']
            if curved:
                # Domain API clips n.x<=o; polygon helper clips n.x>=o.
                pieces=[(poly.clip_halfplane(n,o),prior_metadata['low']),(poly.clip_halfplane(-n,-o),prior_metadata['high'])]
            else:pieces=[(clip_halfplane(poly,-n,-o),prior_metadata['low']),(clip_halfplane(poly,n,o),prior_metadata['high'])]
        for p,rho in pieces:
            if not np.isfinite(rho) or rho<=0:raise ValueError('positive finite prior required')
            if p is not None and (p.area if curved else polygon_area(p))>0:self.parts.append(p);self.part_priors.append(float(rho))
        self.points,self.weights,self.prior=self.quadrature(order)
        A=features(self.points,self.origin,self.scale);self.A=A
        ld=np.longdouble;center=np.asarray(q,dtype=ld)/ld(m);o=self.origin.astype(ld)
        second=np.asarray(Q,dtype=ld)/ld(m)-np.outer(o,center)-np.outer(center,o)+np.outer(o,o)
        self.target_A=np.r_[ld(1),(center-o)/ld(self.scale),second[0,0]/ld(self.scale)**2,second[0,1]/ld(self.scale)**2,second[1,1]/ld(self.scale)**2]
        oldrho=baseline._value_with_prior(self.points,self.prior)
        oldweights=self.weights*oldrho;prob=oldweights/oldweights.sum();self.baseline_prob=prob
        sqrtp=np.sqrt(prob);design=sqrtp[:,None]*np.asarray(A,dtype=float)
        U,R=np.linalg.qr(design,mode='reduced');singular=np.linalg.svd(R,compute_uv=False)
        qr_error=(len(design)+6)*np.finfo(float).eps*singular[0]
        if singular[-1]<=qr_error:raise ReconstructionFailure('unresolved_polynomial_rank',singular=singular.tolist(),error=float(qr_error))
        Phi=np.asarray(potential(self.points.astype(ld)),dtype=ld)
        remainder=Phi-A@self.reference
        beta=np.linalg.solve(R,U.T@(sqrtp*np.asarray(remainder,dtype=float)))
        self.remainder_projection=beta.astype(ld)
        self.projection=self.reference+beta.astype(ld)
        chi=remainder-A@beta.astype(ld)
        self.sigma=ld(np.sqrt(np.sum(prob.astype(ld)*chi*chi)))
        phi_eval_error=64*np.finfo(ld).eps*np.max(abs(Phi)+abs(A)@abs(self.reference))
        projection_error=(len(design)+12)*np.finfo(float).eps*np.linalg.cond(R)*(np.max(abs(remainder))+np.max(abs(A)@abs(beta)))
        self.rank_error=ld(phi_eval_error+projection_error)
        self.polynomial_dependent=bool(polynomial)
        projected_target=self.target_A@self.projection
        self.transformed_target=(self.Psi_target/ld(m)-projected_target)
        if polynomial:
            if np.max(abs(remainder))>phi_eval_error:raise ReconstructionFailure('invalid_declared_polynomial')
            target_bound=64*np.finfo(ld).eps*(abs(self.Psi_target/ld(m))+abs(self.target_A)@abs(self.reference))
            if abs(self.Psi_target/ld(m)-self.target_A@self.reference)>target_bound:
                raise ReconstructionFailure('inconsistent_dependent_potential_target',mismatch=float(self.Psi_target/ld(m)-self.target_A@self.reference),bound=float(target_bound))
            G=np.asarray(A[:,1:],dtype=float);target=np.asarray(self.target_A[1:],dtype=float)
        else:
            if not self.sigma>self.rank_error:raise ReconstructionFailure('unresolved_potential_rank',sigma=float(self.sigma),error=float(self.rank_error))
            G=np.column_stack([A[:,1:],chi/self.sigma]).astype(float)
            target=np.r_[self.target_A[1:],self.transformed_target/self.sigma].astype(float)
        self.target=target;self.design=G
        for k,t in enumerate(target):
            lo=float(G[:,k].min());hi=float(G[:,k].max())
            if not lo<=t<=hi:
                raise ReconstructionFailure('outside_quadrature_moment_hull',feature=k,target=float(t),minimum=lo,maximum=hi,witness='coordinate separating hyperplane')
        theta=np.r_[baseline.theta,0.] if not polynomial else baseline.theta.copy()
        def evaluate(v):
            exponent=G@v+np.log(self.prior);shift=float(exponent.max())
            wm=self.weights*np.exp(exponent-shift);Z=float(wm.sum());p=wm/Z
            residual=p@G-target;delta=G-p@G;H=delta.T@(p[:,None]*delta)
            return residual,H,shift+np.log(Z),p
        self.initial_dual_residual=float(abs(evaluate(theta)[0]).max())
        for iteration in range(60):
            residual,H,logZ,prob=evaluate(theta);norm=float(abs(residual).max())
            if norm<3e-13:break
            try:step=np.linalg.solve(H,residual)
            except np.linalg.LinAlgError:raise ReconstructionFailure('unresolved_singular_dual_hessian',residual=norm,iteration=iteration)
            for damping in range(30):
                trial=theta-step*2.**(-damping)
                if abs(evaluate(trial)[0]).max()<norm:theta=trial;break
            else:raise ReconstructionFailure('unresolved_newton_stall',residual=norm,iteration=iteration)
        else:raise ReconstructionFailure('unresolved_newton_limit',residual=norm)
        self.theta=theta;self.logZ=logZ;self.weighted_mass=self.mass*prob
        if not np.isfinite(prob).all() or np.any(prob<=0):raise ReconstructionFailure('unresolved_density_underflow')
        self.residual=norm;self.iterations=iteration+1;self.condition=float(np.linalg.cond(H))
        self.entropy=self.mass*(np.log(self.mass)+float(theta@target)-logZ)

    def quadrature(self,order):
        values=[p.positive_quadrature(order,self.radii) if hasattr(p,'positive_quadrature') else rule(p,order,self.radii) for p in self.parts]
        return np.concatenate([v[0] for v in values]),np.concatenate([v[1] for v in values]),np.concatenate([np.full(len(v[1]),rho) for v,rho in zip(values,self.part_priors)])

    def _value_with_prior(self,points,prior):
        G=self._evaluation_features(points)
        value=self.mass*np.exp(G@self.theta+np.log(prior)-self.logZ)
        if not np.isfinite(value).all() or np.any(value<=0):raise ReconstructionFailure('unresolved_nonpositive_density')
        return value

    def _evaluation_features(self,points):
        A=features(points,self.origin,self.scale)
        if self.polynomial_dependent:G=np.asarray(A[:,1:],dtype=float)
        else:
            chi=np.asarray(self.potential(np.asarray(points,dtype=np.longdouble)),dtype=np.longdouble)-A@self.reference-A@self.remainder_projection
            G=np.column_stack([A[:,1:],chi/self.sigma]).astype(float)
        return G

    def value(self,points):
        points=np.atleast_2d(points);meta=self.cut_metadata
        prior=np.full(len(points),meta['value']) if meta['flat'] else np.where(points@np.asarray(meta['normal'])>=meta['offset'],meta['high'],meta['low'])
        return self._value_with_prior(points,prior)

    def value_trace(self,points,inward_direction):
        points=np.atleast_2d(points);meta=self.cut_metadata
        if meta['flat']:prior=np.full(len(points),meta['value'])
        else:
            n=np.asarray(meta['normal']);s=points@n-meta['offset'];direction=np.broadcast_to(inward_direction,points.shape)
            if not np.isfinite(direction).all():raise ValueError('finite inward direction required')
            derivative=direction@n
            if np.any((s==0)&(derivative==0)):raise ValueError('tangent trace is undefined')
            prior=np.where((s>0)|((s==0)&(derivative>0)),meta['high'],meta['low'])
        return self._value_with_prior(points,prior)

    def value_trace_batch(self,points,inward_direction,return_branch=False):
        """One-point dot kernels and exact cut predicates; all six features.

        No call to the five-feature parent implementation is permitted here.
        Branch codes: flat0, high1, low-1. An on-cut tangent is undefined.
        """
        points=np.asarray(points,dtype=float)
        if points.ndim!=2 or points.shape[1]!=2:raise ValueError('trace batch requires (nodes,2) points')
        meta=self.cut_metadata
        if meta['flat']:
            prior=np.full(len(points),meta['value']);branch=np.zeros(len(points),dtype=np.int8)
        else:
            normal=np.asarray(meta['normal'])
            signed=(points[:,None,:]@normal)[:,0]-meta['offset']
            direction=np.broadcast_to(np.asarray(inward_direction,dtype=float),points.shape)
            if not np.isfinite(direction).all():raise ValueError('finite inward direction required')
            derivative=(direction[:,None,:]@normal)[:,0]
            oncut=signed==0
            if np.any(oncut&(derivative==0)):raise ValueError('tangent trace is undefined')
            high=(signed>0)|(oncut&(derivative>0))
            prior=np.where(high,meta['high'],meta['low']);branch=np.where(high,1,-1).astype(np.int8)
        G=self._evaluation_features(points)
        value=self.mass*np.exp((G[:,None,:]@self.theta)[:,0]+np.log(prior)-self.logZ)
        if not np.isfinite(value).all() or np.any(value<=0):raise ReconstructionFailure('unresolved_nonpositive_density')
        return (value,branch) if return_branch else value
