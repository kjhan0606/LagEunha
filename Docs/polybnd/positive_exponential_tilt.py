"""Positive finite-quadrature moment fit; no clipping of target moments."""
import numpy as np


def fit(weights,features,target,tolerance=2e-12,max_iterations=80):
    weights=np.asarray(weights,dtype=float);F=np.asarray(features,dtype=float);target=np.asarray(target,dtype=float)
    if weights.ndim!=1 or F.shape!=(len(weights),len(target)) or np.any(weights<=0) or not np.isfinite(weights).all() or not np.isfinite(F).all() or not np.isfinite(target).all():
        raise ValueError('invalid positive moment fitting data')
    logw=np.log(weights);theta=np.zeros(len(target))
    def evaluate(th):
        z=logw+F@th;shift=z.max();p=np.exp(z-shift);Z=p.sum();p/=Z
        mean=p@F;centered=F-mean;H=centered.T@(p[:,None]*centered)
        return mean-target,H,p,shift+np.log(Z)
    for iteration in range(max_iterations):
        residual,H,p,logZ=evaluate(theta);norm=float(np.max(abs(residual)))
        if norm<tolerance:return theta,logZ,p,dict(iterations=iteration,residual=norm)
        step=np.linalg.solve(H,residual)
        for bt in range(32):
            candidate=theta-step*2.**(-bt)
            if np.max(abs(evaluate(candidate)[0]))<norm:theta=candidate;break
        else:raise RuntimeError('positive moment tilt stalled; target not modified')
    raise RuntimeError('positive moment tilt did not converge; target not modified')
