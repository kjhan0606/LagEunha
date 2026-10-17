"""Minimize normalized moment mismatch with a small, always feasible LP."""
import numpy as np
from scipy import sparse
from scipy.optimize import linprog


def box_dual_lower_bound(matrix,rhs,cost,limits,multipliers):
    """Repair inequality duals; eliminate finite box variables analytically."""
    lam=np.minimum(0.,np.asarray(multipliers))
    lam*= (1-1e-12)/max(1.,float(-lam.sum()))
    reduced=np.asarray(cost-matrix.T@lam)
    if reduced[-1]<0:raise RuntimeError('unbounded slack in repaired dual')
    lower=np.array([v[0] for v in limits[:-1]]);upper=np.array([v[1] for v in limits[:-1]])
    terms=np.minimum(reduced[:-1]*lower,reduced[:-1]*upper)
    raw=float(rhs@lam+terms.sum())
    size=len(cost)+len(rhs)+64
    gamma=size*np.finfo(float).eps/(1-size*np.finfo(float).eps)
    coefficient_bound=abs(np.asarray(cost[:-1]))+np.asarray(abs(matrix[:,:-1]).T@abs(lam))
    margin=float(gamma*(abs(rhs)@abs(lam)+(coefficient_bound*np.maximum(abs(lower),abs(upper))).sum()+1))
    return dict(lower_bound=raw-margin,roundoff_margin=margin,slack_reduced_cost=float(reduced[-1]),multipliers=lam)


def minimum_moment_residual(A,b,bounds):
    A=np.asarray(A);b=np.asarray(b);k=A.shape[1];rows=len(b)
    matrix=sparse.vstack([sparse.hstack([sparse.csr_matrix(A),-np.ones((rows,1))]),
        sparse.hstack([-sparse.csr_matrix(A),-np.ones((rows,1))])],format='csr')
    rhs=np.r_[b,-b];cost=np.r_[np.zeros(k),1.]
    limits=list(bounds)+[(0.,None)]
    sol=linprog(cost,A_ub=matrix,b_ub=rhs,bounds=limits,method='highs',
        options={'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10})
    if not sol.success:return dict(classification='numerical_failure',solver_status=int(sol.status),reason=sol.message)
    u=sol.x[:k];actual=float(np.max(abs(A@u-b)))
    lower=np.array([lo for lo,hi in limits]);upper=np.array([np.inf if hi is None else hi for lo,hi in limits])
    dual=float(rhs@sol.ineqlin.marginals+lower@sol.lower.marginals+
        upper[np.isfinite(upper)]@sol.upper.marginals[np.isfinite(upper)])
    stationary=float(np.max(abs(cost-matrix.T@sol.ineqlin.marginals-sol.lower.marginals-sol.upper.marginals)))
    gap=abs(float(sol.fun)-dual)
    sign=float(max(0.,np.max(sol.ineqlin.marginals),-np.min(sol.lower.marginals),np.max(sol.upper.marginals)))
    bound=float(max(0.,np.max(lower-sol.x),np.max(sol.x-upper)))
    repaired=box_dual_lower_bound(matrix,rhs,cost,limits,sol.ineqlin.marginals)
    validated=stationary<1e-8 and gap<1e-8 and sign<1e-8 and bound<1e-8 and actual<=sol.x[-1]+1e-8
    classification='numerical_failure'
    if validated and actual<1e-8:classification='moment_witness'
    elif repaired['lower_bound']>1e-8:classification='finite_basis_infeasible_certificate'
    return dict(classification=classification,minimum_residual=float(sol.fun),actual_moment_residual=actual,
        mass_moment_residual=float(abs(A[0]@u-b[0])),
        dual_lower_bound=repaired['lower_bound'],raw_solver_dual=dual,dual_roundoff_margin=repaired['roundoff_margin'],
        duality_gap=gap,stationarity_error=stationary,dual_sign_error=sign,bound_error=bound,
        density_units=u,inequality_dual=sol.ineqlin.marginals,lower_dual=sol.lower.marginals,upper_dual=sol.upper.marginals)
