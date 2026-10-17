"""Tiny closed-form audit, not a fluid evolution: Rdot=-k R+q I."""
import json
import numpy as np
from factor_moment_evolve import advance_factors

def audit():
    k=2.;q=.3;r0=.2;errors=[]
    d={'cov_A':[np.zeros((2,2))],'cov_K':[k*np.eye(2)],'cov_Qfactor':[np.sqrt(q)*np.eye(2)]}
    for dt in (.02,.01,.005):
        L=advance_factors([np.sqrt(r0)*np.eye(2)],d,dt)[0]
        exact=(np.exp(-k*dt)*r0+q*(-np.expm1(-k*dt))/k)*np.eye(2)
        errors.append(float(np.linalg.norm(L@L.T-exact)))
    orders=[float(np.log2(a/b)) for a,b in zip(errors[:-1],errors[1:])]
    assert all(1.95<x<2.05 for x in orders)
    return dict(local_errors=errors,local_orders=orders,conclusion='Local error O(dt^2): kinetic right-factor loss extension is generally only first order globally. Full coupled problem may mask this defect.')
if __name__=='__main__':print(json.dumps(audit(),indent=2))
