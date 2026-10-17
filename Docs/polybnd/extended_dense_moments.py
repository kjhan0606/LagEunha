"""Evaluate the same saved RK2 moment polynomial with extended arithmetic.

No new moments are inferred and no accepted state is changed. This reduces
within-step quantization of centered moments formed from global raw moments;
it cannot recover precision already absent from the stored source data.
"""
import numpy as np
from factor_moment_evolve import advance_factors_sylvester as advance


def dense_state(y,L,d0,dm,diag,mass,dt,s):
    extended=np.longdouble
    y=np.asarray(y,dtype=extended);d0=np.asarray(d0,dtype=extended);dm=np.asarray(dm,dtype=extended)
    phase=extended(s);path=y+extended(dt)*(phase*d0+phase*phase*(dm-d0))
    factors=advance(L,diag,dt*s) if s else [f.copy() for f in L]
    a=path[:-4].reshape(-1,20)
    for i,m in enumerate(mass):
        m=extended(m);q=a[i,3:5];Q=a[i,5:9].reshape(2,2);p=a[i,9:11];D=a[i,11:15].reshape(2,2)
        S=np.asarray(Q-np.outer(q,q)/m,dtype=float);C=np.asarray(D-np.outer(q,p)/m,dtype=float)
        if np.linalg.eigvalsh((S+S.T)/2).min()<=0:raise RuntimeError('nonpositive spatial covariance')
        B=np.linalg.solve(S,C).T
        Pi=np.outer(p,p)/m+np.asarray(B@S@B.T+factors[i]@factors[i].T,dtype=extended)
        a[i,15:19]=Pi.ravel()
    derivative=d0+2*phase*(dm-d0)
    return path,factors,derivative
