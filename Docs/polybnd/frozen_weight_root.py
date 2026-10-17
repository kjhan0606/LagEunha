"""Active-cell mass constraints with frozen weight rates as Dirichlet data."""
import warnings
import numpy as np
from scipy.sparse.linalg import spsolve,lsmr,MatrixRankWarning
from vector_kinetic_root import record_arrays,evaluate


def solve(records,count,frozen,tolerance=2e-13):
    frozen=np.asarray(frozen,dtype=bool)
    if frozen.shape!=(count,):raise ValueError('invalid frozen-cell mask')
    active=np.flatnonzero(~frozen);data=record_arrays(records);z=np.zeros(count)
    for iteration in range(50):
        residual,J,values=evaluate(data,count,z);norm=float(np.max(abs(residual[active]),initial=0.))
        if norm<tolerance:
            a,b,wn=values;return z,[(u,v,float(s)) for u,v,s in zip(a,b,wn)],norm
        matrix=J[active,:][:,active];right=residual[active]
        with warnings.catch_warnings():
            warnings.simplefilter('error',MatrixRankWarning)
            try:step=spsolve(matrix,right)
            except MatrixRankWarning:step=lsmr(matrix,right,atol=1e-14,btol=1e-14,maxiter=max(200,4*len(active)))[0]
        if not np.all(np.isfinite(step)):raise RuntimeError('nonfinite frozen-boundary weight step')
        for k in range(30):
            trial=z.copy();trial[active]+=step*2.**(-k);r,_,_=evaluate(data,count,trial,False)
            if np.max(abs(r[active]))<norm:z=trial;break
        else:raise RuntimeError('frozen-boundary mass constraint line search failed')
    raise RuntimeError('frozen-boundary mass constraint did not converge')
