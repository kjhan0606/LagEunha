"""Vector/sparse solution of the existing discrete-velocity mass equations.

Same quadrature records and monotone residual, with no density/model fallback.
Only linear algebra and summation order differ from the legacy dense solve.
"""
import warnings
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve,lsmr,MatrixRankWarning
from scipy.sparse.csgraph import connected_components


def record_arrays(records):
    return dict(i=np.asarray([r[0] for r in records],dtype=int),j=np.asarray([r[1] for r in records],dtype=int),
        scale=np.asarray([r[2] for r in records],dtype=float),base=np.asarray([r[3] for r in records],dtype=float),
        distance=np.asarray([r[4] for r in records],dtype=float),
        vl=np.asarray([r[5]@r[9] for r in records]),vr=np.asarray([r[6]@r[9] for r in records]),
        rl=np.asarray([r[7] for r in records]),rr=np.asarray([r[8] for r in records]),
        epsilon=np.asarray([r[15] if len(r)>15 and np.ndim(r[15])==0 else 0. for r in records]))


def evaluate(data,count,z,with_jacobian=True):
    i=data['i'];j=data['j'];internal=j>=0
    inv=np.divide(1.,2*data['distance'],out=np.zeros_like(data['distance']),where=internal)
    wn=data['base']+np.where(internal,(z[i]-z[j])*inv,0.)
    sl=data['vl']-wn[:,None];sr=data['vr']-wn[:,None]
    epsilon=np.asarray(data.get('epsilon',0.))
    smooth=bool(np.any(epsilon))
    if smooth:
        from kinetic_speed_smoothing import split
        epsilon=np.broadcast_to(epsilon,sl.shape[:-1])[:,None]
        lp,_,ld,_=split(sl,epsilon);_,rm,_,rd=split(sr,epsilon)
        wl=.25*data['rl'][:,None]*lp;wr=.25*data['rr'][:,None]*rm
    else:
        wl=.25*data['rl'][:,None]*np.maximum(sl,0.);wr=.25*data['rr'][:,None]*np.minimum(sr,0.)
    flux=data['scale']*(wl.sum(1)+wr.sum(1))
    residual=np.bincount(i,weights=flux,minlength=count)-np.bincount(j[internal],weights=flux[internal],minlength=count)
    jacobian=None
    if with_jacobian:
        rho=.25*(data['rl']*np.where(sl>0,1.,np.where(sl<0,0.,.5)).sum(1)+
            data['rr']*np.where(sr<0,1.,np.where(sr>0,0.,.5)).sum(1))
        if smooth:rho=.25*(data['rl']*ld.sum(1)+data['rr']*rd.sum(1))
        c=(data['scale']*rho*inv)[internal];a=i[internal];b=j[internal]
        rows=np.r_[a,b,a,b];columns=np.r_[a,b,b,a];values=np.r_[c,c,-c,-c]
        jacobian=coo_matrix((values,(rows,columns)),shape=(count,count)).tocsc()
    return residual,jacobian,(wl,wr,wn)


def solve_arrays(data,count,initial=None,tolerance=2e-13):
    z=np.zeros(count) if initial is None else np.asarray(initial,dtype=float).copy();z-=z[-1]
    for iteration in range(50):
        residual,L,values=evaluate(data,count,z);norm=float(np.max(abs(residual)))
        if norm<tolerance:return z,values,norm,iteration
        L.eliminate_zeros();matrix=L[:-1,:-1];rhs=residual[:-1]
        components,_=connected_components(L,directed=False)
        with warnings.catch_warnings():
            warnings.simplefilter('error',MatrixRankWarning)
            try:
                if components>1:raise MatrixRankWarning('separate kinetic weight gauges')
                step=spsolve(matrix,rhs)
            except MatrixRankWarning:step=lsmr(matrix,rhs,atol=1e-14,btol=1e-14,maxiter=max(200,4*count))[0]
        if not np.all(np.isfinite(step)):raise RuntimeError('nonfinite sparse kinetic Newton step')
        dz=np.r_[step,0.]
        for k in range(30):
            trial=z+dz*2.**(-k)
            value,_,_=evaluate(data,count,trial,False)
            if np.max(abs(value))<norm:z=trial;break
        else:raise RuntimeError('vector kinetic mass constraint line search failed')
    raise RuntimeError('vector kinetic mass constraint did not converge')


def solve_records(records,count):
    z,(wl,wr,wn),norm,iterations=solve_arrays(record_arrays(records),count)
    return z,[(a,b,float(v)) for a,b,v in zip(wl,wr,wn)],norm
