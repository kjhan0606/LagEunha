"""Experimental implicit midpoint source solve for the retained moment closure.

Positive-definite covariance is required at each iteration. Failure rejects
the source step; there is no clipping, energy repair, or tolerance fallback.
"""
import os
import numpy as np
import moment_ale_evolve as base
from gravity_source_2d import apply as predictor
from kinetic_factor_rhs import ALPHA
from quadratic_velocity_reconstruction import prepare,velocity_samples


def factors_from_raw(rows,mass,old_factors,frozen):
    result=[]
    for i,(a,m) in enumerate(zip(rows,mass)):
        if frozen[i]:result.append(np.array(old_factors[i],copy=True));continue
        q=a[3:5];Q=a[5:9].reshape(2,2);p=a[9:11];D=a[11:15].reshape(2,2);Pi=a[15:19].reshape(2,2)
        S=Q-np.outer(q,q)/m;C=D-np.outer(q,p)/m;B=np.linalg.solve(S,C).T
        R=Pi-np.outer(p,p)/m-B@S@B.T
        try:result.append(np.linalg.cholesky((R+R.T)/2))
        except np.linalg.LinAlgError as exc:raise RuntimeError('midpoint source covariance is not positive definite') from exc
    return result


def apply(t,y,mass,L,acceleration,dt,frozen_mask=None,density_provider=None,tolerance=2e-14,max_iterations=30):
    mask=np.zeros(len(mass),dtype=bool) if frozen_mask is None else np.asarray(frozen_mask,dtype=bool)
    if mask.shape!=(len(mass),) or not dt>0:raise ValueError('invalid midpoint source arguments')
    old=np.asarray(y);a=old[:-4].reshape(-1,20)
    cells,edges=base.geometry(a[:,:2],a[:,2],t)
    densities=[base.Density(p,m,row[3:5],row[5:9].reshape(2,2)) if density_provider is None else density_provider(i,p,np.r_[m,row[3:]]) for i,(p,m,row) in enumerate(zip(cells,mass,a))]
    forces=[];mixed=[];gravity=[]
    for d in densities:
        g=np.asarray(acceleration(d.points));gravity.append(g)
        forces.append(d.weighted_mass@g);mixed.append(d.points.T@(d.weighted_mass[:,None]*g))
    guess,_,_=predictor(t,old,mass,L,acceleration,dt,mask,density_provider)
    for iteration in range(max_iterations):
        mid=(old+guess)/2;b=mid[:-4].reshape(-1,20);midL=factors_from_raw(b,mass,L,mask)
        states=[base.decode(p,np.r_[m,row[3:]],density=d) for p,m,row,d in zip(cells,mass,b,densities)]
        if os.environ.get('LAGFORCE_QUADRATIC_VELOCITY')=='1':prepare(states,cells,edges,mass,midL)
        updated=old.copy();out=updated[:-4].reshape(-1,20)
        for i,s in enumerate(states):
            if mask[i]:continue
            d=densities[i];v,_=velocity_samples(s,midL[i],mass[i],d.points,ALPHA)
            cross=v.mean(1).T@(d.weighted_mass[:,None]*gravity[i]);delta=dt*(cross+cross.T)
            out[i,9:11]=a[i,9:11]+dt*forces[i];out[i,11:15]=a[i,11:15]+dt*mixed[i].ravel()
            out[i,15:19]=a[i,15:19]+delta.ravel();out[i,19]=a[i,19]+.5*np.trace(delta)
        error=float(np.max(abs(updated-guess)));guess=updated
        if error<tolerance:
            factors=factors_from_raw(out,mass,L,mask)
            thermal=out[:,19]-.5*(out[:,15]+out[:,18])
            if np.any(thermal<=0):raise RuntimeError('midpoint source lost positive thermal energy')
            return updated,factors,dict(iterations=iteration+1,residual=error,gravity_work=float((out[:,19]-a[:,19]).sum()))
    raise RuntimeError('midpoint source fixed point failed to converge')
