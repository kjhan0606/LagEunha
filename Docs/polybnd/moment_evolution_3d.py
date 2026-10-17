"""Positive Sylvester midpoint stepping for the coupled 3D spatial operator.

Original masses are primary and fixed. Raw Pi is derived from the primary
covariance factor, as in the 2D integrator. No covariance eigenvalue clipping.
Dense moving-face path certification is not yet supplied by this driver.
"""
import numpy as np
from factor_moment_evolve import advance_factors_sylvester
from kinetic_moment_rhs_3d import rhs
from moment_state_3d import unpack


def rebuild(raw,factors):
    out=np.array(raw,copy=True)
    for i,(r,L) in enumerate(zip(out,factors)):
        m,q,Q,p,D,Pi,E=unpack(r);S=Q-np.outer(q,q)/m;C=D-np.outer(q,p)/m
        if np.linalg.eigvalsh((S+S.T)/2).min()<=0:raise RuntimeError('nonpositive spatial covariance')
        B=np.linalg.solve(S,C).T;Pi=np.outer(p,p)/m+B@S@B.T+L@L.T
        if E-.5*np.trace(Pi)<=0:raise RuntimeError('nonpositive thermal energy')
        out[i,25:34]=Pi.ravel()
    return out


def step(geometry,t,sites,weights,raw,factors,dt,bounds,boundary_velocity,face_order=3,pressure_reconstruction='constant',return_record=False):
    if not dt>0:raise ValueError('positive step required')
    dr,dx,dw,d0=rhs(geometry,sites,weights,raw,factors,bounds(t),boundary_velocity,face_order,pressure_reconstruction)
    halfL=advance_factors_sylvester(factors,d0,dt/2)
    halfraw=raw+dt/2*dr;halfraw[:,0]=raw[:,0];halfraw=rebuild(halfraw,halfL)
    dh,vh,zh,dm=rhs(geometry,sites+dt/2*dx,weights+dt/2*dw,halfraw,halfL,bounds(t+dt/2),boundary_velocity,face_order,pressure_reconstruction)
    nextL=advance_factors_sylvester(factors,dm,dt)
    nextraw=raw+dt*dh;nextraw[:,0]=raw[:,0];nextraw=rebuild(nextraw,nextL)
    diagnostics=dict(ledger_increment=dt*dm['ledger'],stage_mass_flux=max(d0['mass_flux'],dm['mass_flux']),generator_velocity_defect=0.,raw_Pi_completion=float(np.max(abs(nextraw[:,25:34]-(raw+dt*dh)[:,25:34]))))
    if return_record:
        diagnostics['dense_record']=dict(t=t,dt=dt,sites=sites.copy(),weights=weights.copy(),raw=raw.copy(),factors=np.asarray(factors).copy(),dx0=dx,dw0=dw,dr0=dr,dxmid=vh,dwmid=zh,drmid=dh,diag_mid=dm)
    return sites+dt*vh,weights+dt*zh,nextraw,np.array(nextL),diagnostics


def dense_state(record,s):
    if not 0<=s<=1:raise ValueError('dense path parameter outside step')
    dt=record['dt']
    def interpolate(key,first,middle):return record[key]+dt*(s*record[first]+s*s*(record[middle]-record[first]))
    sites=interpolate('sites','dx0','dxmid');weights=interpolate('weights','dw0','dwmid');raw=interpolate('raw','dr0','drmid')
    raw[:,0]=record['raw'][:,0]
    factors=advance_factors_sylvester(record['factors'],record['diag_mid'],s*dt) if s else record['factors'].copy()
    raw=rebuild(raw,factors)
    dx=record['dx0']+2*s*(record['dxmid']-record['dx0']);dw=record['dw0']+2*s*(record['dwmid']-record['dw0'])
    return sites,weights,raw,np.asarray(factors),dx,dw
