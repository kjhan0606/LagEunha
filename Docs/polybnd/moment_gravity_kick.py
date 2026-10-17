"""Positive quadrature gravity kick, a prerequisite for a coupled Kepler test.

This operates on supplied cell samples in two or three dimensions. It is NOT
wired into the Euler driver: drift/flux/potential-energy accounting must also
be implemented before claiming a fluid disk evolution.
"""
import numpy as np


def moments(points, mass_weights, velocity_samples, internal_energy):
    x=np.asarray(points,dtype=float);dm=np.asarray(mass_weights,dtype=float)
    v=np.asarray(velocity_samples,dtype=float)
    if x.ndim!=2 or v.ndim!=3 or v.shape[0]!=len(x) or v.shape[2]!=x.shape[1]:
        raise ValueError('incompatible position and velocity samples')
    if dm.shape!=(len(x),) or np.any(dm<0) or not dm.sum()>0 or not internal_energy>0:
        raise ValueError('positive mass and internal energy required')
    if not all(np.all(np.isfinite(a)) for a in (x,dm,v)) or not np.isfinite(internal_energy):
        raise ValueError('finite samples required')
    m=dm.sum();q=dm@x;Q=x.T@(dm[:,None]*x);u=v.mean(axis=1)
    p=dm@u;D=x.T@(dm[:,None]*u)
    Pi=np.einsum('n,nki,nkj->ij',dm,v,v)/v.shape[1]
    c=q/m;mean=p/m;dx=x-c
    S=dx.T@(dm[:,None]*dx);C=dx.T@(dm[:,None]*(u-mean))
    B=np.linalg.solve(S,C).T
    residual=v-mean-(dx@B.T)[:,None,:]
    weighted=(np.sqrt(dm/v.shape[1])[:,None,None]*residual).reshape(-1,x.shape[1])
    _,factor=np.linalg.qr(weighted,mode='reduced')
    return dict(m=m,q=q,Q=Q,p=p,D=D,Pi=Pi,E=internal_energy+.5*np.trace(Pi),
                internal_energy=internal_energy,residual_factor=factor.T)


def kick(points,mass_weights,velocity_samples,internal_energy,acceleration,dt):
    x=np.asarray(points,dtype=float);v=np.asarray(velocity_samples,dtype=float)
    g=np.asarray(acceleration(x),dtype=float)
    if g.shape!=x.shape or not np.all(np.isfinite(g)) or not np.isfinite(dt):
        raise ValueError('invalid gravity kick')
    new=v+dt*g[:,None,:]
    return moments(x,mass_weights,new,internal_energy),new
