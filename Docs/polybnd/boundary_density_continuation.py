"""Certified reuse of nearby density witnesses; full LP remains the fallback.

The old dual multipliers yield a lower bound for the new bulk minimum. A
moment-corrected positive witness is reused only below that lower bound plus
the original budget. Every face density is still imposed exactly.
"""
import numpy as np


def prior_witness(p,budget=1e-5):
    """A direct feasible witness certified by the universal L1 lower bound0."""
    k=p['k'];u=p['prior_density_units'].copy();bounds=np.asarray(p['bounds'][:k]);lo,hi=bounds.T;fixed=lo==hi
    u[fixed]=lo[fixed]
    free=np.flatnonzero(~fixed & (u>lo+1e-12) & (u<hi-1e-12))
    if len(free)<6:return None
    u[free]+=np.linalg.lstsq(p['A'][:,free],p['b']-p['A']@u,rcond=1e-12)[0]
    moment=float(np.max(abs(p['A']@u-p['b'])));bound=max(0.,float(np.max(lo-u)),float(np.max(u-hi)))
    objective=p['spatial_objective'](u)
    if moment>1e-10 or bound>1e-12 or not np.isfinite(objective) or objective>budget:return None
    rho=u*p['high']
    fit=dict(feasible=True,objective=objective,moment_error=moment,min_density=float(rho.min()),max_density=float(rho.max()),mass_fractions=rho*p['raw'][:,0]/p['mass'])
    baseline=dict(objective=0.,certificate=dict(equality=np.zeros(6),inequality=np.zeros(3*k)))
    return fit,baseline,u


def lower_bound(p,certificate):
    k=p['k'];y=np.asarray(certificate['equality']);lam=np.minimum(0.,np.asarray(certificate['inequality']).reshape(3,k))
    if y.shape!=(6,) or lam.shape!=(3,k):raise ValueError('certificate dimension mismatch')
    # Each epigraph variable has coefficient -1 in its three inequalities.
    # Scale its nonpositive multipliers to leave a nonnegative reduced cost.
    total=-lam.sum(0);scale=np.divide(p['c'][k:]*(1-1e-12),total,out=np.zeros(k),where=total>0)
    lam=lam*scale[None,:];multipliers=lam.ravel()
    reduced=p['c']-p['Aeq'].T@y-p['Aub'].T@multipliers
    if np.min(reduced[k:]) < -1e-18:raise RuntimeError('invalid epigraph dual bound')
    lower=np.full(k,p['low']*.99/p['high']);upper=np.full(k,1.01)
    bound=float(p['b']@y+p['bub']@multipliers+np.minimum(reduced[:k]*lower,reduced[:k]*upper).sum())
    return bound


def witness(p,previous,vertices,budget=1e-5):
    k=p['k']
    if previous is None or previous['u'].shape!=(k,) or previous['vertices'].shape!=vertices.shape:return None
    if np.max(abs(vertices-previous['vertices'])) > .05*np.sqrt(p['raw'][:,0].sum()):return None
    u=previous['u']*previous['high']/p['high'];u=u.copy()
    bounds=np.asarray(p['bounds'][:k]);lo,hi=bounds.T;fixed=lo==hi
    u[fixed]=lo[fixed]
    free=np.flatnonzero(~fixed & (u>lo+1e-12) & (u<hi-1e-12))
    if len(free)<6:return None
    delta=np.linalg.lstsq(p['A'][:,free],p['b']-p['A']@u,rcond=1e-12)[0];u[free]+=delta
    moment=float(np.max(abs(p['A']@u-p['b'])))
    rho=u*p['high'];bound=max(0.,float(np.max(lo-u)),float(np.max(u-hi)))
    if moment>1e-10 or bound>1e-12:return None
    objective=p['spatial_objective'](u);base_lower=lower_bound(p,previous['certificate'])
    if not np.isfinite(base_lower) or objective>base_lower+budget+1e-10:return None
    fit=dict(feasible=True,objective=objective,moment_error=moment,min_density=float(rho.min()),max_density=float(rho.max()),
        mass_fractions=rho*p['raw'][:,0]/p['mass'])
    return fit,dict(objective=base_lower,certificate=previous['certificate']),u
