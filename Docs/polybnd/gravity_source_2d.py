"""Gravity source for the current 2D positive moment representation.

Supports its affine and covariance-budgeted quadratic velocity samples.
This is a frozen-position kick. A coupled drift and potential-energy flux
ledger are necessary for a full Kepler fluid integration.
"""
import os
import numpy as np
import moment_ale_evolve as base
from kinetic_factor_rhs import ALPHA
from conservative_gravity_kick import kick


def apply(t,flat,mass,factors,acceleration,dt,frozen_mask=None,density_provider=None):
    out=np.array(flat,copy=True);a=out[:-4].reshape(len(mass),20)
    frozen=np.zeros(len(mass),dtype=bool) if frozen_mask is None else np.asarray(frozen_mask,dtype=bool)
    if frozen.shape!=(len(mass),):raise ValueError('invalid frozen-cell mask')
    cells,edges=base.geometry(a[:,:2],a[:,2],t)
    raw=[np.r_[m,row[3:]] for m,row in zip(mass,a)]
    states=[base.decode(poly,r,density=None if density_provider is None else density_provider(i,poly,r)) for i,(poly,r) in enumerate(zip(cells,raw))]
    quadratic=os.environ.get('LAGFORCE_QUADRATIC_VELOCITY','0')=='1'
    from quadratic_velocity_reconstruction import prepare,velocity_samples
    if quadratic:prepare(states,cells,edges,mass,factors)
    next_factors=[];diagnostics=[]
    for i,(s,m,L) in enumerate(zip(states,mass,factors)):
        if frozen[i]:next_factors.append(np.array(L,copy=True));continue
        r=a[i];raw=dict(m=m,q=r[3:5].copy(),Q=r[5:9].reshape(2,2).copy(),p=r[9:11].copy(),D=r[11:15].reshape(2,2).copy(),Pi=r[15:19].reshape(2,2).copy(),E=r[19])
        x=s['density'].points;dm=s['density'].weighted_mass
        v,_=velocity_samples(s,L,m,x,ALPHA)
        new,ln,d=kick(raw,x,dm,v,acceleration,dt)
        r[9:11]=new['p'];r[11:15]=new['D'].ravel();r[15:19]=new['Pi'].ravel();r[19]=new['E']
        next_factors.append(ln);diagnostics.append(d)
    return out,next_factors,dict(gravity_work=sum(d['gravity_work'] for d in diagnostics),max_quadrature_moment_mismatch=max((d['quadrature_moment_mismatch'] for d in diagnostics),default=0.),max_Pi_roundoff_completion=max((d['Pi_roundoff_completion'] for d in diagnostics),default=0.),max_thermal_change=max((abs(d['thermal_change']) for d in diagnostics),default=0.),angular_source=sum(d['angular_source'][0,1] for d in diagnostics))
