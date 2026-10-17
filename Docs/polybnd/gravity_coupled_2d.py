"""Strang gravity kicks around the current positive midpoint fluid step.

Caller supplies the boundary law of each split subflow. Current fixed-mass
velocity/pressure reconstruction is retained through environment settings.
The global potential-energy diagnostic is explicit; no energy correction is
applied. Dense path and frozen Kepler boundary gates remain separate work.
"""
import numpy as np
import moment_ale_evolve as base
from gravity_source_2d import apply
from kinetic_factor_rhs import rhs
from factor_moment_evolve import rebuild,advance_factors_sylvester


def mask_at(frozen,t,y,mass):
    if frozen is None:return np.zeros(len(mass),dtype=bool)
    value=frozen(t,y[:-4].reshape(len(mass),20)[:,:2]) if callable(frozen) else frozen
    result=np.asarray(value,dtype=bool)
    if result.shape!=(len(mass),):raise ValueError('invalid frozen-cell mask')
    return result


def restore_frozen(new,newL,old,oldL,mask,mass):
    if np.any(mask):
        new[:-4].reshape(len(mass),20)[mask]=old[:-4].reshape(len(mass),20)[mask]
        for i in np.flatnonzero(mask):newL[i]=np.array(oldL[i],copy=True)


def hydro_step(t,y,mass,L,dt,frozen=None):
    first=mask_at(frozen,t,y,mass)
    dy,d=rhs(t,y,mass,covariance_factors=L,frozen_mask=first)
    halfL=advance_factors_sylvester(L,d,dt/2);half=rebuild(y+dt*dy/2,mass,halfL)
    restore_frozen(half,halfL,y,L,first,mass)
    middle=mask_at(frozen,t+dt/2,half,mass)
    dh,dm=rhs(t+dt/2,half,mass,covariance_factors=halfL,frozen_mask=middle)
    nextL=advance_factors_sylvester(L,dm,dt);new=rebuild(y+dt*dh,mass,nextL)
    restore_frozen(new,nextL,y,L,middle,mass)
    return new,nextL,dict(stage_mass_flux=max(d['mass_flux'],dm['mass_flux']),generator_velocity_defect=max(d['generator_velocity_defect'],dm['generator_velocity_defect']))


def step(t,y,mass,L,dt,acceleration,frozen=None):
    if not dt>0:raise ValueError('positive split step required')
    kicked,halfL,g0=apply(t,y,mass,L,acceleration,dt/2,mask_at(frozen,t,y,mass))
    drift,driftL,d=hydro_step(t,kicked,mass,halfL,dt,frozen)
    result,nextL,g1=apply(t+dt,drift,mass,driftL,acceleration,dt/2,mask_at(frozen,t+dt,drift,mass))
    return result,nextL,dict(**d,gravity_work=g0['gravity_work']+g1['gravity_work'],angular_source=g0['angular_source']+g1['angular_source'],max_thermal_kick_change=max(g0['max_thermal_change'],g1['max_thermal_change']),max_kick_completion=max(g0['max_Pi_roundoff_completion'],g1['max_Pi_roundoff_completion']))


def potential_energy(t,y,mass,potential):
    a=np.asarray(y)[:-4].reshape(len(mass),20);cells,_=base.geometry(a[:,:2],a[:,2],t);total=0.
    for p,m,row in zip(cells,mass,a):
        density=base.decode(p,np.r_[m,row[3:]])['density'];total+=density.weighted_mass@potential(density.points)
    return float(total)
