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


def hydro_step(t,y,mass,L,dt,frozen=None,return_record=False):
    first=mask_at(frozen,t,y,mass)
    dy,d=rhs(t,y,mass,covariance_factors=L,frozen_mask=first)
    halfL=advance_factors_sylvester(L,d,dt/2);half=rebuild(y+dt*dy/2,mass,halfL)
    restore_frozen(half,halfL,y,L,first,mass)
    middle=mask_at(frozen,t+dt/2,half,mass)
    dh,dm=rhs(t+dt/2,half,mass,covariance_factors=halfL,frozen_mask=middle)
    nextL=advance_factors_sylvester(L,dm,dt);new=rebuild(y+dt*dh,mass,nextL)
    restore_frozen(new,nextL,y,L,middle,mass)
    diagnostics=dict(stage_mass_flux=max(d['mass_flux'],dm['mass_flux']),generator_velocity_defect=max(d['generator_velocity_defect'],dm['generator_velocity_defect']))
    if return_record:
        keys=('cov_linearization_factors','cov_A','cov_K','cov_Qfactor')
        diagnostics['hydro_record']=dict(t=t,dt=dt,y=y.copy(),L=[a.copy() for a in L],
            d0=dy.copy(),dm=dh.copy(),diag={k:dm[k] for k in keys},
            frozen_first=first.copy(),frozen_middle=middle.copy())
    return new,nextL,diagnostics


def step(t,y,mass,L,dt,acceleration,frozen=None,return_record=False):
    if not dt>0:raise ValueError('positive split step required')
    kicked,halfL,g0=apply(t,y,mass,L,acceleration,dt/2,mask_at(frozen,t,y,mass))
    drift,driftL,d=hydro_step(t,kicked,mass,halfL,dt,frozen,return_record)
    result,nextL,g1=apply(t+dt,drift,mass,driftL,acceleration,dt/2,mask_at(frozen,t+dt,drift,mass))
    return result,nextL,dict(**d,gravity_work=g0['gravity_work']+g1['gravity_work'],angular_source=g0['angular_source']+g1['angular_source'],max_thermal_kick_change=max(g0['max_thermal_change'],g1['max_thermal_change']),max_kick_completion=max(g0['max_Pi_roundoff_completion'],g1['max_Pi_roundoff_completion']))


def hydro_dense_state(record,mass,s):
    """Actual hydro subflow extension; mask transitions need a separate audit."""
    from coupled_factor_path_audit import dense_state
    first=np.asarray(record['frozen_first'],dtype=bool)
    middle=np.asarray(record['frozen_middle'],dtype=bool)
    if not np.array_equal(first,middle):
        raise ValueError('frozen mask transition needs a piecewise path audit')
    state,factors,derivative=dense_state(record['y'],record['L'],record['d0'],
        record['dm'],record['diag'],mass,record['dt'],s)
    restore_frozen(state,factors,record['y'],record['L'],first,mass)
    derivative[:-4].reshape(len(mass),20)[first]=0.
    return state,factors,derivative


def potential_energy(t,y,mass,potential):
    a=np.asarray(y)[:-4].reshape(len(mass),20);cells,_=base.geometry(a[:,:2],a[:,2],t);total=0.
    for p,m,row in zip(cells,mass,a):
        density=base.decode(p,np.r_[m,row[3:]])['density'];total+=density.weighted_mass@potential(density.points)
    return float(total)
