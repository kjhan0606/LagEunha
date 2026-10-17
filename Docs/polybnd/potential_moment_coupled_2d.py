"""Experimental coupled evolution with an additional density potential moment.

U=int rho Phi evolves from the same mass flux and physical gravity work.
All fluid stages and gravity kicks use the corresponding enriched density.
No stored mass or gas energy is corrected to force a conservation result.
"""
import os
import numpy as np
import moment_ale_evolve as base
from potential_enriched_density import PotentialDensity
from kinetic_factor_rhs import rhs,ALPHA
from gravity_source_2d import apply
from gravity_coupled_2d import mask_at,restore_frozen
from factor_moment_evolve import rebuild,advance_factors_sylvester
from smooth_mass_path import evaluate_smooth_path
from quadratic_velocity_reconstruction import prepare,velocity_samples


def provider(U,potential,quadrature,frozen):
    cache={}
    def get(i,poly,raw):
        key=(i,np.asarray(poly).tobytes(),np.asarray(raw[:7]).tobytes())
        if key not in cache:
            m,q,Q=raw[0],raw[1:3],raw[3:7].reshape(2,2)
            cache[key]=base.Density(poly,m,q,Q) if frozen[i] else PotentialDensity(poly,m,q,Q,U[i],potential,quadrature)
        return cache[key]
    return get


def initial_potential(t,y,mass,potential,quadrature):
    a=y[:-4].reshape(-1,20);cells,_=base.geometry(a[:,:2],a[:,2],t);values=[]
    for p,m,row in zip(cells,mass,a):
        d=base.decode(p,np.r_[m,row[3:]])['density'];x,w=quadrature(p)
        values.append(float((w*d.value(x))@potential(x)))
    return np.asarray(values)


def stage(t,y,mass,L,U,potential,acceleration,quadrature,mask):
    get=provider(U,potential,quadrature,mask)
    dy,diag=rhs(t,y,mass,covariance_factors=L,frozen_mask=mask,density_provider=get)
    rates=dy[:-4].reshape(-1,20)
    flux,exchange=evaluate_smooth_path(t,y,L,mass,(rates[:,:2],rates[:,2]),
        frozen_mask=mask,potential=potential,density_provider=get)
    if np.max(abs(flux[~mask]),initial=0.)>=2e-13:raise RuntimeError('enriched stage mass root mismatch')
    a=y[:-4].reshape(-1,20);cells,edges=base.geometry(a[:,:2],a[:,2],t)
    raw=np.column_stack([mass,a[:,3:]])
    states=[base.decode(p,r,density=get(i,p,r)) for i,(p,r) in enumerate(zip(cells,raw))]
    if os.environ.get('LAGFORCE_QUADRATIC_VELOCITY')=='1':prepare(states,cells,edges,mass,L)
    work=np.zeros(len(mass))
    for i,s in enumerate(states):
        if mask[i]:continue
        d=s['density'];v,_=velocity_samples(s,L[i],mass[i],d.points,ALPHA)
        work[i]=d.weighted_mass@np.sum(v.mean(1)*acceleration(d.points),axis=1)
    du=exchange['potential_mass_rates']-work;du[mask]=0.
    diag['potential_boundary_rate']=exchange['active_potential_mass_exchange']
    diag['gravity_work_rate']=float(work.sum())
    return dy,du,diag


def step(t,y,mass,L,U,dt,potential,acceleration,quadrature,frozen=None,return_record=False,gravity_scheme='finite'):
    if dt<=0:raise ValueError('positive timestep required')
    if gravity_scheme=='finite':kick=apply
    elif gravity_scheme=='midpoint':
        from midpoint_gravity_source import apply as kick
    else:raise ValueError('unknown gravity integration scheme')
    mask=mask_at(frozen,t,y,mass);oldU=np.asarray(U).copy()
    kicked,halfL,g0=kick(t,y,mass,L,acceleration,dt/2,mask,provider(U,potential,quadrature,mask))
    dy,du,d0=stage(t,kicked,mass,halfL,U,potential,acceleration,quadrature,mask)
    middleL=advance_factors_sylvester(halfL,d0,dt/2)
    middle=rebuild(kicked+dt*dy/2,mass,middleL);middleU=U+dt*du/2
    restore_frozen(middle,middleL,kicked,halfL,mask,mass);middleU[mask]=oldU[mask]
    if not np.array_equal(mask_at(frozen,t+dt/2,middle,mass),mask):raise ValueError('frozen mask transition')
    dh,duh,dm=stage(t+dt/2,middle,mass,middleL,middleU,potential,acceleration,quadrature,mask)
    driftL=advance_factors_sylvester(halfL,dm,dt);drift=rebuild(kicked+dt*dh,mass,driftL);newU=U+dt*duh
    restore_frozen(drift,driftL,kicked,halfL,mask,mass);newU[mask]=oldU[mask]
    if not np.array_equal(mask_at(frozen,t+dt,drift,mass),mask):raise ValueError('frozen mask transition')
    out,newL,g1=kick(t+dt,drift,mass,driftL,acceleration,dt/2,mask,provider(newU,potential,quadrature,mask))
    diag=dict(stage_mass_flux=max(d0['mass_flux'],dm['mass_flux']),
        generator_velocity_defect=max(d0['generator_velocity_defect'],dm['generator_velocity_defect']),
        potential_boundary_exchange=dt*dm['potential_boundary_rate'],
        gravity_kick_work=g0['gravity_work']+g1['gravity_work'],
        hydro_potential_work=dt*dm['gravity_work_rate'],gravity_scheme=gravity_scheme)
    if return_record:
        keys=('cov_linearization_factors','cov_A','cov_K','cov_Qfactor')
        diag['hydro_record']=dict(t=t,dt=dt,y=kicked,L=halfL,d0=dy,dm=dh,
            diag={k:dm[k] for k in keys},frozen_first=mask,frozen_middle=mask,
            U=oldU,du0=du,dum=duh)
    return out,newL,newU,diag


def dense_state(record,mass,s):
    from gravity_coupled_2d import hydro_dense_state
    if not 0<=s<=1:raise ValueError('dense parameter outside step')
    y,L,dy=hydro_dense_state(record,mass,s)
    U=np.asarray(record['U']);d0=np.asarray(record['du0']);dm=np.asarray(record['dum'])
    potential=U+record['dt']*(s*d0+s*s*(dm-d0))
    rate=d0+2*s*(dm-d0)
    mask=np.asarray(record['frozen_first'],dtype=bool)
    potential[mask]=U[mask];rate[mask]=0.
    return y,L,dy,potential,rate
