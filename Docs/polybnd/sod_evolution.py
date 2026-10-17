"""Sod benchmark adapter for the existing fixed-mass positive moment RHS."""
import os
import numpy as np
from adaptive_moment_batch import setup
from factor_moment_evolve import rebuild
import moment_ale_evolve as model
import kinetic_factor_rhs as kinetic
from sod_exact import physical_raw


def initial(n,amplitude):
    if n%2:raise ValueError('aligned Sod mesh requires even n')
    _,state=setup('wave',n,0.);a=state[:-4].reshape(-1,20)
    xx,yy=np.meshgrid((np.arange(n)+.5)/n,(np.arange(n)+.5)/n)
    a[:,:2]=np.column_stack([xx.ravel(),yy.ravel()])
    cells,_=model.geometry(a[:,:2],a[:,2],0.);mass=[]
    for row,poly in zip(a,cells):
        raw,*_=physical_raw(poly,0.,amplitude);mass.append(raw[0]);row[3:]=raw[1:]
    mass=np.asarray(mass);L=[np.zeros((2,2)) for _ in mass];new=rebuild(state,mass,L)
    return mass,new,L,float(np.max(abs(new-state)))


def score(t,state,mass,amplitude,factors=None):
    a=np.asarray(state)[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],t)
    states=[model.decode(poly,np.r_[m,row[3:]]) for m,row,poly in zip(mass,a,cells)]
    quadratic=os.environ.get('LAGFORCE_QUADRATIC_VELOCITY','0')=='1'
    if quadratic:
        from quadratic_velocity_reconstruction import prepare,correction
        prepare(states,cells,edges,mass,factors)
    centers,grad=kinetic.pressure_gradients(states,edges)
    from quadratic_pressure_reconstruction import enabled,prepare as prepare_pressure,value
    pressure_poly=prepare_pressure(states,cells)[0] if enabled() else None
    errors={k:0. for k in ['pressure_cell','velocity_cell','density_cell','pressure_field','velocity_field','density_field','energy_cell']}
    refmass=quad=0.;minrho=minpressure=float('inf');entropy=0.
    for i,(m,row,poly,s) in enumerate(zip(mass,a,cells,states)):
        raw,points,w,rho,v,P=physical_raw(poly,t,amplitude,16);fine,*_=physical_raw(poly,t,amplitude,24)
        quad+=float(np.max(abs(fine-raw)));refmass+=raw[0];area=w.sum()
        vn=s['v']+(points-s['c'])@s['B'].T
        if quadratic:vn+=correction(s,points)
        pn=value(pressure_poly[i],points) if pressure_poly is not None else s['P']+(points-centers[i])@grad[i]
        rn=s['density'].value(points)
        errors['density_cell']+=abs(m-raw[0]);errors['velocity_cell']+=area*np.linalg.norm(s['v']-raw[7:9]/raw[0])
        errors['pressure_cell']+=abs(area*s['P']-w@P);errors['energy_cell']+=abs(row[19]-raw[17])
        errors['density_field']+=w@abs(rn-rho);errors['velocity_field']+=w@np.linalg.norm(vn-v,axis=1);errors['pressure_field']+=w@abs(pn-P)
        minrho=min(minrho,float(rn.min()));minpressure=min(minpressure,float(pn.min()))
        entropy+=float((w*rn)@np.log(pn/rn**model.GAMMA))
    return dict(t=t,amplitude=amplitude,errors_over_amplitude={k:float(v) for k,v in errors.items()},
                reference_quadrature_difference=quad,reference_global_mass_relative=float(abs(refmass-mass.sum())/mass.sum()),
                density_min=minrho,pressure_min=minpressure,integrated_entropy=entropy,
                reference_valid=bool(quad<1e-5 and abs(refmass-mass.sum())/mass.sum()<1e-10))
