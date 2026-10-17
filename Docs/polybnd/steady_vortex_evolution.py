"""Stationary 2D Euler vortex benchmark of the numerical moment closure."""
import os
import numpy as np
from adaptive_moment_batch import setup
from factor_moment_evolve import rebuild
from steady_vortex_exact import physical_raw,parameters
from simple_wave_evolution import install_domain
import moment_ale_evolve as model
import kinetic_factor_rhs as kinetic


def initial(n,amplitude):
    _,skeleton=setup('wave',n,0.);a=skeleton[:-4].reshape(-1,20)
    cells,_=model.geometry(a[:,:2],a[:,2],0.);masses=[];factors=[]
    for i,poly in enumerate(cells):
        raw,points,weights,rho,velocity,pressure=physical_raw(poly,0.,amplitude,16)
        m=raw[0];c=raw[1:3]/m;mean=raw[7:9]/m
        S=raw[3:7].reshape(2,2)-np.outer(raw[1:3],raw[1:3])/m
        C=raw[9:13].reshape(2,2)-np.outer(raw[1:3],raw[7:9])/m;B=np.linalg.solve(S,C).T
        residual=velocity-mean-(points-c)@B.T
        _,R=np.linalg.qr(np.sqrt(weights*rho)[:,None]*residual,mode='reduced')
        factors.append(R.T);masses.append(m);a[i,3:]=raw[1:]
    mass=np.asarray(masses);state=rebuild(skeleton,mass,factors)
    return mass,state,factors,float(np.max(abs(state-skeleton)))


def score(t,state,mass,amplitude,factors=None):
    a=np.asarray(state)[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],t)
    states=[model.decode(poly,np.r_[m,row[3:]]) for m,row,poly in zip(mass,a,cells)]
    quadratic=os.environ.get('LAGFORCE_QUADRATIC_VELOCITY','0')=='1'
    if quadratic:
        if factors is None:raise ValueError('vortex quadratic score needs primary factors')
        from quadratic_velocity_reconstruction import prepare,correction
        prepare(states,cells,edges,mass,factors)
    centers,grad=kinetic.pressure_gradients(states,edges);_,drop=parameters(amplitude)
    from quadratic_pressure_reconstruction import enabled as quadratic_pressure_enabled
    pressure_polynomials=None
    if quadratic_pressure_enabled():
        from quadratic_pressure_reconstruction import prepare as prepare_pressure,value as pressure_value
        pressure_polynomials,_=prepare_pressure(states,cells)
    errors={k:0. for k in ['pressure_cell','velocity_cell','density_cell','pressure_field','velocity_field','density_field','energy_cell']}
    exact_mass=reference_quad=0.;minrho=minpressure=float('inf')
    for i,(m,row,poly,s) in enumerate(zip(mass,a,cells,states)):
        raw,points,weights,rho,v,P=physical_raw(poly,t,amplitude,16);fine,*_=physical_raw(poly,t,amplitude,24)
        reference_quad+=float(np.max(abs(fine-raw)))/min(amplitude,drop)
        area=weights.sum();exact_mass+=raw[0]
        errors['pressure_cell']+=abs(area*s['P']-weights@P)/drop
        errors['velocity_cell']+=area*np.linalg.norm(s['v']-raw[7:9]/raw[0])/amplitude
        errors['density_cell']+=abs(m-raw[0])
        errors['energy_cell']+=abs(row[19]-raw[17])/drop
        vn=s['v']+(points-s['c'])@s['B'].T
        if quadratic:vn+=correction(s,points)
        pn=s['P']+(points-centers[i])@grad[i];rn=s['density'].value(points)
        if pressure_polynomials is not None:pn=pressure_value(pressure_polynomials[i],points)
        errors['velocity_field']+=weights@np.linalg.norm(vn-v,axis=1)/amplitude
        errors['pressure_field']+=weights@abs(pn-P)/drop
        errors['density_field']+=weights@abs(rn-rho)
        minrho=min(minrho,float(rn.min()));minpressure=min(minpressure,float(pn.min()))
    return dict(t=t,amplitude=amplitude,pressure_scale=drop,density_scale=1.,errors_over_amplitude={k:float(v) for k,v in errors.items()},
        density_min=minrho,pressure_min=minpressure,reference_quadrature_difference=reference_quad,
        reference_global_mass_relative=float(abs(exact_mass-mass.sum())/mass.sum()),
        reference_valid=bool(reference_quad<1e-5 and abs(exact_mass-mass.sum())/mass.sum()<1e-10),
        note='Velocity errors normalized by peak speed, pressure/energy by pressure drop, density by background1. Stationary exact fields only initialize/score; cells and reference cell integrals move.')
