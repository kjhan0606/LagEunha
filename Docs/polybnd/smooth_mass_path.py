"""Mass-only audit of a numerical path using its continuous moment density."""
import numpy as np
import json
import os
from pathlib import Path
import moment_ale_evolve as model
from kinetic_factor_rhs import ALPHA

_record_key=None
_record_cache=None


def evaluation_backend_name():return 'cell_batched' if os.environ.get('LAGFORCE_VECTOR_SMOOTH_PATH','0')=='1' else 'scalar'


def initialize_worker():
    from adaptive_moment_batch import D0
    model.domain=lambda t:D0.copy()
    model.boundary_velocity=lambda p:np.zeros(2)


def sample_record(record,mass,s):
    from coupled_factor_path_audit import dense_state
    from adaptive_coupled_mass import restore
    r=restore(record)
    if 'geometry_lift' in r:
        from geometry_mass_lift import lifted_state
        state,factors,derivative=lifted_state(r,mass,s)
    elif r.get('dense_moment_arithmetic')=='longdouble':
        from extended_dense_moments import dense_state as extended
        state,factors,derivative=extended(r['y'],r['L'],r['d0'],r['dm'],r['diag'],mass,r['dt'],s)
    else:
        state,factors,derivative=dense_state(r['y'],r['L'],r['d0'],r['dm'],r['diag'],mass,r['dt'],s)
    rates=derivative[:-4].reshape(-1,20)
    flux,diag=evaluate_smooth_path(r['t']+s*r['dt'],state,factors,mass,(rates[:,:2],rates[:,2]))
    return dict(flux=flux,diag=diag)


def sample_record_file(path,expected_digest,mass,s):
    from predictive_evolution_batch import record_digest
    global _record_key,_record_cache
    key=(str(Path(path).resolve()),expected_digest)
    if _record_key!=key:
        record=json.loads(Path(path).read_text())['record']
        if record_digest(record)!=expected_digest:raise ValueError('smooth path record hash mismatch')
        _record_cache=record;_record_key=key
    return sample_record(_record_cache,mass,s)


def evaluate_smooth_path(t,state,factors,mass,geometry_rates,order=4,frozen_mask=None,potential=None):
    a=np.asarray(state)[:-4].reshape(-1,20);mass=np.asarray(mass);N=len(mass)
    x=np.asarray(a[:,:2],dtype=float);weights=np.asarray(a[:,2],dtype=float)
    cells,edges=model.geometry(x,weights,t)
    raw=np.column_stack([mass,a[:,3:]]);states=[model.decode(v,r) for v,r in zip(cells,raw)]
    quadratic=os.environ.get('LAGFORCE_QUADRATIC_VELOCITY','0')=='1'
    closure={}
    if quadratic:
        from quadratic_velocity_reconstruction import prepare,velocity_samples
        closure=prepare(states,cells,edges,mass,factors)
    gen,z=geometry_rates;material=np.asarray([model.velocity(s,x) for s,x in zip(states,x)])
    if frozen_mask is not None:
        frozen=np.asarray(frozen_mask,dtype=bool)
        if frozen.shape!=(N,):raise ValueError('invalid frozen-cell mask')
        if np.any(np.asarray(gen)[frozen]!=0) or np.any(np.asarray(z)[frozen]!=0):
            raise ValueError('frozen geometry has nonzero prescribed rates')
        material[frozen]=0.
    if os.environ.get('LAGFORCE_VECTOR_SMOOTH_PATH','0')=='1' and potential is None:
        from vector_smooth_mass_path import evaluate
        flux=evaluate(states,edges,x,factors,mass,gen,z,order,model.boundary_velocity)
        return flux,dict(generator_defect=float(np.max(abs(gen-material))),
            min_pressure=min(s['P'] for s in states),
            min_factor_R=min(float(np.linalg.eigvalsh(L@L.T).min()/m) for L,m in zip(factors,mass)),**closure)
    from kinetic_speed_smoothing import parameters,face_speed,split
    smoothing,smoothing_length=parameters()
    active=np.ones(N,dtype=bool) if frozen_mask is None else ~np.asarray(frozen_mask,dtype=bool)
    potential_exchange=0.
    g,w=np.polynomial.legendre.leggauss(order);g=(g+1)/2;w=w/2;residual=np.zeros(N)
    for i,j,start,end,n,length in edges:
        points=start+g[:,None]*(end-start)
        left=states[i];vl=left['v']+(points-left['c'])@left['B'].T
        samples_l=vl[:,None,:]+(ALPHA@factors[i].T/np.sqrt(mass[i]))[None,:,:]
        if quadratic:samples_l,_=velocity_samples(left,factors[i],mass[i],points,ALPHA)
        rl=left['density'].value(points)
        if j>=0:
            right=states[j];vr=right['v']+(points-right['c'])@right['B'].T
            samples_r=vr[:,None,:]+(ALPHA@factors[j].T/np.sqrt(mass[j]))[None,:,:]
            if quadratic:samples_r,_=velocity_samples(right,factors[j],mass[j],points,ALPHA)
            rr=right['density'].value(points);delta=x[j]-x[i];distance=np.linalg.norm(delta)
            wn=(delta@gen[j]-(points-x[i])@(gen[j]-gen[i])+.5*(z[i]-z[j]))/distance
        else:
            wn=np.asarray([model.boundary_velocity(p)@n for p in points]);rr=rl
            samples_r=samples_l-2*(samples_l@n-wn[:,None])[:,:,None]*n
        sl=samples_l@n-wn[:,None];sr=samples_r@n-wn[:,None]
        epsilon=face_speed(states,i,j,smoothing,smoothing_length)
        if epsilon:
            lp,_,_,_=split(sl,epsilon);_,rm,_,_=split(sr,epsilon)
            flux=.25*(rl*lp.sum(1)+rr*rm.sum(1))
        else:flux=.25*(rl*np.maximum(sl,0.).sum(1)+rr*np.minimum(sr,0.).sum(1))
        value=float(length*(w@flux));residual[i]+=value
        if j>=0:residual[j]-=value
        if potential is not None:
            sign=-int(active[i])+(int(active[j]) if j>=0 else 0)
            if sign:
                phi=np.asarray(potential(points))
                if phi.shape!=(len(points),) or not np.isfinite(phi).all():raise ValueError('invalid potential samples')
                potential_exchange+=sign*float(length*(w@(flux*phi)))
    diag=dict(generator_defect=float(np.max(abs(gen-material))),
        min_pressure=min(s['P'] for s in states),
        min_factor_R=min(float(np.linalg.eigvalsh(L@L.T).min()/m) for L,m in zip(factors,mass)),**closure)
    if potential is not None:diag['active_potential_mass_exchange']=potential_exchange
    return -residual,diag
