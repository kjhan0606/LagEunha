"""Integrate only material sites and mass-constrained weights on a given moment path.

The raw RK2 moment path and positive covariance factors are unchanged. This is
an alternative second-order geometry update, not generator steering. Independent
mass/path audits use the derivative of the saved DOP853 interpolant, not a new
constraint solve substituted into the audit.
"""
import os
import numpy as np
from scipy.integrate import solve_ivp
from coupled_factor_path_audit import dense_state
from dense_path_derivative import segment_derivative
from kinetic_factor_rhs import ALPHA
from vector_kinetic_root import solve_arrays
import moment_ale_evolve as model


def rates(t,state,factors,mass,root_tolerance=2e-13):
    mass=np.asarray(mass);N=len(mass);a=np.asarray(state)[:-4].reshape(N,20)
    x=np.asarray(a[:,:2],dtype=float);weights=np.asarray(a[:,2],dtype=float)
    cells,edges=model.geometry(x,weights,t)
    states=[model.decode(p,np.r_[m,row[3:]]) for p,m,row in zip(cells,mass,a)]
    quadratic=os.environ.get('LAGFORCE_QUADRATIC_VELOCITY','0')=='1'
    if quadratic:
        from quadratic_velocity_reconstruction import prepare,velocity_samples
        prepare(states,cells,edges,mass,factors)
    gen=np.asarray([model.velocity(s,p) for s,p in zip(states,x)])
    g,w=np.polynomial.legendre.leggauss(4);g=(g+1)/2;w=w/2
    from kinetic_speed_smoothing import parameters,face_speed
    smoothing,smoothing_length=parameters()
    rows=[]
    for i,j,start,end,n,length in edges:
        points=start+g[:,None]*(end-start)
        s=states[i];vl=s['v']+(points-s['c'])@s['B'].T
        samples_l=vl[:,None,:]+(ALPHA@factors[i].T/np.sqrt(mass[i]))[None,:,:]
        if quadratic:samples_l,_=velocity_samples(s,factors[i],mass[i],points,ALPHA)
        rl=s['density'].value(points)
        if j>=0:
            s=states[j];vr=s['v']+(points-s['c'])@s['B'].T
            samples_r=vr[:,None,:]+(ALPHA@factors[j].T/np.sqrt(mass[j]))[None,:,:]
            if quadratic:samples_r,_=velocity_samples(s,factors[j],mass[j],points,ALPHA)
            rr=s['density'].value(points);delta=x[j]-x[i];distance=np.linalg.norm(delta)
            base=(delta@gen[j]-(points-x[i])@(gen[j]-gen[i]))/distance
        else:
            base=np.asarray([model.boundary_velocity(p)@n for p in points]);distance=0.;rr=rl
            samples_r=samples_l-2*(samples_l@n-base[:,None])[:,:,None]*n
        rows.append((i,j,length*w,distance,base,samples_l@n,samples_r@n,rl,rr,face_speed(states,i,j,smoothing,smoothing_length)))
    data=dict(i=np.repeat([r[0] for r in rows],4),j=np.repeat([r[1] for r in rows],4),
        scale=np.concatenate([r[2] for r in rows]),distance=np.repeat([r[3] for r in rows],4),
        base=np.concatenate([r[4] for r in rows]),vl=np.concatenate([r[5] for r in rows]),vr=np.concatenate([r[6] for r in rows]),
        rl=np.concatenate([r[7] for r in rows]),rr=np.concatenate([r[8] for r in rows]))
    if smoothing:data['epsilon']=np.repeat([r[9] for r in rows],4)
    if root_tolerance==2e-13:z,_,norm,iterations=solve_arrays(data,N)
    else:z,_,norm,iterations=solve_arrays(data,N,tolerance=root_tolerance)
    return np.column_stack([gen,z]),dict(mass_flux=norm,root_iterations=iterations)


def raw_path(record,mass,s):
    if record.get('dense_moment_arithmetic')=='longdouble':
        from extended_dense_moments import dense_state as extended
        return extended(record['y'],record['L'],record['d0'],record['dm'],record['diag'],mass,record['dt'],s)
    return dense_state(record['y'],record['L'],record['d0'],record['dm'],record['diag'],mass,record['dt'],s)


def evaluate_polynomial(segment,t):
    x=(t-segment['t_old'])/segment['h'];value=np.zeros_like(np.asarray(segment['y_old']),dtype=np.result_type(t,float))
    for i,f in enumerate(reversed(np.asarray(segment['F']))):
        value+=f
        value*=x if i%2==0 else 1-x
    return value+segment['y_old']


def polynomial_derivative(segment,t):
    # The helper differentiates the saved polynomial, independently of rates().
    from types import SimpleNamespace
    return segment_derivative(SimpleNamespace(t_old=segment['t_old'],h=segment['h'],y_old=np.asarray(segment['y_old']),F=np.asarray(segment['F'])),t)


def lifted_state(record,mass,s):
    state,factors,derivative=raw_path(record,mass,s)
    segments=record['geometry_lift']['segments'];t=record['dt']*s
    ends=np.asarray([v['t'] for v in segments]);index=int(np.clip(np.searchsorted(ends,t,side='left'),0,len(segments)-1))
    segment=segments[index];a=state[:-4].reshape(-1,20);da=derivative[:-4].reshape(-1,20)
    a[:,:3]=evaluate_polynomial(segment,t).reshape(-1,3)
    da[:,:3]=polynomial_derivative(segment,t).reshape(-1,3)
    return state,factors,derivative


def integrate(record,mass,rtol=1e-12,atol=1e-14,root_tolerance=2e-13):
    if 'geometry_lift' in record:raise ValueError('geometry lift already present')
    initial=np.asarray(record['y'])[:-4].reshape(-1,20)[:,:3].ravel().copy();dt=record['dt'];peak=0.;evaluations=0
    def fun(t,coordinates):
        nonlocal peak,evaluations
        state,L,_=raw_path(record,mass,t/dt);state[:-4].reshape(-1,20)[:,:3]=coordinates.reshape(-1,3)
        derivative,diag=rates(record['t']+t,state,L,mass,root_tolerance=root_tolerance);peak=max(peak,diag['mass_flux']);evaluations+=1
        if evaluations%128==0:print('geometry lift',record['t'],dt,'evaluations',evaluations,'phase',t/dt,'mass root',diag['mass_flux'],flush=True)
        return derivative.ravel()
    solution=solve_ivp(fun,(0.,dt),initial,method='DOP853',rtol=rtol,atol=atol,dense_output=True)
    if not solution.success:raise RuntimeError('material/mass geometry lift failed: '+solution.message)
    segments=[dict(t_old=v.t_old,t=v.t,h=v.h,y_old=v.y_old.copy(),F=v.F.copy()) for v in solution.sol.interpolants]
    return dict(segments=segments,rtol=rtol,atol=atol,mass_root_tolerance=root_tolerance,evaluations=solution.nfev,maximum_mass_flux=peak,
        note='Geometry ODE from numerical mean fluid velocity and the original kinetic mass root. The stored polynomial derivative is independently audited.')
