"""Smooth nonlinear isentropic Euler simple wave before characteristic crossing.

One-dimensional exact solution embedded in a 2D reflecting box. Reference
fields are used only for initialization and scoring by the evolution driver.
"""
import numpy as np

GAMMA=5/3
P0=.6
C0=np.sqrt(GAMMA*P0)
LEFT=.2
RIGHT=.6


def initial_velocity(x,amplitude):
    x=np.asarray(x);s=np.pi*(x-LEFT)/(RIGHT-LEFT);inside=(x>LEFT)&(x<RIGHT)
    value=np.where(inside,amplitude*np.sin(s)**4,0.)
    derivative=np.where(inside,4*amplitude*np.pi/(RIGHT-LEFT)*np.sin(s)**3*np.cos(s),0.)
    return value,derivative


def fields(points,t,amplitude=.02,with_time_derivatives=False):
    points=np.asarray(points);x=points[...,0];xi=x-C0*t;k=(GAMMA+1)/2
    # This conservative bound ensures a globally invertible characteristic map.
    derivative_bound=4*abs(amplitude)*np.pi/(RIGHT-LEFT)
    if t<0 or 1-k*t*derivative_bound<=0:raise ValueError('simple-wave reference may have crossed characteristics')
    if RIGHT+(C0+k*max(amplitude,0.))*t>=1.:raise ValueError('simple wave reached reflecting outer wall')
    for _ in range(12):
        value,derivative=initial_velocity(xi,amplitude)
        residual=xi+(C0+k*value)*t-x
        xi=xi-residual/(1+k*derivative*t)
    value,_=initial_velocity(xi,amplitude)
    residual=xi+(C0+k*value)*t-x
    if np.max(abs(residual))>2e-13:raise RuntimeError('simple-wave characteristic inverse failed')
    c=C0+(GAMMA-1)*value/2
    rho=(c/C0)**(2/(GAMMA-1));pressure=P0*rho**GAMMA
    velocity=np.zeros(points.shape);velocity[...,0]=value
    if with_time_derivatives:
        _,slope=initial_velocity(xi,amplitude);ux=slope/(1+k*t*slope)
        ut=-(C0+k*value)*ux;rho_t=rho*ut/c;pressure_t=GAMMA*pressure*ut/c
        velocity_t=np.zeros_like(velocity);velocity_t[...,0]=ut
        return rho,velocity,pressure,rho_t,velocity_t,pressure_t
    return rho,velocity,pressure


def polygon_rule(poly,order=12):
    poly=np.asarray(poly);center=poly.mean(axis=0);g,w=np.polynomial.legendre.leggauss(order)
    g=(g+1)/2;w=w/2;points=[];weights=[]
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        u=a-center;v=b-center;det=abs(u[0]*v[1]-u[1]*v[0])
        points.append(center+g[:,None,None]*u+(1-g[:,None,None])*g[None,:,None]*v)
        weights.append(det*(1-g[:,None])*w[:,None]*w[None,:])
    return np.concatenate([v.reshape(-1,2) for v in points]),np.concatenate([v.ravel() for v in weights])


def physical_raw(poly,t,amplitude,order=12):
    from weighted_mass_patch import clip_halfplane
    # Split exactly at the two compact-pulse fronts before smooth quadrature.
    left=LEFT+C0*t;right=RIGHT+C0*t
    pieces=[clip_halfplane(np.asarray(poly),[1.,0.],left),
        clip_halfplane(clip_halfplane(np.asarray(poly),[-1.,0.],-left),[1.,0.],right),
        clip_halfplane(np.asarray(poly),[-1.,0.],-right)]
    rules=[polygon_rule(v,order) for v in pieces if len(v)>=3]
    points=np.concatenate([v[0] for v in rules]);weights=np.concatenate([v[1] for v in rules])
    rho,velocity,pressure=fields(points,t,amplitude)
    dm=weights*rho;Pi=velocity.T@(dm[:,None]*velocity)
    raw=np.r_[dm.sum(),dm@points,(points.T@(dm[:,None]*points)).ravel(),dm@velocity,
        (points.T@(dm[:,None]*velocity)).ravel(),Pi.ravel(),weights@pressure/(GAMMA-1)+.5*np.trace(Pi)]
    return raw,points,weights,rho,velocity,pressure
