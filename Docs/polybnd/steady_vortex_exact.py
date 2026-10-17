"""Compact stationary compressible Euler vortex with a smooth quiet exterior.

rho=1, u=Omega q**3 J(x-center), P=1-Omega**2 R**2 q**7/14,
q=max(0,1-r**2/R**2). Positive pressure, zero wall velocity, no gravity.
This reference is for initialization/scoring only.
"""
import numpy as np
from simple_wave_exact import polygon_rule,GAMMA

CENTER=np.array([.5,.5]);RADIUS=.35;POWER=3


def parameters(peak=.5):
    omega=peak*np.sqrt(2*POWER+1)/RADIUS*((2*POWER+1)/(2*POWER))**POWER
    drop=omega**2*RADIUS**2/(2*(2*POWER+1))
    if peak<=0 or drop>=1:raise ValueError('vortex requires positive speed and pressure')
    return omega,drop


def fields(points,t=0.,amplitude=.5,derivatives=False):
    points=np.asarray(points);d=points-CENTER;omega,drop=parameters(amplitude)
    q=np.maximum(0.,1-np.sum(d*d,axis=-1)/RADIUS**2);rotation=np.stack([-d[...,1],d[...,0]],axis=-1)
    f=omega*q**POWER;velocity=f[...,None]*rotation
    rho=np.ones(q.shape);pressure=1-drop*q**(2*POWER+1)
    if not derivatives:return rho,velocity,pressure
    df=-2*POWER*omega/RADIUS**2*q[...,None]**(POWER-1)*d
    J=np.array([[0.,-1.],[1.,0.]])
    gradient=f[...,None,None]*J+rotation[..., :,None]*df[...,None,:]
    pressure_gradient=omega**2*q[...,None]**(2*POWER)*d
    return rho,velocity,pressure,gradient,pressure_gradient


def physical_raw(poly,t,amplitude=.5,order=16):
    points,weights=polygon_rule(poly,order);rho,velocity,pressure=fields(points,t,amplitude)
    dm=weights*rho;Pi=velocity.T@(dm[:,None]*velocity)
    raw=np.r_[dm.sum(),dm@points,(points.T@(dm[:,None]*points)).ravel(),dm@velocity,
        (points.T@(dm[:,None]*velocity)).ravel(),Pi.ravel(),weights@pressure/(GAMMA-1)+.5*np.trace(Pi)]
    return raw,points,weights,rho,velocity,pressure
