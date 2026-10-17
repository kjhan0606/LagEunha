"""Exact non-vacuum Euler contact pressure/speed for shared face traction.

An explicit alternative to Davis/HLLC traction, never an automatic fallback.
Normal velocity determines the scalar Riemann problem in either 2D or 3D.
"""
import numpy as np
from scipy.optimize import brentq


def star(rho_left,velocity_left,pressure_left,rho_right,velocity_right,pressure_right,normal,gamma=5/3):
    rl,rr,pl,pr=map(float,(rho_left,rho_right,pressure_left,pressure_right))
    if min(rl,rr,pl,pr)<=0 or not all(np.isfinite(v) for v in (rl,rr,pl,pr,gamma)) or gamma<=1:
        raise ValueError('positive finite Euler states required')
    ul=float(np.dot(velocity_left,normal));ur=float(np.dot(velocity_right,normal))
    cl=np.sqrt(gamma*pl/rl);cr=np.sqrt(gamma*pr/rr)
    if not np.isfinite(ul+ur):raise ValueError('finite velocity required')
    if ur-ul>=2*(cl+cr)/(gamma-1):raise ValueError('Euler vacuum requires a separate vacuum model')
    def wave(p,rho,pk,c):
        if p>pk:return (p-pk)*np.sqrt(2/((gamma+1)*rho)/(p+(gamma-1)/(gamma+1)*pk))
        return 2*c/(gamma-1)*((p/pk)**((gamma-1)/(2*gamma))-1)
    scale=max(pl,pr)
    def residual(q):return wave(q*scale,rl,pl,cl)+wave(q*scale,rr,pr,cr)+ur-ul
    upper=1.
    while residual(upper)<0:
        upper*=2
        if not np.isfinite(upper*scale):raise RuntimeError('Euler star pressure bracket overflow')
    q=brentq(residual,0.,upper,xtol=np.finfo(float).tiny,rtol=8*np.finfo(float).eps,maxiter=200)
    p=q*scale
    if not p>0:raise RuntimeError('unresolved positive Euler star pressure')
    speed=.5*(ul+ur+wave(p,rr,pr,cr)-wave(p,rl,pl,cl))
    return p,speed
