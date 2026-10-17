"""Canonical spline gravity and its consistent potential, GM=1 by default.

Matches the spline mass polynomial in Exam/exam.c. Explicit parameters avoid
process-global environment caches. This supplies a source component, not the
coupled Kepler fluid disk benchmark.
"""
import numpy as np


def enclosed_mass(q):
    q=np.asarray(q,dtype=float)
    if np.any(q<0) or not np.all(np.isfinite(q)):raise ValueError('finite nonnegative radius required')
    out=np.ones_like(q);inner=q<=.5;outer=(q>.5)&(q<1)
    x=q[inner];out[inner]=x**3*(32/3-192/5*x*x+32*x**3)
    x=q[outer];out[outer]=-1/15+x**3*(64/3-48*x+192/5*x*x-32/3*x**3)
    return out


def coordinates(points,center,rc,gm):
    x=np.asarray(points,dtype=float)
    if x.ndim<1 or x.shape[-1] not in (2,3) or not np.all(np.isfinite(x)):
        raise ValueError('finite 2D or 3D points required')
    if not np.isfinite(rc) or rc<=0 or not np.isfinite(gm) or gm<=0:raise ValueError('positive gravity parameters required')
    center=np.zeros(x.shape[-1]) if center is None else np.asarray(center,dtype=float)
    if center.shape!=(x.shape[-1],) or not np.all(np.isfinite(center)):raise ValueError('invalid gravity center')
    delta=x-center;return delta,np.linalg.norm(delta,axis=-1)


def acceleration(points,center=None,rc=.5,gm=1.):
    delta,r=coordinates(points,center,rc,gm);q=r/rc
    coefficient=np.zeros_like(r);inner=q<=.5;outer=~inner
    x=q[inner];coefficient[inner]=-(gm/rc**3)*(32/3-192/5*x*x+32*x**3)
    coefficient[outer]=-gm*enclosed_mass(q[outer])/r[outer]**3
    return coefficient[...,None]*delta


def potential(points,center=None,rc=.5,gm=1.):
    _,r=coordinates(points,center,rc,gm);q=r/rc
    def low(x):return 16/3*x*x-48/5*x**4+32/5*x**5
    def high(x):return 1/(15*x)+32/3*x*x-16*x**3+48/5*x**4-32/15*x**5
    out=np.empty_like(q);inner=q<=.5;middle=(q>.5)&(q<1);outer=q>=1
    out[inner]=low(q[inner])-low(.5)+high(.5)-high(1.)-1.
    out[middle]=high(q[middle])-high(1.)-1.
    out[outer]=-1/q[outer]
    return gm/rc*out
