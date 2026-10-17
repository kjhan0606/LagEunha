"""Optional positive smooth flux splitting with a vanishing O(h^2) speed.

This changes the numerical transport flux; it is an explicit research option,
not a mass-root tolerance adjustment. Stored moments and material sites remain
unchanged. Zero coefficient reproduces the original upwind splitting.
"""
import os
import numpy as np
from fixedmass_hllc import GAMMA


def parameters():
    coefficient=float(os.environ.get('LAGFORCE_KINETIC_SMOOTHING','0'))
    length=float(os.environ.get('LAGFORCE_KINETIC_SMOOTHING_LENGTH','1'))
    if not np.isfinite(coefficient) or coefficient<0 or not np.isfinite(length) or length<=0:raise ValueError('invalid kinetic smoothing parameters')
    return coefficient,length


def face_speed(states,i,j,coefficient,length=1.):
    if coefficient==0:return 0.
    si=states[i];sj=states[j] if j>=0 else si
    ai=si['density'].area;aj=sj['density'].area
    ci=np.sqrt(GAMMA*si['P']*ai/float(si['raw'][0]));cj=np.sqrt(GAMMA*sj['P']*aj/float(sj['raw'][0]))
    return float(coefficient*.5*(ci+cj)*.5*(ai+aj)/length**2)


def split(speed,epsilon):
    speed=np.asarray(speed);epsilon=np.broadcast_to(epsilon,speed.shape)
    if not np.any(epsilon):
        plus=np.maximum(speed,0.);minus=np.minimum(speed,0.)
        derivative=np.where(speed>0,1.,np.where(speed<0,0.,.5))
        return plus,minus,derivative,1-derivative
    root=np.hypot(speed,epsilon)
    # Stable tiny positive tails, without subtraction of root and |speed|.
    denominator=root+abs(speed)
    tail=np.divide(.5*epsilon**2,denominator,out=np.zeros_like(root),where=denominator>0)
    plus=np.where(speed>=0,.5*(speed+root),tail)
    minus=-np.where(speed<=0,.5*(-speed+root),tail)
    dp=np.divide(plus,root,out=np.full_like(root,.5),where=root>0)
    dm=np.divide(-minus,root,out=np.full_like(root,.5),where=root>0)
    return plus,minus,dp,dm
