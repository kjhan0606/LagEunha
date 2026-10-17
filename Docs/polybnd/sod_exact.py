"""Sod states with gamma=5/3, matching this moment solver's equation of state.

Reflecting unit square, discontinuity x=.5. Reference is valid before either
outer wave reaches a wall. Used only for initialization and error measurement.
"""
import numpy as np
from scipy.optimize import brentq
from functools import lru_cache
from simple_wave_exact import GAMMA, polygon_rule
from weighted_mass_patch import clip_halfplane


@lru_cache(None)
def star():
    g=GAMMA; cl=np.sqrt(g); cr=np.sqrt(g*.1/.125)
    def wave(p,r,pk,c):
        if p>pk:
            return (p-pk)*np.sqrt(2/((g+1)*r)/(p+(g-1)/(g+1)*pk))
        return 2*c/(g-1)*((p/pk)**((g-1)/(2*g))-1)
    p=brentq(lambda p:wave(p,1.,1.,cl)+wave(p,.125,.1,cr),.1,1.,xtol=1e-14)
    u=.5*(wave(p,.125,.1,cr)-wave(p,1.,1.,cl))
    rl=p**(1/g); ratio=(g-1)/(g+1)
    rr=.125*(p/.1+ratio)/(ratio*p/.1+1)
    shock=cr*np.sqrt((g+1)/(2*g)*p/.1+(g-1)/(2*g))
    return dict(p=p,u=u,rhol=rl,rhor=rr,head=-cl,tail=u-np.sqrt(g*p/rl),shock=shock)


def fronts(t):
    s=star()
    return [.5+t*s[k] for k in ('head','tail','u','shock')]


def fields(points,t,amplitude=1.):
    if amplitude!=1.:raise ValueError('Sod amplitude is fixed at one')
    if t<0 or (t>0 and (fronts(t)[0]<=0 or fronts(t)[-1]>=1)):
        raise ValueError('Sod reference reached reflecting wall')
    x=np.asarray(points)[...,0];rho=np.where(x<.5,1.,.125);p=np.where(x<.5,1.,.1)
    v=np.zeros(np.asarray(points).shape)
    if t==0:return rho,v,p
    s=star();xi=(x-.5)/t;g=GAMMA;cl=np.sqrt(g)
    fan=(xi>=s['head'])&(xi<s['tail'])
    left=(xi>=s['tail'])&(xi<s['u']);right=(xi>=s['u'])&(xi<s['shock'])
    uf=2/(g+1)*(cl+xi[fan]);cf=cl-(g-1)*uf/2
    rho[fan]=(cf/cl)**(2/(g-1));p[fan]=(cf/cl)**(2*g/(g-1));v[...,0][fan]=uf
    rho[left]=s['rhol'];rho[right]=s['rhor'];p[left|right]=s['p'];v[...,0][left|right]=s['u']
    return rho,v,p


def physical_raw(poly,t,amplitude=1.,order=16):
    cuts=[.5] if t==0 else fronts(t)
    remaining=np.asarray(poly);pieces=[]
    for x in cuts:
        part=clip_halfplane(remaining,[1.,0.],x)
        if len(part)>=3:pieces.append(part)
        remaining=clip_halfplane(remaining,[-1.,0.],-x)
        if len(remaining)<3:break
    if len(remaining)>=3:pieces.append(remaining)
    rules=[polygon_rule(p,order) for p in pieces]
    points=np.concatenate([p for p,w in rules]);weights=np.concatenate([w for p,w in rules])
    rho,v,P=fields(points,t,amplitude);dm=weights*rho;Pi=v.T@(dm[:,None]*v)
    raw=np.r_[dm.sum(),dm@points,(points.T@(dm[:,None]*points)).ravel(),dm@v,
              (points.T@(dm[:,None]*v)).ravel(),Pi.ravel(),weights@P/(GAMMA-1)+.5*np.trace(Pi)]
    return raw,points,weights,rho,v,P
