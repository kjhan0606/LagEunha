"""Gravity kick on existing raw moments using positive velocity quadrature.

Keeps original m,q,Q. Primary momentum and mixed moments receive quadrature
force increments. QR forms the new residual factor; no eigenvalue clipping.
Any mismatch between quadrature and existing moments is checked explicitly.
"""
import numpy as np
from moment_gravity_kick import moments


def kick(raw,points,mass_weights,velocity_samples,acceleration,dt,tolerance=1e-10):
    x=np.asarray(points);dm=np.asarray(mass_weights);v=np.asarray(velocity_samples)
    m=raw['m'];q=np.asarray(raw['q']);Q=np.asarray(raw['Q']);p=np.asarray(raw['p']);D=np.asarray(raw['D']);Pi=np.asarray(raw['Pi']);E=raw['E']
    internal=E-.5*np.trace(Pi);represented=moments(x,dm,v,internal)
    scale=max(m,np.linalg.norm(q),np.linalg.norm(Q),np.linalg.norm(p),np.linalg.norm(D),np.linalg.norm(Pi),abs(E))
    mismatch=max(float(np.max(abs(np.asarray(represented[k])-raw[k]))) for k in ('m','q','Q','p','D','Pi'))
    if mismatch>tolerance*scale:raise ValueError('gravity quadrature does not represent the stored moments')
    g=np.asarray(acceleration(x))
    if g.shape!=x.shape or not np.all(np.isfinite(g)) or not np.isfinite(dt):raise ValueError('invalid gravity acceleration or duration')
    force=dm@g;mixed=x.T@(dm[:,None]*g);mean=v.mean(axis=1)
    cross=mean.T@(dm[:,None]*g);squared=g.T@(dm[:,None]*g)
    dPi=dt*(cross+cross.T)+dt*dt*squared;work=.5*np.trace(dPi)
    pn=p+dt*force;Dn=D+dt*mixed;vn=v+dt*g[:,None,:]
    c=q/m;S=Q-np.outer(q,q)/m;C=Dn-np.outer(q,pn)/m;B=np.linalg.solve(S,C).T
    residual=vn-pn/m-((x-c)@B.T)[:,None,:]
    weighted=(np.sqrt(dm/v.shape[1])[:,None,None]*residual).reshape(-1,x.shape[1]);_,R=np.linalg.qr(weighted,mode='reduced');L=R.T
    Pin=np.outer(pn,pn)/m+B@S@B.T+L@L.T
    completion=float(np.max(abs(Pin-(Pi+dPi))))
    if completion>tolerance*scale:raise RuntimeError('gravity positive-factor completion exceeds moment tolerance')
    En=E+work;thermal_change=float((En-.5*np.trace(Pin))-internal)
    if En-.5*np.trace(Pin)<=0:raise RuntimeError('nonpositive thermal energy after gravity kick')
    out=dict(m=m,q=q.copy(),Q=Q.copy(),p=pn,D=Dn,Pi=Pin,E=En)
    return out,L,dict(quadrature_moment_mismatch=mismatch,Pi_roundoff_completion=completion,thermal_change=thermal_change,gravity_work=float(work),angular_source=dt*(mixed-mixed.T))
