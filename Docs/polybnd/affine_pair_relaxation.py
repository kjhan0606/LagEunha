"""Experimental fixed-geometry pair relaxation of full affine contact velocity.

Not a complete Euler closure. S is the mass-weighted central spatial moment;
R is residual velocity covariance integral and is retained, not thermalized.
"""
import numpy as np

def kinetic(a):
    return .5*a['m']*np.dot(a['v'],a['v'])+.5*np.trace(a['G']@a['S']@a['G'].T+a['R'])

def angular(a):
    C=a['S']@a['G'].T
    return a['m']*(a['c'][0]*a['v'][1]-a['c'][1]*a['v'][0])+C[0,1]-C[1,0]

def relax(a,b,point,tangent,gamma,dt):
    if gamma<0 or dt<0 or not np.isfinite(gamma+dt):raise ValueError('invalid damping or duration')
    if not np.isfinite(point).all() or not np.isfinite(tangent).all():raise ValueError('nonfinite geometry')
    if abs(np.dot(tangent,tangent)-1)>1e-12:raise ValueError('unit tangent required')
    old=[a,b];arms=[];duals=[];mobilities=[];contact=[]
    for s in old:
        if not (s['m']>0 and s['U']>0):raise ValueError('positive mass and thermal energy required')
        if not all(np.isfinite(s[k]).all() for k in ['c','v','S','G','R']):raise ValueError('nonfinite state')
        if np.max(abs(s['S']-s['S'].T))>1e-12 or np.linalg.eigvalsh(s['S']).min()<=0:raise ValueError('positive symmetric S required')
        arm=point-s['c'];dual=np.linalg.solve(s['S'],arm)
        arms.append(arm);duals.append(dual);mobilities.append(1/s['m']+arm@dual)
        contact.append(s['v']+s['G']@arm)
    H=sum(mobilities);slip=float((contact[0]-contact[1])@tangent)
    fraction=-np.expm1(-gamma*H*dt)
    impulse=-fraction*slip/H*tangent
    heat=.5*slip**2/H*fraction*(2-fraction)
    new=[]
    for k,s in enumerate(old):
        j=impulse if k==0 else -impulse
        z={key:value.copy() if isinstance(value,np.ndarray) else value for key,value in s.items()}
        z['v']=s['v']+j/s['m'];z['G']=s['G']+np.outer(j,duals[k]);z['U']=s['U']+.5*heat
        new.append(z)
    return *new,dict(slip_before=slip,slip_after_expected=(1-fraction)*slip,
                     generated_heat=heat,impulse=impulse.tolist(),mobility=H)
