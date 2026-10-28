"""Two-body off-center torque/heat ledger; not a fluid or Kepler solver.

Fixed inertias, evolving centers, common contact point. Optional tangential
slip friction is an explicit diagnostic constitutive law, not adopted physics.
"""
import json
import numpy as np

m=np.array([1.,1.3]); inertia=np.array([.12,.19])
offset=np.array([.12,.23])
def cross(a,b):return a[...,0]*b[...,1]-a[...,1]*b[...,0]
def turn(a):return np.stack((-a[...,1],a[...,0]),axis=-1)

def rates(y,gamma):
    x=y[:4].reshape(2,2); v=y[4:8].reshape(2,2); spin=y[8:10]/inertia
    c=x.mean(axis=0)+offset; arm=c-x
    direction=x[1]-x[0]; normal=direction/np.linalg.norm(direction)
    tangent=turn(normal)
    contact=v+spin[:,None]*turn(arm)
    slip=float((contact[0]-contact[1])@tangent)
    force=-gamma*slip*tangent
    forces=np.array([force,-force]); torque=cross(arm,forces)
    dy=np.zeros_like(y);dy[:4]=v.ravel();dy[4:8]=(forces/m[:,None]).ravel()
    dy[8:10]=torque;dy[10:12]=.5*gamma*slip**2
    mechanical=float(np.sum(forces*v)+np.dot(spin,torque))
    residual=mechanical+float(np.sum(dy[10:12]))
    lrate=float(np.sum(cross(x,forces))+np.sum(torque))
    return dy,dict(slip=slip,power_residual=residual,angular_rate_residual=lrate,
                   rotational_power=float(spin@torque),heat_rate=gamma*slip**2)

def quantities(y):
    x=y[:4].reshape(2,2);v=y[4:8].reshape(2,2)
    translation=float(.5*np.sum(m[:,None]*v*v))
    rotation=float(.5*np.sum(y[8:10]**2/inertia))
    heat=float(sum(y[10:12]))
    return dict(translation=translation,rotation=rotation,heat=heat,
                energy=translation+rotation+heat,
                angular_momentum=float(sum(cross(x,m[:,None]*v))+sum(y[8:10])),
                momentum=np.sum(m[:,None]*v,axis=0).tolist())

initial=np.array([-.5,0,.5,0, 0,.1,0,-.1, .12,-.076, 0,0.])
# An instantaneous rigid rotation has zero contact slip, including off-axis arms.
rigid=initial.copy();omega=.7
rigid[4:8]=(omega*turn(rigid[:4].reshape(2,2))).ravel()
rigid[8:10]=inertia*omega
_,rigid_diag=rates(rigid,.5)
assert abs(rigid_diag['slip'])<1e-14

results=[];finals=[]
for gamma,dt in [(0.,.01),(.5,.02),(.5,.01),(.5,.005)]:
    y=initial.copy();q0=quantities(y);max_power=max_lrate=0.;history=[]
    for step in range(round(10/dt)):
        stages=[]
        k1,d=rates(y,gamma);stages.append(d)
        k2,d=rates(y+.5*dt*k1,gamma);stages.append(d)
        k3,d=rates(y+.5*dt*k2,gamma);stages.append(d)
        k4,d=rates(y+dt*k3,gamma);stages.append(d)
        y+=dt*(k1+2*k2+2*k3+k4)/6
        max_power=max(max_power,*[abs(z['power_residual']) for z in stages])
        max_lrate=max(max_lrate,*[abs(z['angular_rate_residual']) for z in stages])
        assert np.isfinite(y).all()
        if (step+1)%round(1/dt)==0:
            history.append(dict(time=(step+1)*dt,**quantities(y),**rates(y,gamma)[1]))
    q=quantities(y)
    results.append(dict(gamma=gamma,dt=dt,initial=q0,final=q,history=history,
                        energy_error=q['energy']-q0['energy'],
                        angular_error=q['angular_momentum']-q0['angular_momentum'],
                        max_stage_power_residual=max_power,max_stage_angular_rate_residual=max_lrate))
    if gamma:finals.append(y)
    assert max_power<1e-13 and max_lrate<1e-13
    assert abs(q['energy']-q0['energy'])<1e-7
    assert abs(q['angular_momentum']-q0['angular_momentum'])<1e-7
errors=[float(np.linalg.norm(finals[i]-finals[i+1])) for i in range(2)]
assert errors[1]<errors[0]
print(json.dumps(dict(rigid_rotation=rigid_diag,cases=results,
    rk4_successive_differences=errors,ratio=errors[0]/errors[1],
    scope=__doc__+' Relaxation conserves total angular momentum, not each spin. '
    'Rotational energy need not decrease monotonically. Fixed I only; evolving '
    'cells require inertia work and full moment closure. No production force changed.'),indent=2))
