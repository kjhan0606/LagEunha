"""Actual shared-face assembly of affine pair relaxation, frozen 2D geometry.

This isolates the relaxation substep; it does not evolve Euler/Kepler geometry.
"""
import json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_moving_cell_spin_20261018')
sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_affine_pair_relaxation_20261018')
from affine_moment_patch import bounded_mesh,integrate,P0,GAMMA
from affine_pair_relaxation import relax,kinetic,angular

def setup(n,curved):
    xx,yy=np.meshgrid((np.arange(n)+.5)/n,(np.arange(n)+.5)/n)
    x=np.column_stack((xx.ravel(),yy.ravel()))+np.random.default_rng(731).uniform(-.15,.15,(n*n,2))/n
    polygons,edges=bounded_mesh(x,np.zeros(n*n))
    velocity=lambda p:np.column_stack((.5*np.sin(2*np.pi*p[:,1]) if curved else .2+1.1*p[:,1],np.zeros(len(p))))
    states=[]
    for poly in polygons:
        s=integrate(poly,velocity);G=np.linalg.solve(s['S'],s['C']).T
        states.append(dict(m=s['area'],c=s['center'],v=s['mean'],S=s['S'],G=G,
                           R=s['T']-G@s['S']@G.T,U=P0*s['area']/(GAMMA-1)))
    contacts=[];nodes,weights=np.polynomial.legendre.leggauss(3)
    for i,j,a,b,normal,length in edges:
        if j<0:continue
        for z,w in zip(nodes,weights):
            contacts.append((i,j,.5*(a+b)+.5*z*(b-a),(b-a)/length,.5*length*w))
    return states,contacts

def vector(states):
    return np.concatenate([np.r_[s['v'],s['G'].ravel(),s['U']] for s in states])
def totals(states):
    return np.r_[sum(kinetic(s)+s['U'] for s in states),sum(angular(s) for s in states),
                 sum((s['m']*s['v'] for s in states),np.zeros(2)),sum(s['U'] for s in states)]
def run(initial,contacts,dt,reverse=False):
    states=list(initial);order=list(reversed(contacts)) if reverse else contacts
    for step in range(round(.02/dt)):
        # Palindromic composition of exact pair maps, second order for frozen geometry.
        for sequence in (order,list(reversed(order))):
            for i,j,p,t,weight in sequence:
                states[i],states[j],_=relax(states[i],states[j],p,t,weight,dt/2)
    return states

rows=[]
for n in [4,8]:
    for curved in [False,True]:
        states,contacts=setup(n,curved);initial=totals(states);rate=0.
        for i,j,p,t,weight in contacts:
            a,b=states[i],states[j];slip=(a['v']+a['G']@(p-a['c'])-b['v']-b['G']@(p-b['c']))@t
            rate+=weight*slip**2
        row=dict(n=n,field='curved' if curved else 'affine',contacts=len(contacts),initial_heat_rate=rate)
        if not curved:
            final=run(states,contacts,.01)
            row['null_max_state_change']=float(np.max(abs(vector(final)-vector(states))))
            assert row['null_max_state_change']<1e-12
        else:
            finals=[];order_errors=[];checks=[]
            for dt in [.01,.005,.0025]:
                final=run(states,contacts,dt);backward=run(states,contacts,dt,True)
                q=totals(final);err=q-initial
                assert np.max(abs(err[:4]))<1e-12 and err[4]>=0
                assert all(np.array_equal(a['R'],b['R']) for a,b in zip(states,final))
                finals.append(vector(final));order_errors.append(float(np.linalg.norm(vector(final)-vector(backward))))
                checks.append(dict(dt=dt,energy_error=float(err[0]),angular_error=float(err[1]),
                                   momentum_error=float(np.max(abs(err[2:4]))),generated_heat=float(err[4])))
            differences=[float(np.linalg.norm(finals[k]-finals[k+1])) for k in range(2)]
            row.update(checks=checks,successive_state_differences=differences,
                       convergence_ratio=differences[0]/differences[1],reverse_order_differences=order_errors)
            assert differences[1]<differences[0] and order_errors[-1]<order_errors[0]
        rows.append(row)
print(json.dumps(dict(scope=__doc__,cases=rows,
    limits=['Fixed geometry and internal faces only; no boundary traction, advection or gravity.',
    'Unit face damping is uncalibrated. Two spatial resolutions are not full convergence proof.',
    'Curved smooth fields can still be artificially heated. No Kepler stability claim.',
    'Heat is retained, residual covariance is not converted to heat.']),indent=2))
