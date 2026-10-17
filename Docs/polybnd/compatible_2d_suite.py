"""Numerical-state 2D acoustic/Gresho convergence with compatible moment RHS.
Stationary reflecting outer box; Gresho support lies inside. No reference flux.
"""
import json,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import moment_ale_evolve as model
from compatible_moment_rhs import rhs
from adaptive_moment_batch import setup,diagnostic
from periodic2d_hllc import quadrature


def gresho(points):
    d=points-.5;r=np.linalg.norm(d,axis=1);vphi=np.where(r<.2,5*r,np.where(r<.4,2-5*r,0.))
    vel=np.column_stack([-d[:,1],d[:,0]])*(vphi/np.maximum(r,1e-100))[:,None]
    P=np.full(len(r),3+4*np.log(2.));inner=r<.2;mid=(r>=.2)&(r<.4)
    P[inner]=5+12.5*r[inner]**2
    P[mid]=9+12.5*r[mid]**2-20*r[mid]+4*np.log(r[mid]/.2)
    return vel,P


def initial(kind,n):
    mass,y=setup('wave',n,.001 if kind=='wave' else 0.)
    if kind=='wave':return mass,y
    a=y[:-4].reshape(n*n,20);polys,_=model.geometry(a[:,:2],a[:,2],0.)
    for i,poly in enumerate(polys):
        r,w=quadrature(poly);v,P=gresho(r);Pi=v.T@(w[:,None]*v)
        raw=np.r_[sum(w),w@r,(r.T@(w[:,None]*r)).ravel(),w@v,
                  (r.T@(w[:,None]*v)).ravel(),Pi.ravel(),sum(w*P)/(model.GAMMA-1)+.5*np.trace(Pi)]
        mass[i]=raw[0];a[i,3:]=raw[1:]
    return mass,y


def compare(kind,t,y,mass):
    a=y[:-4].reshape(len(mass),20);polys,_=model.geometry(a[:,:2],a[:,2],t)
    errors=np.zeros(3);kinetic=0.;volume=0.
    for i,poly in enumerate(polys):
        raw=np.r_[mass[i],a[i,3:]];s=model.decode(poly,raw);r,w=quadrature(poly);area=sum(w);volume+=area
        if kind=='wave':
            omega=2*np.pi*np.sqrt(2);eps=.001
            c=np.cos(2*np.pi*r);sn=np.sin(2*np.pi*r)
            pressure=.6+eps*c[:,0]*c[:,1]*np.cos(omega*t);rho=1+pressure-.6
            v=eps/np.sqrt(2)*np.column_stack([sn[:,0]*c[:,1],c[:,0]*sn[:,1]])*np.sin(omega*t)
            scale=eps
        else:v,pressure=gresho(r);rho=np.ones(len(w));scale=1.
        refmass=w@rho;refv=(w*rho)@v/refmass
        errors[0]+=area*abs(s['P']-w@pressure/area)/scale
        errors[1]+=area*np.linalg.norm(s['v']-refv)/scale
        errors[2]+=area*abs(mass[i]/area-refmass/area)/scale
        kinetic+=.5*np.trace(raw[13:17].reshape(2,2))
    return dict(pressure_L1=float(errors[0]/volume),velocity_L1=float(errors[1]/volume),density_L1=float(errors[2]/volume),kinetic=kinetic)


def run(kind,n,end,max_step=.0025,rtol=1e-9):
    mass,y0=initial(kind,n);start=time.monotonic();count=0;stage_R=0.;last=[0.,y0.copy()]
    def f(t,y):
        nonlocal count,stage_R
        dy,d=rhs(t,y,mass);count+=1;stage_R=min(stage_R,d['residual_covariance_min']);last[:]=[t,y.copy()];return dy
    tag=f'{kind}_n{n}_dt{max_step:g}'
    try:
        sol=solve_ivp(f,(0,end),y0,method='DOP853',rtol=rtol,atol=rtol*.01,max_step=max_step,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        history=[diagnostic(t,sol.sol(t),mass,y0) for t in np.linspace(0,end,11)]
        final=history[-1];err=compare(kind,end,sol.y[:,-1],mass);err0=compare(kind,0.,y0,mass)
        result=dict(kind=kind,n=n,end=end,max_step=max_step,rtol=rtol,status='completed',rhs_calls=count,seconds=time.monotonic()-start,
                    min_stage_R=stage_R,min_sampled_R=min(d['min_R'] for d in history),final=final,errors=err,
                    kinetic_fraction=err['kinetic']/err0['kinetic'] if err0['kinetic'] else None)
        Path(tag+'.history.json').write_text(json.dumps(history,indent=2))
        Path(tag+'.checkpoint.json').write_text(json.dumps(dict(t=end,mass=mass.tolist(),state=sol.y[:,-1].tolist()),indent=2))
    except Exception as e:
        result=dict(kind=kind,n=n,end=end,status='failed',error=repr(e),rhs_calls=count,seconds=time.monotonic()-start,min_stage_R=stage_R,last_rhs_time=float(last[0]))
        Path(tag+'.failure.json').write_text(json.dumps(dict(t=float(last[0]),mass=mass.tolist(),state=last[1].tolist()),indent=2))
    return result

if __name__=='__main__':
    import sys
    if '--residual-stress' in sys.argv or '--linear-pressure' in sys.argv:
        from functools import partial
        rhs=partial(rhs,residual_stress='--residual-stress' in sys.argv,linear_pressure='--linear-pressure' in sys.argv)
    for kind,end in [('wave',.05),('gresho',.02)]:
        for n in (4,6,8):print(json.dumps(run(kind,n,end)),flush=True)
    print(json.dumps(run('wave',4,.05,max_step=.00125,rtol=1e-10)),flush=True)
