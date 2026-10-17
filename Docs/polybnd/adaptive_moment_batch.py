"""Diagnostic batch of numerical-state moment evolution with adaptive DOP853.
No exact interior solution is injected after initialization. This deliberately
retains the old R=0 affine moment flux closure to expose its failure modes.
It does NOT claim to implement the implicit capacity path for general density.
"""
import json,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import moment_ale_evolve as model
from periodic2d_hllc import quadrature

D0=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])

def install_affine_shear_domain():
    """Install only the boundary law; do not construct a different interior IC."""
    def velocity(p):return np.array([.2+1.1*p[1],0.])
    def domain(t):return D0+np.array([t*velocity(p) for p in D0])
    model.domain=domain;model.boundary_velocity=velocity

def setup(kind,n,amplitude):
    def bvelocity(p):
        if kind=='moving_contact':return np.array([amplitude,-.5*amplitude])
        return np.array([.2+1.1*p[1],0.]) if kind in ('shear','shear_contact') else np.zeros(2)
    def domain(t):return D0+np.array([t*bvelocity(p) for p in D0])
    model.domain=domain;model.boundary_velocity=bvelocity
    xx,yy=np.meshgrid((np.arange(n)+.5)/n,(np.arange(n)+.5)/n)
    x=np.column_stack([xx.ravel(),yy.ravel()])+np.random.default_rng(731).uniform(-.15,.15,(n*n,2))/n
    polygons,_=model.geometry(x,np.zeros(n*n),0);rows=[]
    for site,poly in zip(x,polygons):
        r,w=quadrature(poly)
        v=np.array([bvelocity(p) for p in r]);P=np.full(len(w),.6);rho=np.ones(len(w))
        if kind in ('contact','moving_contact'):rho[:]=1. if site[0]<.5 else 100.
        if kind=='shear_contact':rho[:]=1. if site[0]<.5 else amplitude
        if kind=='wave':
            delta=amplitude*np.cos(2*np.pi*r[:,0])*np.cos(2*np.pi*r[:,1]);P+=delta;rho+=delta
        dm=w*rho;Pi=v.T@(dm[:,None]*v)
        rows.append(np.r_[sum(dm),dm@r,(r.T@(dm[:,None]*r)).ravel(),dm@v,
                         (r.T@(dm[:,None]*v)).ravel(),Pi.ravel(),sum(w*P)/(model.GAMMA-1)+.5*np.trace(Pi)])
    raw=np.array(rows);mass=raw[:,0].copy()
    return mass,np.r_[np.column_stack([x,np.zeros(n*n),raw[:,1:]]).ravel(),np.zeros(4)]


def diagnostic(t,y,mass,y0):
    a=y[:-4].reshape(len(mass),20);a0=y0[:-4].reshape(len(mass),20)
    raw=np.column_stack([mass,a[:,3:]]);r0=np.column_stack([mass,a0[:,3:]])
    polys,edges=model.geometry(a[:,:2],a[:,2],t);states=[model.decode(p,r) for p,r in zip(polys,raw)]
    ledger=np.r_[raw[:,7:9].sum(0)-r0[:,7:9].sum(0),sum(raw[:,17]-r0[:,17]),sum(raw[:,10]-raw[:,11]-r0[:,10]+r0[:,11])]-y[-4:]
    return dict(t=float(t),min_R=min(float(np.linalg.eigvalsh(s['R']).min()/m) for s,m in zip(states,mass)),
                max_R=max(float(np.max(abs(s['R']))/m) for s,m in zip(states,mass)),
                min_pressure=min(s['P'] for s in states),max_pressure=max(s['P'] for s in states),
                pressure_departure=max(abs(s['P']-.6) for s in states),
                ledger_error=float(max(abs(ledger))),topology=sorted((i,j) for i,j,*_ in edges if j>=0))


def run(kind,n,amplitude,rtol,end=.02):
    mass,y0=setup(kind,n,amplitude);calls=0;min_stage_R=0.;max_mass_flux=0.;start=time.monotonic()
    last=[0.,y0.copy()]
    def fun(t,y):
        nonlocal calls,min_stage_R,max_mass_flux
        dy,d=model.rhs(t,y,mass);calls+=1
        min_stage_R=min(min_stage_R,d['residual_covariance_min']);max_mass_flux=max(max_mass_flux,d['mass_flux'])
        last[:]=[float(t),y.copy()];return dy
    tag=f'{kind}_n{n}_a{amplitude:g}_tol{rtol:g}'
    try:
        sol=solve_ivp(fun,(0.,end),y0,method='DOP853',rtol=rtol,atol=rtol*1e-2,max_step=.005,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        history=[diagnostic(t,sol.sol(t),mass,y0) for t in np.linspace(0,end,21)]
        final=diagnostic(end,sol.y[:,-1],mass,y0)
        Path(tag+'.history.json').write_text(json.dumps(history,indent=2))
        Path(tag+'.checkpoint.json').write_text(json.dumps(dict(time=end,mass=mass.tolist(),state=sol.y[:,-1].tolist()),indent=2))
        return dict(kind=kind,n=n,amplitude=amplitude,rtol=rtol,end=end,status='completed',calls=calls,
                    accepted_steps=len(sol.t)-1,seconds=time.monotonic()-start,min_stage_R=min_stage_R,
                    min_sampled_R=min(h['min_R'] for h in history),max_mass_flux=max_mass_flux,final=final)
    except Exception as e:
        Path(tag+'.failure_state.json').write_text(json.dumps(dict(time=last[0],mass=mass.tolist(),state=last[1].tolist())))
        return dict(kind=kind,n=n,amplitude=amplitude,rtol=rtol,status='failed',error=repr(e),calls=calls,seconds=time.monotonic()-start,min_stage_R=min_stage_R)

if __name__=='__main__':
    import sys
    if '--compatible' in sys.argv:
        from compatible_moment_rhs import rhs
        model.rhs=rhs
    cases=[('shear',n,0.,tol) for n in (3,4) for tol in (1e-6,1e-8,1e-10)]
    cases += [('contact',4,0.,1e-10)]
    cases += [('wave',4,a,tol) for a in (.001,.1) for tol in (1e-6,1e-8,1e-10)]
    for args in cases:print(json.dumps(run(*args)),flush=True)
