"""Instantaneous mass-flux constrained weights, frozen numerical predictor.
Generator path unchanged. Independent dense-path derivative audit included.
"""
import json,time
from pathlib import Path
import numpy as np
from numpy.polynomial import Polynomial as P
from scipy.optimize import root
from scipy.integrate import solve_ivp
from predicted_spacetime_flux import FluxAudit
from dense_path_derivative import path_derivative

def run(kind,h,tol,n=2,analytic=False,max_step=1/16,root_xtol=1e-9):
    if not 0 < max_step <= 1:raise ValueError("max_step must be in (0,1]")
    started=time.monotonic();audit=FluxAudit(kind,n,h);p=audit.path;N=p.N
    w0=np.array([w(0.) for w in p.w]);calls=0;max_solve=0.
    def set_weights(w):p.w=[P([float(v)]) for v in w]
    def constrained(t,w):
        nonlocal calls,max_solve
        set_weights(w)
        def residual(z):
            nonlocal calls;calls+=1
            return audit.rate(t,np.r_[z,0.])[:-4].reshape(N,4)[:-1,0]
        if analytic:
            cache={}
            def evaluate(z):
                nonlocal calls
                key=z.tobytes()
                if cache.get('key')!=key:
                    calls+=1
                    flux,jac=audit.rate(t,np.r_[z,0.],return_jacobian=True)
                    cache.update(key=key,res=flux[:-4].reshape(N,4)[:-1,0],jac=jac[:-1,:-1])
                return cache
            solution=root(lambda z:evaluate(z)['res'],np.zeros(N-1),jac=lambda z:evaluate(z)['jac'],method='hybr',options={'xtol':root_xtol})
        else:
            solution=root(residual,np.zeros(N-1),method='hybr',options={'xtol':root_xtol})
        z=np.r_[solution.x,0.];flux=audit.rate(t,z);error=float(max(abs(flux[:-4].reshape(N,4)[:,0])))
        if not np.all(np.isfinite(z)) or error>2e-14:raise RuntimeError('instantaneous weight solve failed '+str(error))
        max_solve=max(max_solve,error);return z,flux
    def fun(t,y):
        z,f=constrained(t,y[:N]);return np.r_[z,f]
    y0=np.r_[w0,np.zeros(4*N+4)]
    sol=solve_ivp(fun,(0.,1.),y0,method='DOP853',rtol=tol,atol=tol*.01,max_step=max_step,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    samples=[];fd_scales=(2e-5,1e-5)
    for t in np.linspace(.025,.975,21):
        y=sol.sol(t);set_weights(y[:N]);rates=[]
        for eps in fd_scales:
            z=(sol.sol(t+eps)[:N]-sol.sol(t-eps)[:N])/(2*eps)
            f=audit.rate(t,z);rates.append(float(max(abs(f[:-4].reshape(N,4)[:,0])/p.mass)))
        # Fresh constraint solve only supplies an additional derivative comparison.
        exact_z,_=constrained(t,y[:N]);fd=(sol.sol(t+1e-5)[:N]-sol.sol(t-1e-5)[:N])/(2e-5)
        dz=path_derivative(sol.sol,t)[:N]
        denseflux,jac=audit.rate(t,dz,return_jacobian=True)
        dense_error=float(max(abs(denseflux[:-4].reshape(N,4)[:,0])/p.mass))
        samples.append(dict(dense_polynomial_relative_rate=dense_error,fd_to_polynomial_difference=float(max(abs(fd-dz))),jacobian_condition=float(np.linalg.cond(jac[:-1,:-1])),t=float(t),relative_mass_rate_normalized_time=rates,weight_derivative_error=float(max(abs(fd-exact_z)))))
    end=sol.y[:,-1];flux=end[N:-4].reshape(N,4);boundary=end[-4:]
    max_fd=max(s['relative_mass_rate_normalized_time'][-1] for s in samples)
    result=dict(kind=kind,n=n,analytic_jacobian=analytic,max_step=max_step,root_xtol=root_xtol,duration=h,rtol=tol,status='completed',seconds=time.monotonic()-started,calls=calls,accepted_steps=len(sol.t)-1,
        max_constraint_flux=max_solve,integrated_relative_mass_defect=float(max(abs(flux[:,0])/p.mass)),
        independent_path_relative_rate=max_fd,boundary_ledger_error=float(max(abs(flux.sum(0)-boundary))),
        max_weight_derivative_error=max(s['weight_derivative_error'] for s in samples))
    result['dense_polynomial_relative_rate']=max(s['dense_polynomial_relative_rate'] for s in samples)
    result['dense_polynomial_gate']=result['dense_polynomial_relative_rate']<1e-8
    result['sampled_path_gate']=bool(max_fd<1e-8 and result['integrated_relative_mass_defect']<1e-10 and result['boundary_ledger_error']<1e-10)
    tag=f'{kind}_n{n}_h{h:g}_tol{tol:g}_analytic{int(analytic)}'
    if max_step!=1/16:tag+=f'_maxstep{max_step:g}'
    if root_xtol!=1e-9:tag+=f'_root{root_xtol:g}'
    Path(tag+'.dense.json').write_text(json.dumps(dict(ts=sol.sol.ts.tolist(),segments=[dict(t_old=s.t_old,h=s.h,y_old=s.y_old.tolist(),F=s.F.tolist()) for s in sol.sol.interpolants])))
    Path(tag+'.path.json').write_text(json.dumps(dict(result=result,samples=samples,accepted_t=sol.t.tolist(),accepted_weights=sol.y[:N].T.tolist(),final_flux=flux.tolist(),boundary=boundary.tolist(),generator_polynomials=[[q.coef.tolist() for q in x] for x in p.x]),indent=2))
    return result
if __name__=='__main__':
    for kind in ('shear','wave'):
        for h in (.01,.05):
            for tol in (1e-8,1e-10):
                try:r=run(kind,h,tol)
                except Exception as e:r=dict(kind=kind,duration=h,rtol=tol,status='failed',error=repr(e))
                print(json.dumps(r),flush=True)
