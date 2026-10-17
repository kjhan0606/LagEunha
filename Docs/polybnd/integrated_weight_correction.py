"""Endpoint integrated-mass constraint for a frozen numerical predictor.
Sites unchanged; weights gain tau^2 delta. Not an all-time mass constraint.
"""
import json,time
import numpy as np
from numpy.polynomial import Polynomial as P
from predicted_spacetime_flux import FluxAudit
from moving_power_spacetime import events

def run(kind,h):
    started=time.monotonic();audit=FluxAudit(kind,2,h);p=audit.path;original=[w.copy() for w in p.w];calls=0
    original_sites=[[q.coef.copy() for q in x] for x in p.x]
    def evaluate(delta,order=32):
        nonlocal calls;calls+=1
        p.w=[w+P([0,0,float(v)]) for w,v in zip(original,delta)]
        x,w=p.state(0);ts,_=events(x,np.zeros_like(x),w,np.zeros_like(w),constraint_builder=p.constraints)
        value=audit.integral(ts,order);return value[:-4].reshape(p.N,4)[:,0],value,ts
    delta=np.zeros(p.N);r,v,ts=evaluate(delta);before=float(max(abs(r)/p.mass));history=[]
    for iteration in range(12):
        error=float(max(abs(r)/p.mass));history.append(error)
        if error<2e-12:break
        J=np.zeros((p.N-1,p.N-1));eps=1e-6
        for j in range(p.N-1):
            d=np.zeros(p.N);d[j]=eps
            J[:,j]=(evaluate(delta+d)[0][:-1]-evaluate(delta-d)[0][:-1])/(2*eps)
        step=np.r_[np.linalg.solve(J,-r[:-1]),0.]
        for k in range(18):
            trial=delta+step*2.**(-k)
            try:rr,vv,tt=evaluate(trial)
            except (RuntimeError,ValueError):continue
            if max(abs(rr)/p.mass)<error:delta=trial;r,v,ts=rr,vv,tt;break
        else:raise RuntimeError('integrated mass Newton line search failed')
    else:raise RuntimeError('integrated mass Newton not converged')
    r32,v32,ts=evaluate(delta,32);r64,v64,_=evaluate(delta,64)
    relative=float(max(abs(r64)/p.mass));refinement=float(max(abs(v64-v32)))
    rates=[audit.rate(t)[:-4].reshape(p.N,4)[:,0]/h for t in (.25,.5,.75)]
    for x,old in zip(p.x,original_sites):
        for q,c in zip(x,old):np.testing.assert_array_equal(q.coef,c)
    flux=v64[:-4].reshape(p.N,4);ledger=float(max(abs(flux.sum(0)-v64[-4:])))
    return dict(kind=kind,n=2,duration=h,status='completed',seconds=time.monotonic()-started,evaluations=calls,
        before_relative_mass_defect=before,after_relative_mass_defect=relative,weight_correction=delta.tolist(),newton_history=history,
        quadrature_difference=refinement,boundary_ledger_error=ledger,
        max_sampled_instantaneous_relative_mass_rate=float(max(np.max(abs(r)/p.mass) for r in rates)),
        endpoint_gate=bool(relative<1e-10 and refinement<1e-9 and ledger<1e-11),generator_path_unchanged=True,
        limitation='Frozen predictor, endpoint integrated constraint only; intermediate fixed mass and full hydro evolution remain unsolved.')
if __name__=='__main__':
    for kind in ('shear','wave'):
        for h in (.01,.05):
            try:r=run(kind,h)
            except Exception as e:r=dict(kind=kind,duration=h,status='failed',error=repr(e))
            print(json.dumps(r),flush=True)
