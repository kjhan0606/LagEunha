"""Evolving moment/factor midpoint steps with independent dense-path defects.

All evolution stages use kinetic mass solves and unsteered material sites.
Audit-only geometry rates are the derivatives of the actual quadratic step
path, never a replacement constraint solve. This measures whether stagewise
mass constraints survive interpolation when the fluid state also evolves.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from adaptive_moment_batch import setup, diagnostic
from kinetic_factor_rhs import rhs
from factor_moment_evolve import rebuild,advance_factors_sylvester as advance


def dense_state(y,L,d0,dm,diag_mid,mass,dt,s):
    path=y+dt*(s*d0+s*s*(dm-d0))
    factors=advance(L,diag_mid,dt*s) if s else [f.copy() for f in L]
    path=rebuild(path,mass,factors)
    derivative=d0+2*s*(dm-d0)
    return path,factors,derivative


def run(steps,end=.05,n=3,amplitude=.01):
    folder=Path(f'steps_{steps}');folder.mkdir(exist_ok=False)
    start=time.monotonic();dt=end/steps
    mass,y=setup('wave',n,amplitude);L=[np.zeros((2,2)) for _ in mass]
    y=rebuild(y,mass,L);y0=y.copy();N=len(mass)
    history=[];max_rate=0.;max_generator=0.;max_stage=0.
    min_pressure=float('inf');min_gram=0.;cumulative=np.zeros(N);coarse_cumulative=np.zeros(N);max_prefix=0.
    last_t=0.
    try:
        for k in range(steps):
            t=k*dt
            d0,diag0=rhs(t,y,mass,covariance_factors=L)
            halfL=advance(L,diag0,dt/2)
            half=rebuild(y+dt*d0/2,mass,halfL)
            dm,diag_mid=rhs(t+dt/2,half,mass,covariance_factors=halfL)
            newL=advance(L,diag_mid,dt)
            new=rebuild(y+dt*dm,mass,newL)
            max_stage=max(max_stage,diag0['mass_flux'],diag_mid['mass_flux'])
            # Rebuild Pi along dense output using Gram factors; x,w,q,Q,p,D,E
            # use the quadratic midpoint extension. No eigenvalue clipping.
            step_rows=[]
            step_flux=[]
            for j,s in enumerate(np.linspace(0.,1.,33)):
                state,factors,derivative=dense_state(y,L,d0,dm,diag_mid,mass,dt,s)
                rates=derivative[:-4].reshape(N,20)
                _,d=rhs(t+s*dt,state,mass,covariance_factors=factors,
                        audit_geometry_rates=(rates[:,:2],rates[:,2]))
                physical=float(np.max(abs(d['mass_rates'])/mass))
                normalized=end*physical
                max_rate=max(max_rate,normalized)
                max_generator=max(max_generator,d['generator_velocity_defect'])
                min_pressure=min(min_pressure,d['min_pressure'])
                min_gram=min(min_gram,min(float(np.linalg.eigvalsh(f@f.T).min()/m) for f,m in zip(factors,mass)))
                step_rows.append(dict(node=j,s=float((t+s*dt)/end),mass_rate=normalized,
                    generator_velocity_defect=d['generator_velocity_defect']))
                step_flux.append(d['mass_rates'])
            # Nested composite Simpson estimates; signed cell defects matter.
            flux=np.asarray(step_flux)
            def simpson(values):
                return dt/(3*(len(values)-1))*(values[0]+values[-1]+4*values[1:-1:2].sum(0)+2*values[2:-1:2].sum(0))
            fine=simpson(flux);coarse=simpson(flux[::2]);cumulative+=fine;coarse_cumulative+=coarse
            max_prefix=max(max_prefix,float(np.max(abs(cumulative)/mass)))
            hn=diagnostic(t+dt,new,mass,y0)
            np.testing.assert_allclose(dense_state(y,L,d0,dm,diag_mid,mass,dt,1.)[0],new,rtol=0,atol=1e-14)
            row=dict(step=k+1,t=t+dt,samples=step_rows,
                max_rate=max(q['mass_rate'] for q in step_rows),
                quadrature_mass_difference=float(np.max(abs(fine-coarse)/mass)),
                ledger_error=hn['ledger_error'],min_pressure=hn['min_pressure'])
            history.append(row)
            with (folder/'history.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
            # Accepted dynamical state only; no acceptance claim for dense path.
            y,L=new,newL;last_t=t+dt
            checkpoint=dict(time=last_t,mass=mass.tolist(),state=y.tolist(),
                factors=[f.tolist() for f in L],step=k+1,steps=steps,
                dense_mass_defect=cumulative.tolist())
            temp=folder/'checkpoint.tmp';temp.write_text(json.dumps(checkpoint));temp.replace(folder/'checkpoint.json')
            print(json.dumps(dict(steps=steps,step=k+1,time=last_t,max_mass_path_residual=max_rate)),flush=True)
        final=diagnostic(end,y,mass,y0)
        result=dict(status='completed',steps=steps,dt=dt,duration=end,n=n,
            max_path_mass_rate_normalized=max_rate,max_generator_velocity_defect=max_generator,
            max_stage_mass_flux=max_stage,min_pressure=min_pressure,min_factor_R=min_gram,
            integrated_relative_mass_defect=float(np.max(abs(cumulative)/mass)),
            max_prefix_mass_defect=max_prefix,
            quadrature_mass_difference=float(np.max(abs(cumulative-coarse_cumulative)/mass)),
            max_step_quadrature_mass_difference=max(h['quadrature_mass_difference'] for h in history),
            max_ledger_error=max(h['ledger_error'] for h in history),final=final,
            seconds=time.monotonic()-start,samples=steps*33,
            note='evolving fluid; stagewise unsteered/fixed-mass scheme; dense path consistency tested separately')
        result['state_gate']=bool(min_pressure>0 and min_gram>=-1e-20 and result['max_ledger_error']<1e-10)
        result['mass_path_gate']=bool(max_rate<1e-8 and result['integrated_relative_mass_defect']<1e-10
                                     and result['quadrature_mass_difference']<1e-11)
    except Exception as exc:
        result=dict(status='failed',steps=steps,time=last_t,error=repr(exc),seconds=time.monotonic()-start)
    (folder/'result.json').write_text(json.dumps(result,indent=2))
    with Path('results.jsonl').open('a') as f:f.write(json.dumps(result)+'\n')
    print(json.dumps(result),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('steps',type=int,nargs='+');args=parser.parse_args()
    results=[run(count) for count in args.steps]
    summary=dict(results=results)
    if len(results)==3 and all(r['status']=='completed' for r in results):
        states=[np.asarray(json.loads(Path(f'steps_{c}/checkpoint.json').read_text())['state'])[:-4].reshape(-1,20) for c in args.steps]
        groups={'sites':slice(0,2),'weights':slice(2,3),'q':slice(3,5),'Q':slice(5,9),
                'p':slice(9,11),'D':slice(11,15),'Pi':slice(15,19),'E':slice(19,20)}
        convergence={}
        for name,sl in groups.items():
            diffs=[float(np.linalg.norm(x[:,sl]-y[:,sl])) for x,y in zip(states[:-1],states[1:])]
            convergence[name]=dict(successive_L2_differences=diffs,
                observed_order=float(np.log2(diffs[0]/diffs[1])) if min(diffs)>1e-13 else None)
        summary['self_convergence']=convergence
    Path('summary.json').write_text(json.dumps(summary,indent=2))
