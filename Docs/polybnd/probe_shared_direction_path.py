"""Local vs numerical global direction from one preserved N8 state."""
import json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import numpy as np
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import diagnostic
from predictive_evolution_batch import initialize_case,density_options,initialize_path_worker,assess
from predictive_boundary_density import BoundaryDensity
from predictive_coupling_batch import metrics
from factor_moment_evolve import rebuild,advance_factors_sylvester as advance
from shear_contact_check import affine_errors
import kinetic_factor_rhs as kinetic


def run():
    source=json.loads(Path('captured_checkpoint.json').read_text());initial=json.loads(Path('initial_config_checkpoint.json').read_text())
    cfg=initial['config'];mass,_=initialize_case(cfg);np.testing.assert_array_equal(mass,source['mass'])
    y=np.asarray(source['state']);L=[np.asarray(f) for f in source['factors']];t=source['t'];dt=1.25e-6;rows=[]
    for mode in ('moment','global_moment'):
        started=time.monotonic();options=density_options(cfg,True);options['mode']=mode
        options['refine_moment_direction']=mode=='global_moment';provider=BoundaryDensity(**options)
        args=dict(face_density_provider=provider,face_segments=64,face_order=4,exact_bin_transport=True)
        d0,dg0=kinetic.rhs(t,y,mass,covariance_factors=L,**args)
        halfL=advance(L,dg0,dt/2);half=rebuild(y+dt*d0/2,mass,halfL)
        dm,dgm=kinetic.rhs(t+dt/2,half,mass,covariance_factors=halfL,**args)
        newL=advance(L,dgm,dt);new=rebuild(y+dt*dm,mass,newL)
        keys=('cov_linearization_factors','cov_A','cov_K','cov_Qfactor')
        record=dict(t=t,dt=dt,y=y,L=L,d0=d0,dm=dm,diag={k:dgm[k] for k in keys})
        atomic(mode+'_trial.json',dict(record=record,new=new,newL=newL,options=options))
        samples={}
        with ProcessPoolExecutor(max_workers=16,mp_context=mp.get_context('spawn'),initializer=initialize_path_worker,
                initargs=(cfg['n'],cfg['amplitude'],True,True,options)) as executor:
            control=assess(record,mass,1e-4,provider,17,samples,mode+'_nodes.json',mode,executor)
            audit=assess(record,mass,1e-4,provider,65,samples,mode+'_nodes.json',mode,executor)
        corrected=new.copy();corrected[-4:]-=y[-4:];d=diagnostic(t+dt,corrected,mass,y)
        affine=affine_errors(new,mass);m0=metrics(d0,dg0,mass);mm=metrics(dm,dgm,mass)
        passed=audit['ratio']<=1 and d['ledger_error']<1e-10 and d['pressure_departure']<1e-8 and max(affine.values())<1e-8
        passed=passed and max(m0['mass_flux'],mm['mass_flux'])<1e-10 and audit['generator']<1e-8 and audit['pressure']>0 and audit['factor_R']>=-1e-20
        rows.append(dict(mode=mode,passed=bool(passed),source_time=t,dt=dt,options=options,control=control,audit=audit,
            base=m0,midpoint=mm,diagnostic=d,affine=affine,seconds=time.monotonic()-started))
        atomic('shared_direction_progress.json',dict(results=rows));print('shared direction',mode,passed,audit['ratio'],flush=True)
    report=dict(completed=True,passed=rows[-1]['passed'],results=rows,
        note='Explicit single global direction hypothesis for the planar benchmark, inferred from numerical global m/q/Q on the numerical convex domain. Per-cell offset still inferred from numerical mean; local m/q/Q, positivity and L1-cap enforced. No reference plane in RHS. Independent trial; ongoing local-moment baseline unchanged.')
    atomic('shared_direction_report.json',report);raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()
