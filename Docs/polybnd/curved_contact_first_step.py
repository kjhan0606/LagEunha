"""One guarded curved-contact step, reusing the approved initial RHS."""
import json,time,hashlib,os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import numpy as np
from shapely.geometry import Polygon
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import setup,diagnostic
from predictive_evolution_batch import assess,initialize_path_worker
from predictive_path_refinement import assess_adaptive
from predictive_boundary_density import BoundaryDensity
from factor_moment_evolve import rebuild,advance_factors_sylvester as advance
from shear_contact_check import affine_errors
from moment_aligned_density_audit import polygon_moments
import moment_ale_evolve as model
import kinetic_factor_rhs as kinetic


def load_saved_rhs(data):
    result=dict(data);result['dy']=np.asarray(result['dy']);diag=dict(result['diag'])
    for key in ('cov_linearization_factors','cov_A','cov_K','cov_Qfactor'):
        diag[key]=[np.asarray(v) for v in diag[key]]
    result['diag']=diag
    return result


def run():
    started=time.monotonic();cfg=json.loads(Path('curved_step_config.json').read_text())
    source=json.loads(Path('preserved/initial_state.json').read_text())
    rhs=load_saved_rhs(json.loads(Path('preserved/positive_rhs.json').read_text()))
    gate=json.loads(Path('preserved/source_rhs_report.json').read_text())
    if not gate['passed'] or gate['n']!=source['n']:raise RuntimeError('source RHS gate mismatch')
    mass=np.asarray(source['mass']);y=np.asarray(source['state']);L=[np.asarray(v) for v in source['factors']]
    n=source['n'];dt=float(cfg['dt']);horizon=1e-4
    setup('shear_contact',n,100.)
    options=dict(mode='moment',segments=64,radial=1e-7,core_level=3,align_core=True,
        vector_geometry=True,construct_prior_witness=True,lp_retry=True,smooth_phase_floor=True,phase_fraction_floor=1e-7)
    provider=BoundaryDensity(**options)
    d0=np.asarray(rhs['dy']);diag0=rhs['diag']
    halfL=advance(L,diag0,dt/2);half=rebuild(y+dt*d0/2,mass,halfL)
    dm,diag=kinetic.rhs(dt/2,half,mass,covariance_factors=halfL,face_density_provider=provider,
        face_segments=64,face_order=4,exact_bin_transport=True)
    newL=advance(L,diag,dt);new=rebuild(y+dt*dm,mass,newL)
    record=dict(t=0.,dt=dt,y=y,L=L,d0=d0,dm=dm,
        diag={k:diag[k] for k in ('cov_linearization_factors','cov_A','cov_K','cov_Qfactor')})
    atomic('pending.json',dict(t=0.,dt=dt,record=record,new=new,newL=newL,
        initial_rhs_reused=True,source_state_sha256=hashlib.sha256(y.tobytes()).hexdigest()))
    samples={};workers=int(os.environ.get('LAGFORCE_PATH_WORKERS','8'))
    with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize_path_worker,
            initargs=(n,100.,True,True,options)) as executor:
        for nodes in (17,65):
            checked=assess(record,mass,horizon,provider,nodes,samples,Path('pending_nodes.json'),
                'curved_first_step',executor,record_path=Path('pending.json'))
        checked=assess_adaptive(record,mass,horizon,provider,samples,Path('pending_refined_nodes.json'),
            'curved_first_step',executor,record_path=Path('pending.json'),max_nodes=1025)
    d=diagnostic(dt,new,mass,y);affine=affine_errors(new,mass)
    stage_mass=max(diag0['mass_flux'],diag['mass_flux'])
    passed=np.isfinite(checked['ratio']) and checked['ratio']<=1 and checked['pressure']>0 and checked['factor_R']>=-1e-20 and d['min_pressure']>0 and d['ledger_error']<1e-10 and max(affine.values())<1e-8 and d['pressure_departure']<1e-8 and stage_mass<1e-10
    result=dict(completed=True,passed=bool(passed),n=n,dt=dt,path_control_horizon=horizon,
        checked=checked,diagnostic=d,affine_errors=affine,stage_mass_flux=stage_mass,
        initial_rhs_reused=True,seconds=time.monotonic()-started,
        note='First curved-contact trial only. Reuses approved initial state/RHS without re-evaluation. Unsteered generators, fixed original mass array, original local-moment prior/epsilon1e-7/S64, refined core3. All original path budgets retained; no long-time/PDE convergence claim.')
    if passed:
        atomic('accepted_checkpoint.json',dict(t=dt,mass=mass,state=new,factors=newL,n=n,density_provider=options))
        reference=json.loads(Path('preserved/reference_moments.json').read_text());boundary=np.asarray(reference['reference_boundary'])
        boundary[:,0]+=dt*(.2+1.1*boundary[:,1]);high=Polygon(boundary)
        a=new[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],dt)
        _,closure=provider(cells,edges,np.column_stack([mass,a[:,3:]]));error=0.;exact=[]
        for cell,(parts,rhos) in zip(cells,provider.last_cells):
            inside=Polygon(cell).intersection(high)
            overlap=np.zeros(7) if inside.area<1e-28 else polygon_moments(np.asarray(inside.exterior.coords[:-1]))
            exact.append(polygon_moments(cell)+99*overlap)
            for triangle,rho in zip(parts,rhos):
                poly=Polygon(triangle);area=poly.intersection(high).area
                error+=abs(rho-100.)*area+abs(rho-1.)*max(0.,poly.area-area)
        exact=np.asarray(exact)
        result['final_score']=dict(density_L1_relative=error/mass.sum(),
            cell_mass_L1_relative=float(abs(exact[:,0]-mass).sum()/mass.sum()),
            exact_global_mass_relative_error=float(abs(exact[:,0].sum()-mass.sum())/mass.sum()),closure=closure)
        if result['final_score']['exact_global_mass_relative_error']>1e-11:raise RuntimeError('reference mass partition failed')
    atomic('curved_first_step_report.json',result);print('curved first step',passed,'ratio',checked['ratio'],flush=True)
    raise SystemExit(0 if passed else 2)


if __name__=='__main__':run()
