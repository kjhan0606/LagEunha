"""Full instantaneous coupling gates for the saved curved-contact moments."""
import json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import numpy as np
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import setup,diagnostic
from moment_aligned_density_audit import polygon_moments
from predictive_boundary_density import BoundaryDensity
from predictive_coupling_batch import metrics
from continue_common_contact import valid_metric
from shear_contact_check import affine_errors
import moment_ale_evolve as model
import kinetic_factor_rhs as kinetic


def affine_state_from_position_moments(skeleton,raw,geometry,pressure=.6):
    """Exact raw moments for u=(.2+1.1*y,0), with zero kinetic covariance."""
    raw=np.asarray(raw);mass=raw[:,0].copy();state=np.array(skeleton,copy=True)
    a=state[:-4].reshape(-1,20);q=raw[:,1:3];Q=raw[:,3:7].reshape(-1,2,2)
    u=np.array([.2,0.]);B=np.array([[0.,1.1],[0.,0.]])
    momentum=mass[:,None]*u+q@B.T
    mixed=np.einsum('ni,j->nij',q,u)+Q@B.T
    bq=q@B.T
    kinetic_raw=mass[:,None,None]*np.outer(u,u)+np.einsum('ni,j->nij',bq,u)+np.einsum('i,nj->nij',u,bq)+B@Q@B.T
    a[:,3:9]=raw[:,1:7];a[:,9:11]=momentum
    a[:,11:15]=mixed.reshape(-1,4);a[:,15:19]=kinetic_raw.reshape(-1,4)
    a[:,19]=pressure*np.asarray(geometry)[:,0]/(model.GAMMA-1)+.5*np.trace(kinetic_raw,axis1=1,axis2=2)
    return mass,state,[np.zeros((2,2)) for _ in mass]


def inspect(entry):
    started=time.monotonic();n=entry['n'];saved=json.loads(Path(entry['source']).read_text())
    raw=np.asarray(saved['raw']);geometry=np.asarray(saved['geometry'])
    _,skeleton=setup('shear_contact',n,100.);a=skeleton[:-4].reshape(-1,20)
    polys,edges=model.geometry(a[:,:2],a[:,2],0.)
    np.testing.assert_allclose(np.asarray([polygon_moments(v) for v in polys]),geometry,rtol=1e-12,atol=1e-14)
    mass,state,factors=affine_state_from_position_moments(skeleton,raw,geometry)
    initial=affine_errors(state,mass)
    atomic('n'+str(n)+'_initial_state.json',dict(mass=mass,state=state,factors=factors,t=0.,
        geometry=geometry,n=n,source=entry['source'],
        note='Physical m/q/Q copied unchanged; remaining raw moments initialize the exact affine velocity and uniform pressure.'))
    provider=BoundaryDensity(mode='moment',segments=64,radial=1e-7,core_level=3,align_core=True,
        vector_geometry=True,construct_prior_witness=True,lp_retry=True,smooth_phase_floor=True,phase_fraction_floor=1e-7)
    try:
        args=dict(face_density_provider=provider,face_segments=64,face_order=4)
        full,df=model.rhs(0.,state,mass,**args)
        positive,dk=kinetic.rhs(0.,state,mass,covariance_factors=factors,exact_bin_transport=True,**args)
        low,_=kinetic.rhs(0.,state,mass,covariance_factors=factors,exact_bin_transport=True,**dict(args,face_order=2))
        difference=float(np.max(abs(positive-low)/np.maximum(1.,abs(positive))))
        vals=dict(full_hllc=metrics(full,df,mass),positive_kinetic=metrics(positive,dk,mass))
        passed=all(valid_metric(v) for v in vals.values()) and difference<1e-8 and max(initial.values())<1e-10
        result=dict(n=n,passed=bool(passed),initial_affine_errors=initial,metrics=vals,
            split_order2_vs4_difference=difference,seconds=time.monotonic()-started)
        atomic('n'+str(n)+'_positive_rhs.json',dict(dy=positive,diag=dk))
    except Exception as exc:result=dict(n=n,passed=False,error=repr(exc),seconds=time.monotonic()-started)
    atomic('n'+str(n)+'_rhs_report.json',result);print('curved RHS',n,result['passed'],flush=True);return result


def run():
    cfg=json.loads(Path('curved_rhs_config.json').read_text())
    with ProcessPoolExecutor(max_workers=2,mp_context=mp.get_context('spawn')) as executor:
        results=list(executor.map(inspect,cfg['cases']))
    report=dict(completed=True,passed=all(r['passed'] for r in results),results=results,
        note='Instantaneous full-HLLC and positive-kinetic RHS only. Weights are solved from numerical moments; no exact interface is supplied. No time evolution or PDE convergence claim.')
    atomic('curved_rhs_report.json',report)
    raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()
