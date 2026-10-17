"""Capture curved-contact failed bulk LP and separate basis/numerical causes."""
import json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import numpy as np
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import setup
from predictive_boundary_density import BoundaryDensity
import predictive_boundary_density as boundary_module
from bounded_moment_feasibility import minimum_moment_residual
from boundary_layer_feasibility_audit import boundary_parts
from face_moment_feasibility_audit import build
from vector_boundary_geometry import triangle_moments_many,phase_averages
from moment_aligned_density_audit import polygon_moments
import moment_ale_evolve as model


def inspect(entry):
    n=entry['n'];saved=json.loads(Path(entry['source']).read_text());raw=np.asarray(saved['raw']);geometry=np.asarray(saved['geometry'])
    _,y=setup('shear_contact',n,100.);a=y[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],0.)
    np.testing.assert_allclose(np.asarray([polygon_moments(v) for v in cells]),geometry,rtol=1e-12,atol=1e-14)
    captured=[];original_solve=boundary_module.solve
    def recording_solve(problem,*args,**kwargs):
        result=original_solve(problem,*args,**kwargs)
        if not result['feasible'] and not captured:captured.append((problem,result))
        return result
    boundary_module.solve=recording_solve
    provider=BoundaryDensity(mode='moment',segments=64,radial=1e-7,core_level=2,align_core=True,
        vector_geometry=True,construct_prior_witness=True,lp_retry=True,smooth_phase_floor=True,phase_fraction_floor=1e-7)
    error=None
    try:provider(cells,edges,raw)
    except RuntimeError as exc:error=repr(exc)
    finally:boundary_module.solve=original_solve
    if not captured:raise RuntimeError('original failure was not reproduced')
    if getattr(provider,'last_failed_stage',None)!='bulk':raise RuntimeError('expected original bulk failure, got '+str(error))
    problem,original_result=captured[0];cell=provider.last_failed_cell;prior=provider.last_attempted_priors[cell]
    means=raw[:,0]/geometry[:,0];low=float(means.min());high=float(means.max());rows=[]
    for level in (2,3,4):
        if level==2:p=problem
        else:
            parts=boundary_parts(cells[cell],1e-7,64,core_level=level,core_prior=prior,moment_function=triangle_moments_many)
            geo=triangle_moments_many(parts);values=None if prior['flat'] else phase_averages(parts,prior,low,high)
            p=build(cells[cell],raw[cell,:7],parts,prior,means[cell],low,high,geometry_moments=geo,area_prior=values)
        result=minimum_moment_residual(p['A'],p['b'],p['bounds'][:p['k']])
        np.savez('n'+str(n)+'_level'+str(level)+'_certificate.npz',A=p['A'],b=p['b'],bounds=np.asarray(p['bounds'][:p['k']]),
            **{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
        row=dict(n=n,cell=cell,core_level=level,subcells=p['k'],**{k:v for k,v in result.items() if not isinstance(v,np.ndarray)})
        rows.append(row);atomic('n'+str(n)+'_progress.json',dict(rows=rows))
        print('basis',n,cell,level,row['classification'],row.get('minimum_residual'),flush=True)
    report=dict(completed=True,n=n,failed_cell=cell,original_error=error,original_result={k:v for k,v in original_result.items() if k!='mass_fractions'},
        rows=rows,note='Only the first failed bulk cell is diagnosed. Reduced LP minimizes moment mismatch at unchanged density bounds; primal witnesses or dual lower bounds distinguish feasibility from original objective-solver failure. No face constraints, physical time evolution or moment changes.')
    atomic('n'+str(n)+'_diagnosis.json',report);return report


def run():
    cfg=json.loads(Path('basis_config.json').read_text())
    with ProcessPoolExecutor(max_workers=2,mp_context=mp.get_context('spawn')) as executor:results=list(executor.map(inspect,cfg['cases']))
    atomic('curved_basis_diagnosis.json',dict(completed=True,results=results))


if __name__=='__main__':run()
