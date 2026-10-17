"""Same curved polygon contact across grids: static reconstruction only."""
import json,time,os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import multiprocessing as mp
import numpy as np
from shapely.geometry import Polygon
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import setup,D0
from moment_aligned_density_audit import polygon_moments
from density_moment_bounds import density_moment_bounds
from predictive_boundary_density import BoundaryDensity
from face_density_selection_audit import one_sided_interval
import moment_ale_evolve as model


def inspect(n,core_level=2,saved_source=None,segments=64):
    started=time.monotonic();_,y=setup('shear_contact',n,100.);a=y[:-4].reshape(-1,20)
    cells,edges=model.geometry(a[:,:2],a[:,2],0.)
    theta=np.arange(128)*2*np.pi/128
    boundary=.5+.22*np.column_stack([np.cos(theta),np.sin(theta)]);high=Polygon(boundary)
    if saved_source is None:
        raw=[];geometry=[]
        for vertices in cells:
            geo=polygon_moments(vertices);overlap=Polygon(vertices).intersection(high)
            moment=np.zeros(7) if overlap.area<1e-28 else polygon_moments(np.asarray(overlap.exterior.coords[:-1]))
            geometry.append(geo);raw.append(geo+99*moment)
        raw=np.asarray(raw);geometry=np.asarray(geometry)
    else:
        saved=json.loads(Path(saved_source).read_text());raw=np.asarray(saved['raw']);geometry=np.asarray(saved['geometry'])
        np.testing.assert_allclose(np.asarray([polygon_moments(v) for v in cells]),geometry,rtol=1e-12,atol=1e-14)
        np.testing.assert_array_equal(boundary,saved['reference_boundary'])
    tag='n'+str(n)+('_core'+str(core_level) if core_level!=2 else '')+('_s'+str(segments) if segments!=64 else '')
    target=polygon_moments(D0)+99*polygon_moments(boundary)
    reference_error=float(np.max(abs(raw.sum(0)-target)/np.maximum(1.,abs(target))))
    if reference_error>1e-11:raise RuntimeError('curved initial moment partition failed')
    means=raw[:,0]/geometry[:,0];low=float(means.min());upper=1.01*float(means.max());lower=.99*low
    certificates=[density_moment_bounds(g,r) for g,r in zip(geometry,raw)]
    violations=[]
    for i,c in enumerate(certificates):
        if c['density_eigen_max']>upper+1e-7*max(1.,upper) or c['density_eigen_min']<lower-1e-7*max(1.,upper):
            violations.append(dict(cell=i,**c,model_lower=lower,model_upper=upper))
    atomic(tag+'_initial_moments.json',dict(raw=raw,geometry=geometry,reference_boundary=boundary,
        note='True phase values1/100 initialize and score only; provider receives numerical moments and cell means.'))
    base=dict(n=n,core_level=core_level,segments=segments,saved_source=saved_source,total_mass=float(raw[:,0].sum()),reference_moment_error=reference_error,
        power_cell_backend='active_all_constraints' if os.environ.get('LAGFORCE_ACTIVE_POWER_CELLS','0')=='1' else 'exhaustive',
        inferred_phase_levels=[low,float(means.max())],model_density_bounds=[lower,upper],
        required_upper_from_moments=max(c['density_eigen_max'] for c in certificates),
        bound_violations=violations)
    if violations:
        result=dict(base,passed=False,status='density_bounds_incompatible_with_physical_moments',seconds=time.monotonic()-started)
        atomic(tag+'_report.json',result);print('curved',n,result['status'],flush=True);return result
    provider=BoundaryDensity(mode='moment',segments=segments,radial=1e-7,core_level=core_level,align_core=True,
        vector_geometry=True,construct_prior_witness=True,lp_retry=True,smooth_phase_floor=True,phase_fraction_floor=1e-7)
    try:
        _,diag=provider(cells,edges,raw);error=face_error=face_norm=0.
        for parts,rhos in provider.last_cells:
            for triangle,rho in zip(parts,rhos):
                poly=Polygon(triangle);inside=poly.intersection(high).area
                error+=abs(rho-100.)*inside+abs(rho-1.)*max(0.,poly.area-inside)
        for edge_id,(i,j,start,end,normal,length) in enumerate(edges):
            for cell,outward in ((i,normal),(j,-normal)):
                if cell<0:continue
                interval=one_sided_interval(boundary,start,end,-outward)
                values=provider.last_profiles[(cell,edge_id)];S=len(values)
                for k,rho in enumerate(values):
                    inside=0. if interval is None else max(0.,min((k+1)/S,interval[1])-max(k/S,interval[0]))
                    outside=max(0.,1/S-inside)
                    face_error+=length*(abs(rho-100.)*inside+abs(rho-1.)*outside)
                    face_norm+=length*(100.*inside+outside)
        result=dict(base,passed=bool(diag['closure_moment_error']<1e-8 and diag['closure_min_density']>0 and diag['closure_volume_cap_excess']<1e-8),
            status='reconstructed',density_L1_relative=error/raw[:,0].sum(),face_pointwise_L1_relative=face_error/face_norm,
            closure=diag,seconds=time.monotonic()-started)
    except Exception as exc:result=dict(base,passed=False,status='closure_unavailable',error=repr(exc),seconds=time.monotonic()-started)
    atomic(tag+'_report.json',result);print('curved',n,core_level,result['status'],result.get('density_L1_relative'),flush=True);return result


def run():
    results=[]
    path=Path('curve_sweep_config.json')
    entries=json.loads(path.read_text())['cases'] if path.exists() else [dict(n=n,core_level=2) for n in (4,8,16)]
    with ProcessPoolExecutor(max_workers=2,mp_context=mp.get_context('spawn')) as executor:
        futures=[executor.submit(inspect,e['n'],e.get('core_level',2),e.get('source'),e.get('segments',64)) for e in entries]
        for f in as_completed(futures):results.append(f.result());atomic('curved_progress.json',dict(results=results))
    report=dict(completed=True,all_reconstructions_passed=all(r['passed'] for r in results),results=results,
        note='Same physical128-sided curved contact and phase ratio100 across grids at t0. Static density reconstruction only, no Euler evolution or spatial PDE convergence claim. Failed necessary bounds have squared-polynomial moment certificates independent of the finite LP basis.')
    atomic('curved_contact_report.json',report)


if __name__=='__main__':run()
