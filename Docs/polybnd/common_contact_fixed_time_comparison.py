"""Same-time common-IC grid table and separate frozen-state face-bin audit."""
import json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import multiprocessing as mp
import numpy as np
from shapely.geometry import Polygon
from adaptive_coupled_mass import atomic
from common_interface_contact import initialize
from common_contact_evolution import reference_moments
from predictive_boundary_density import BoundaryDensity
from predictive_evolution_batch import density_options
from predictive_coupling_batch import metrics
from predictive_resolution_audit import trace_error
import moment_ale_evolve as model
import kinetic_factor_rhs as kinetic


def validate_cohort_time(actual,expected):
    if not np.isfinite(expected) or expected<0 or not np.isfinite(actual) or abs(actual-expected)>1e-15:
        raise RuntimeError('cohort times differ')


def inspect(entry):
    started=time.monotonic();initial=json.loads(Path(entry['initial']).read_text());cfg=initial['config']
    cp=json.loads(Path(entry['checkpoint']).read_text());mass,_,refs,phases=initialize(cfg['n'],cfg['amplitude'],cfg['angle'])
    np.testing.assert_array_equal(cp['mass'],mass)
    validate_cohort_time(cp['t'],float(entry.get('expected_time',1e-5)))
    options=density_options(cfg,True);options.setdefault('refine_moment_direction',False)
    if cfg['kind']!='common_planar_contact' or options['mode']!='moment' or not options['align_core'] or options['core_level']!=2 or options['phase_fraction_floor']!=1e-7 or options['refine_moment_direction'] or options['radial']!=1e-7 or options['segments']!=64 or not options['smooth_phase_floor']:
        raise RuntimeError('cohort reconstruction differs')
    y=np.asarray(cp['state']);a=y[:-4].reshape(-1,20);L=[np.asarray(f) for f in cp['factors']]
    cells,edges=model.geometry(a[:,:2],a[:,2],cp['t']);raw=np.column_stack([mass,a[:,3:]])
    rows=[];coupling=None
    for segments in (64,256):
        provider=BoundaryDensity(**dict(options,segments=segments));_,closure=provider(cells,edges,raw)
        if segments==64:
            args=dict(covariance_factors=L,face_density_provider=provider,face_segments=64,exact_bin_transport=True)
            dy2,d2=kinetic.rhs(cp['t'],y,mass,face_order=2,**args)
            dy4,d4=kinetic.rhs(cp['t'],y,mass,face_order=4,**args)
            coupling=dict(order2_vs4_difference=float(np.max(abs(dy2-dy4)/np.maximum(1.,abs(dy4)))),
                order2=metrics(dy2,d2,mass),order4=metrics(dy4,d4,mass))
        exact,coverage,advected=reference_moments(cells,refs,phases,cp['t']);high=Polygon(advected[1]);error=0.
        for parts,rhos in provider.last_cells:
            for triangle,rho in zip(parts,rhos):
                p=Polygon(triangle);Ah=p.intersection(high).area
                error+=abs(rho-1.)*max(0.,p.area-Ah)+abs(rho-cfg['amplitude'])*Ah
        volume=error/float(mass.sum());bound=float(np.sum(abs(mass-exact[:,0]))/mass.sum());face=trace_error(provider,edges,advected,phases)
        global_error=float(abs(exact[:,0].sum()-mass.sum())/mass.sum())
        valid=coverage<1e-11 and global_error<1e-11 and volume+1e-8>=bound and closure['closure_moment_error']<1e-8 and closure['closure_min_density']>0 and closure['closure_volume_cap_excess']<1e-8
        rows.append(dict(segments=segments,reference_valid=bool(valid),density_L1_relative=volume,
            cell_average_L1_lower_bound=bound,face_pointwise_L1_relative=face['relative_L1'],
            face_absolute_L1=face['absolute_L1'],face_normalization=face['normalization'],closure=closure))
        print('fixed time',entry['name'],segments,'volume',volume,'face',face['relative_L1'],flush=True)
    compatible=coupling['order2_vs4_difference']<1e-10 and all(m['mass_flux']<1e-10 and m['pressure_rate']<1e-8 and max(m['affine_rhs_errors'].values())<1e-8 for m in (coupling['order2'],coupling['order4']))
    result=dict(case=entry['name'],n=cfg['n'],angle=cfg['angle'],t=cp['t'],passed=bool(compatible and all(r['reference_valid'] for r in rows)),
        coupling=coupling,rows=rows,seconds=time.monotonic()-started)
    atomic(entry['name']+'_comparison.json',result);return result


def run():
    cfg=json.loads(Path('fixed_time_config.json').read_text())
    if not json.loads(Path(cfg['completion_gate']).read_text())['passed']:raise RuntimeError('N8 continuation incomplete/failed')
    results=[]
    with ProcessPoolExecutor(max_workers=2,mp_context=mp.get_context('spawn')) as executor:
        futures=[executor.submit(inspect,dict(entry,expected_time=cfg.get('expected_time',1e-5))) for entry in cfg['cases']]
        for f in as_completed(futures):
            results.append(f.result());atomic('fixed_time_progress.json',dict(results=results))
    report=dict(completed=True,passed=all(r['passed'] for r in results),results=results,
        expected_time=cfg.get('expected_time',1e-5),
        note='Same physical planar IC and specified common time across N3/N4/N8. No oracle in RHS/reconstruction. S256 is a read-only reconstruction of the SAME S64-evolved moments, not an S256 evolution. Gauss2/4 endpoint compatibility checked. All-face normalization varies with n; table is not a spatial convergence order. Reference/compatibility gates are not an accuracy threshold.')
    atomic('fixed_time_comparison_report.json',report);raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()
