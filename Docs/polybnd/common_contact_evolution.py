"""Common physical planar IC: static gates, chained evolution, then scoring.

The exact initial moments define a new experiment at t=0. Reference material
polygons are used only for read-only scoring after the numerical evolution.
"""
import json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import numpy as np
from shapely.geometry import Polygon
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import diagnostic
from common_interface_contact import initialize
from predictive_evolution_batch import advance_interval,density_options
from predictive_boundary_density import BoundaryDensity
from predictive_coupling_batch import metrics
from predictive_resolution_audit import trace_error
from shear_interface_audit import moments
from shear_contact_check import affine_errors
import moment_ale_evolve as model
import kinetic_factor_rhs as kinetic


CASES=[(f'n{n}_angle{angle}',n,float(angle)) for n in (3,4) for angle in (0,23)]


def selected_cases():
    path=Path('experiment_config.json')
    if not path.exists():return CASES
    entries=json.loads(path.read_text()).get('cases')
    if entries is None:return CASES
    cases=[(str(c['name']),int(c['n']),float(c['angle'])) for c in entries]
    if not cases or len({name for name,_,_ in cases})!=len(cases):raise ValueError('empty or duplicate cases')
    if any(n<2 or '/' in name or name in ('.','..') for name,n,_ in cases):raise ValueError('invalid case')
    return cases


def reference_moments(cells,refs,phases,t):
    # Here there are two phase polygons, not one material polygon per cell.
    advected=[]
    for vertices in refs:
        r=vertices.copy();r[:,0]+=t*(.2+1.1*r[:,1]);advected.append(r)
    material=[Polygon(r) for r in advected];values=[];coverage=[]
    for vertices in cells:
        cell=Polygon(vertices);overlaps=[cell.intersection(ref) for ref in material]
        values.append(sum((rho*moments(p) for p,rho in zip(overlaps,phases)),np.zeros(7)))
        coverage.append(abs(sum(p.area for p in overlaps)-cell.area))
    return np.asarray(values),float(max(coverage)),advected


def prepare_case(task):
    name,n,angle=task;started=time.monotonic();path=Path('inputs')/name/'checkpoint.json'
    if path.exists():raise RuntimeError('new common experiment refuses existing input: '+str(path))
    path.parent.mkdir(parents=True,exist_ok=True)
    cfg=dict(kind='common_planar_contact',n=n,amplitude=100.,angle=angle,
        smooth_phase_floor=True,path_control_horizon=1e-4,
        density_provider=dict(mode='moment',segments=64,radial=1e-7,core_level=2,
            align_core=True,construct_prior_witness=True,lp_retry=True))
    experiment=Path('experiment_config.json')
    if experiment.exists():
        selected=json.loads(experiment.read_text())
        cfg['density_provider']['phase_fraction_floor']=float(selected['phase_fraction_floor'])
    mass,y,_,_=initialize(n,100.,angle);L=[np.zeros((2,2)) for _ in mass]
    atomic(path,dict(config=cfg,t=0.,mass=mass,state=y,factors=L,
        provenance='New common physical planar IC at t=0, not a restart of the older cellwise-constant experiment.'))
    try:
        provider=BoundaryDensity(**density_options(cfg,True))
        args=dict(face_density_provider=provider,face_segments=provider.segments,face_order=4)
        dy_full,df=model.rhs(0.,y,mass,**args)
        dy,dk=kinetic.rhs(0.,y,mass,covariance_factors=L,exact_bin_transport=True,**args)
        low,_=kinetic.rhs(0.,y,mass,covariance_factors=L,exact_bin_transport=True,**dict(args,face_order=2))
        difference=float(np.max(abs(dy-low)/np.maximum(1.,abs(dy))))
        vals=dict(full_hllc=metrics(dy_full,df,mass),positive_kinetic=metrics(dy,dk,mass))
        passed=all(v['mass_flux']<1e-10 and v['closure_moment_error']<1e-8 and v['density_min']>0
            and v['closure_volume_cap_excess']<1e-8 and v['pressure_rate']<1e-8
            and max(v['affine_rhs_errors'].values())<1e-8 and v['ledger_rate_error']<1e-10
            and v['generator_velocity_defect']<1e-12 for v in vals.values()) and difference<1e-8
        result=dict(case=name,passed=bool(passed),metrics=vals,split_order2_vs4_difference=difference,
            total_initial_mass=float(mass.sum()),density_provider=density_options(cfg,True),seconds=time.monotonic()-started)
    except Exception as exc:result=dict(case=name,passed=False,error=repr(exc),seconds=time.monotonic()-started)
    atomic('static_'+name+'.json',result);print('static',name,result['passed'],flush=True);return result


def score_case(task):
    name,n,angle=task;source=json.loads((Path('inputs')/name/'checkpoint.json').read_text());cfg=source['config']
    mass,y0,refs,phases=initialize(n,cfg['amplitude'],angle)
    np.testing.assert_array_equal(mass,source['mass']);np.testing.assert_array_equal(y0,source['state'])
    short=json.loads((Path('evolution')/name/'short/checkpoint.json').read_text())
    final=json.loads((Path('evolution')/name/'extended/checkpoint.json').read_text())
    duration=final['t']-source['t'];cps=(short,final)
    absolute=sum(np.asarray(cp['absolute']) for cp in cps);signed=sum(np.asarray(cp['signed']) for cp in cps)
    quad=sum(np.asarray(cp['quad_delta']) for cp in cps)
    rate=max(cp['peak']['rate']*duration/cfg['path_control_horizon'] for cp in cps)
    d=diagnostic(final['t'],np.asarray(final['state']),mass,y0);affine=affine_errors(np.asarray(final['state']),mass)
    combined=dict(duration=duration,normalized_peak_mass_rate=rate,
        integrated_absolute_mass_defect=float(np.max(absolute/mass)),integrated_signed_mass_defect=float(np.max(abs(signed)/mass)),
        quadrature_mass_difference=float(np.max(abs(quad)/mass)),ledger_error=d['ledger_error'],
        affine_errors=affine,pressure_error=d['pressure_departure'],generator_defect=max(cp['peak']['generator'] for cp in cps))
    combined['passed']=bool(rate<1e-8 and combined['integrated_absolute_mass_defect']<1e-10 and combined['quadrature_mass_difference']<1e-11
        and d['ledger_error']<1e-10 and max(affine.values())<1e-8 and d['pressure_departure']<1e-8 and combined['generator_defect']<1e-8)
    rows=[]
    for interval,cp in (('initial',source),('short',short),('extended',final)):
        np.testing.assert_array_equal(cp['mass'],mass);t=cp['t'];a=np.asarray(cp['state'])[:-4].reshape(-1,20)
        cells,edges=model.geometry(a[:,:2],a[:,2],t);raw=np.column_stack([mass,a[:,3:]])
        provider=BoundaryDensity(**density_options(cfg,True));_,closure=provider(cells,edges,raw)
        # No reference quantities are passed to the provider or RHS.
        exact,coverage,advected=reference_moments(cells,refs,phases,t)
        high=Polygon(advected[1]);error=0.
        for parts,rhos in provider.last_cells:
            for triangle,rho in zip(parts,rhos):
                p=Polygon(triangle);A=p.area;Ah=p.intersection(high).area
                error+=abs(rho-1.)*max(0.,A-Ah)+abs(rho-cfg['amplitude'])*Ah
        bound=float(np.sum(abs(mass-exact[:,0]))/mass.sum());volume=error/float(mass.sum())
        global_mass=float(abs(exact[:,0].sum()-mass.sum())/mass.sum());face=trace_error(provider,edges,advected,phases)
        valid=coverage<1e-11 and global_mass<1e-11 and volume+1e-8>=bound and closure['closure_moment_error']<1e-8
        valid=valid and closure['closure_min_density']>0 and closure['closure_volume_cap_excess']<1e-8
        row=dict(interval=interval,t=t,reference_valid=bool(valid),density_L1_relative=volume,
            cell_average_L1_lower_bound=bound,face_pointwise_L1_relative=face['relative_L1'],
            face_absolute_L1=face['absolute_L1'],face_normalization=face['normalization'],
            q_L1_over_total_mass=float(np.sum(abs(raw[:,1:3]-exact[:,1:3]))/mass.sum()),
            Q_L1_over_total_mass=float(np.sum(abs(raw[:,3:7]-exact[:,3:7]))/mass.sum()),
            max_partition_area_error=coverage,exact_global_mass_relative_error=global_mass,
            closure=closure,inferred_priors=provider.last_priors)
        rows.append(row);atomic(name+'_density_progress.json',dict(rows=rows))
        print('score',name,interval,'volume',volume,'face',face['relative_L1'],flush=True)
    result=dict(case=name,passed=combined['passed'] and all(r['reference_valid'] for r in rows),combined_evolution=combined,rows=rows)
    atomic(name+'_density_report.json',result);return result


def run(resume=False):
    cases=selected_cases()
    if resume:
        gate=json.loads(Path('static_report.json').read_text())
        if {r['case'] for r in gate['results']}!={name for name,_,_ in cases}:raise RuntimeError('resume static cases differ')
        for name,n,angle in cases:
            cfg=json.loads((Path('inputs')/name/'checkpoint.json').read_text())['config']
            if (cfg['kind'],cfg['n'],cfg['angle'])!=('common_planar_contact',n,angle):raise RuntimeError('resume input differs')
    else:
        with ProcessPoolExecutor(max_workers=2,mp_context=mp.get_context('spawn')) as executor:
            static=list(executor.map(prepare_case,cases))
        gate=dict(passed=all(r['passed'] for r in static),results=static);atomic('static_report.json',gate)
    if not gate['passed']:raise SystemExit(2)
    all_results=[]
    for interval,duration,steps in (('short',1e-5,1),('extended',1e-4,4)):
        tasks=[(name,interval,duration,steps) for name,_,_ in cases]
        with ProcessPoolExecutor(max_workers=2,mp_context=mp.get_context('spawn')) as executor:
            results=list(executor.map(advance_interval,tasks))
        all_results.extend(results);atomic('evolution_progress.json',dict(results=all_results))
        if not all(r['passed'] for r in results):break
    passed=len(all_results)==2*len(cases) and all(r['passed'] for r in all_results)
    atomic('evolution_report.json',dict(completed=len(all_results)==2*len(cases),passed=passed,results=all_results))
    if not passed:raise SystemExit(2)
    with ProcessPoolExecutor(max_workers=2,mp_context=mp.get_context('spawn')) as executor:
        scored=list(executor.map(score_case,cases))
    report=dict(completed=True,passed=all(r['passed'] for r in scored),cases={r['case']:r for r in scored},
        note='Common physical IC across grids; inferred moment normal/aligned core2 in both evolution and all path workers. Fixed mass and unsteered sites. Positive midpoint/RK2, not RK4. 17/65 prescribed-path sample gates. References only after reconstruction; reference validity is not an accuracy/convergence gate. No spatial/time convergence order or general existence proof claimed.')
    atomic('common_evolution_density_report.json',report)
    raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':
    import sys
    run(resume='--resume' in sys.argv)
