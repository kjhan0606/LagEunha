"""Common-contact continuation from explicit accepted checkpoints."""
import json,hashlib,time,math
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import multiprocessing as mp
import numpy as np
from shapely.geometry import Polygon
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import diagnostic
from predictive_evolution_batch import initialize_case,density_options,advance_interval
from predictive_boundary_density import BoundaryDensity
from predictive_coupling_batch import metrics
from common_interface_contact import initialize
from common_contact_evolution import reference_moments
from predictive_resolution_audit import trace_error
from shear_contact_check import affine_errors
import kinetic_factor_rhs as kinetic
import moment_ale_evolve as model
from preserved_experiment_input import install_input


def valid_metric(m):
    return m['mass_flux']<1e-10 and m['closure_moment_error']<1e-8 and m['density_min']>0 and m['closure_volume_cap_excess']<1e-8 and m['pressure_rate']<1e-8 and max(m['affine_rhs_errors'].values())<1e-8 and m['ledger_rate_error']<1e-10 and m['generator_velocity_defect']<1e-8


def score(cp,cfg,mass,label):
    _,_,refs,phases=initialize(cfg['n'],cfg['amplitude'],cfg['angle']);t=cp['t'];a=np.asarray(cp['state'])[:-4].reshape(-1,20)
    cells,edges=model.geometry(a[:,:2],a[:,2],t);raw=np.column_stack([mass,a[:,3:]])
    provider=BoundaryDensity(**density_options(cfg,True));_,closure=provider(cells,edges,raw)
    exact,coverage,advected=reference_moments(cells,refs,phases,t);high=Polygon(advected[1]);error=0.
    for parts,rhos in provider.last_cells:
        for triangle,rho in zip(parts,rhos):
            p=Polygon(triangle);A=p.area;Ah=p.intersection(high).area
            error+=abs(rho-1.)*max(0.,A-Ah)+abs(rho-cfg['amplitude'])*Ah
    volume=error/float(mass.sum());bound=float(np.sum(abs(mass-exact[:,0]))/mass.sum());face=trace_error(provider,edges,advected,phases)
    global_error=float(abs(exact[:,0].sum()-mass.sum())/mass.sum())
    valid=coverage<1e-11 and global_error<1e-11 and volume+1e-8>=bound and closure['closure_moment_error']<1e-8 and closure['closure_min_density']>0 and closure['closure_volume_cap_excess']<1e-8
    return dict(label=label,t=t,reference_valid=bool(valid),density_L1_relative=volume,cell_average_L1_lower_bound=bound,
        q_L1_over_total_mass=float(np.sum(abs(raw[:,1:3]-exact[:,1:3]))/mass.sum()),Q_L1_over_total_mass=float(np.sum(abs(raw[:,3:7]-exact[:,3:7]))/mass.sum()),
        face_pointwise_L1_relative=face['relative_L1'],reference_coverage=coverage,exact_global_mass_relative_error=global_error,closure=closure)


def run_case(entry):
    name=entry['case'];source=json.loads(Path(entry['source_copy']).read_text());ic=json.loads(Path(entry['initial_config_copy']).read_text())
    if source['config']['start']!=ic['t'] or source['config']['source_state_sha256']!=hashlib.sha256(np.asarray(ic['state']).tobytes()).hexdigest():
        raise RuntimeError('preserved prefix does not start at the supplied initial state')
    mode=entry.get('mode','moment')
    cfg=dict(ic['config']);cfg['density_provider']=dict(cfg['density_provider'],mode=mode,refine_moment_direction=mode=='global_moment')
    cfg['refined_path_audit']=True;cfg['nonmonotone_step_search']=True
    cfg['adaptive_path_quadrature']=bool(entry.get('adaptive_path_quadrature',False))
    if 'max_path_quadrature_nodes' in entry:cfg['max_path_quadrature_nodes']=int(entry['max_path_quadrature_nodes'])
    mass,_=initialize_case(cfg);np.testing.assert_array_equal(source['mass'],mass);np.testing.assert_array_equal(ic['mass'],mass)
    y=np.asarray(source['state']);L=[np.asarray(f) for f in source['factors']];provider=BoundaryDensity(**density_options(cfg,True))
    args=dict(covariance_factors=L,face_density_provider=provider,face_segments=64,exact_bin_transport=True)
    dy4,d4=kinetic.rhs(source['t'],y,mass,face_order=4,**args)
    dy2,d2=kinetic.rhs(source['t'],y,mass,face_order=2,**args)
    m4=metrics(dy4,d4,mass);m2=metrics(dy2,d2,mass)
    diff=float(np.max(abs(dy4-dy2)/np.maximum(1.,abs(dy4))))
    if not valid_metric(m4):raise RuntimeError('continuation source gate failed: '+name)
    cfg['face_order']=2 if valid_metric(m2) and diff<1e-10 else 4
    gate=dict(case=name,passed=True,face_order=cfg['face_order'],order2_vs4_difference=diff,order4=m4,order2=m2)
    atomic(name+'_source_gate.json',gate)
    path=Path('inputs')/name/'checkpoint.json'
    input_cp=dict(source,config=cfg,provenance=dict(source_copy=entry['source_copy'],source_original=entry['source_original'],
        source_sha256=hashlib.sha256(Path(entry['source_copy']).read_bytes()).hexdigest(),
        note='Accepted state/mass preserved; reconstruction mode is explicit. Face order2 used only after order4 comparison; refined path audit and bounded dt search retain all original gates.'))
    install_input(path,input_cp)
    target=float(entry.get('target',1e-5));short_target=min(1e-5,target);results=[];checkpoints=[]
    max_dt=float(entry.get('max_dt',1.25e-6))
    if source['t']<short_target-1e-15:
        duration=short_target-source['t'];steps=max(1,math.ceil(duration/max_dt))
        result=advance_interval((name,'short',duration,steps));results.append(result)
        if result['passed']:
            checkpoints.append(('short',json.loads((Path('evolution')/name/'short/checkpoint.json').read_text())))
            if target>short_target+1e-15:
                results.append(advance_interval((name,'extended',target-short_target,max(1,math.ceil((target-short_target)/max_dt)))))
                if results[-1]['passed']:checkpoints.append(('extended',json.loads((Path('evolution')/name/'extended/checkpoint.json').read_text())))
    elif source['t']<target-1e-15:
        duration=target-source['t'];steps=max(1,math.ceil(duration/max_dt))
        result=advance_interval((name,'short',duration,steps));results.append(result)
        if result['passed']:checkpoints.append(('target',json.loads((Path('evolution')/name/'short/checkpoint.json').read_text())))
    if not all(r['passed'] for r in results):
        report=dict(case=name,passed=False,results=results,source_gate=gate);atomic(name+'_continuation_report.json',report);return report
    final=checkpoints[-1][1] if checkpoints else source;duration=final['t']-source['t'];cps=[cp for _,cp in checkpoints]
    for cp in cps:np.testing.assert_array_equal(cp['mass'],mass)
    signed=sum((np.asarray(cp['signed']) for cp in cps),np.zeros(len(mass)))
    absolute=sum((np.asarray(cp['absolute']) for cp in cps),np.zeros(len(mass)))
    quad=sum((np.asarray(cp['quad_delta']) for cp in cps),np.zeros(len(mass)))
    rate=max((cp['peak']['rate']*duration/cfg['path_control_horizon'] for cp in cps),default=0.)
    corrected=np.asarray(final['state']);corrected[-4:]-=y[-4:];d=diagnostic(final['t'],corrected,mass,y);affine=affine_errors(np.asarray(final['state']),mass)
    combined=dict(source_time=source['t'],target_time=final['t'],duration=duration,normalized_peak_mass_rate=rate,
        integrated_absolute_mass_defect=float(np.max(absolute/mass)),integrated_signed_mass_defect=float(np.max(abs(signed)/mass)),
        quadrature_mass_difference=float(np.max(abs(quad)/mass)),ledger_error=d['ledger_error'],pressure_error=d['pressure_departure'],affine_errors=affine)
    combined['passed']=bool(rate<1e-8 and combined['integrated_absolute_mass_defect']<1e-10 and combined['quadrature_mass_difference']<1e-11 and d['ledger_error']<1e-10 and d['pressure_departure']<1e-8 and max(affine.values())<1e-8)
    full_duration=final['t']-ic['t'];full_abs=np.asarray(source['absolute'])+absolute;full_quad=abs(np.asarray(source['quad_delta']))+abs(quad)
    full_rate=max([source['peak']['rate']]+[cp['peak']['rate'] for cp in cps])*full_duration/cfg['path_control_horizon']
    full_y=np.asarray(final['state']);full_y[-4:]-=np.asarray(ic['state'])[-4:]
    fd=diagnostic(final['t'],full_y,mass,np.asarray(ic['state']))
    whole=dict(start=ic['t'],end=final['t'],prefix_steps=source['steps'],new_steps=sum(r['steps'] for r in results),
        prefix_same_reconstruction=mode==ic['config']['density_provider']['mode'],normalized_peak_mass_rate=full_rate,
        integrated_absolute_mass_defect=float(np.max(full_abs/mass)),quadrature_error_screen=float(np.max(full_quad/mass)),
        ledger_error=fd['ledger_error'],pressure_error=fd['pressure_departure'])
    whole['passed']=bool(full_rate<1e-8 and whole['integrated_absolute_mass_defect']<1e-10 and whole['quadrature_error_screen']<1e-11 and whole['ledger_error']<1e-10 and whole['pressure_error']<1e-8)
    rows=[score(input_cp,cfg,mass,'source')]+[score(cp,cfg,mass,label) for label,cp in checkpoints]
    report=dict(case=name,passed=combined['passed'] and whole['passed'] and all(r['reference_valid'] for r in rows),results=results,
        source_gate=gate,combined_evolution=combined,whole_history=whole,rows=rows,
        note='Accepted prefix preserved. New steps have257+focused65 audit; prefix retains its existing65-node evidence. Reconstruction mode and any hypothesis change are explicit.')
    atomic(name+'_continuation_report.json',report);return report


def run():
    config=json.loads(Path('continuation_config.json').read_text());audit=json.loads(Path(config['refined_audit']).read_text())
    if not audit['completed']:raise RuntimeError('refined trial audit incomplete')
    for entry in config['cases']:
        mode=entry.get('mode','moment')
        if not next(r for r in audit['results'] if r['mode']==mode)['passed']:raise RuntimeError('refined trial gate failed: '+mode)
    results=[]
    with ProcessPoolExecutor(max_workers=2,mp_context=mp.get_context('spawn')) as executor:
        futures=[executor.submit(run_case,entry) for entry in config['cases']]
        for f in as_completed(futures):
            results.append(f.result());atomic('continuation_progress.json',dict(results=results))
    report=dict(completed=True,passed=all(r['passed'] for r in results),cases={r['case']:r for r in results},
        note='Independent per-case chains. Existing local-moment model is retained unless explicitly changed. Mass fixed and sites unsteered. Method(3), general convergence, Kepler and ITS not claimed complete.')
    atomic('common_continuation_report.json',report);raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()
