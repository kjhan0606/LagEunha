"""Physical-error-controlled nonlinear Euler simple-wave moment benchmark.

Fixed numerical mass, material generators, positive midpoint covariance factors.
This benchmark explicitly allocates time/path errors below the physical error
goal. It does not change the older strict affine-contact path gates.
"""
import hashlib,json,time,os,multiprocessing as mp
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from adaptive_coupled_mass import atomic,simpson,restore
from adaptive_moment_batch import setup,D0,diagnostic
from factor_moment_evolve import rebuild,advance_factors_sylvester as advance
from coupled_factor_path_audit import dense_state
from curved_contact_first_step import load_saved_rhs
from smooth_mass_path import evaluate_smooth_path,sample_record,sample_record_file,initialize_worker
from simple_wave_exact import physical_raw
from saved_curved_contact import checked_json
from smooth_path_quadrature import node_indices,assess_values,requires_refinement
import moment_ale_evolve as model
import kinetic_factor_rhs as kinetic


def pressure_diagnostics(diag):
    return {k:v for k,v in diag.items() if k.startswith('pressure_') and (isinstance(v,str) or np.ndim(v)==0)}


def update_pressure_peak(peak,rows):
    if not any(rows):return
    peak['pressure_diagnostics_steps_observed']=peak.get('pressure_diagnostics_steps_observed',0)+1
    for row in rows:
        for key,value in row.items():
            if key.endswith('_limited_cells'):peak['pressure_limited_cells_max']=max(peak.get('pressure_limited_cells_max',0),int(value))
            if key.endswith('_min_positivity_scale'):peak['pressure_positivity_scale_min']=min(peak.get('pressure_positivity_scale_min',1.),float(value))
            if key.endswith('_neighbor_count_deficient_cells'):peak['pressure_deficient_stencil_cells_max']=max(peak.get('pressure_deficient_stencil_cells_max',0),int(value))


def choose_next_dt(cfg,dt,max_dt,checked):
    if cfg.get('fixed_dt',False):return max_dt
    ratio=checked['physical_ratio'] if cfg.get('geometry_step_controller')=='physical_path' else checked['ratio']
    return min(max_dt,dt*float(np.clip(.9/np.sqrt(max(ratio,1e-12)),.5,1.5)))


def complete_geometry_trial(trial,mass,cfg):
    from geometry_mass_lift import integrate,lifted_state
    record=restore(trial['record'])
    record['geometry_lift']=integrate(record,mass,rtol=cfg.get('geometry_rtol',1e-12),atol=cfg.get('geometry_atol',1e-14),root_tolerance=cfg.get('geometry_mass_root_tolerance',2e-13))
    endpoint,_,_=lifted_state(record,mass,1.)
    new=np.asarray(trial['new']);new[:-4].reshape(-1,20)[:,:3]=endpoint[:-4].reshape(-1,20)[:,:3]
    trial.update(record=record,new=new);trial.pop('geometry_pending',None)
    return trial


def install_domain():
    model.domain=lambda t:D0.copy()
    model.boundary_velocity=lambda p:np.zeros(2)


def preserved_initial(cfg):
    source=checked_json(cfg['preserved_initial_state'],cfg['preserved_initial_state_sha256'])
    expected=dict(kind=cfg.get('problem','nonlinear_simple_wave'),amplitude=cfg['amplitude'],n=cfg['n'])
    if source.get('physical_initial_condition')!=expected:
        if 'preserved_source_config' not in cfg:raise ValueError('preserved wave source needs physical identity')
        old=checked_json(cfg['preserved_source_config'],cfg['preserved_source_config_sha256'])
        if old['n']!=cfg['n'] or old['amplitude']!=cfg['amplitude'] or source.get('config_sha256')!=cfg['preserved_source_config_sha256']:
            raise ValueError('preserved wave source configuration mismatch')
    mass=np.asarray(source['mass']);state=np.asarray(source['state']);L=[np.asarray(v) for v in source['factors']]
    if mass.shape!=(cfg['n']**2,) or state.shape!=(20*len(mass)+4,) or len(L)!=len(mass):
        raise ValueError('preserved wave source shape mismatch')
    if np.any(mass<=0) or not np.all(np.isfinite(mass)) or not np.all(np.isfinite(state)) or not 0<=source['t']<cfg['end']:
        raise ValueError('invalid preserved wave source')
    return source,mass,state,L


def finish_pending_commit(pending_path,checkpoint_path,sample_path,identity):
    """Complete a saved decision; never recompute an approved dynamical step."""
    if not pending_path.exists():return
    trial=json.loads(pending_path.read_text())
    if 'commit' not in trial:return
    commit=trial['commit'];cp=commit['checkpoint']
    if cp['config_sha256']!=identity:raise ValueError('pending commit configuration mismatch')
    if checkpoint_path.exists():
        existing=json.loads(checkpoint_path.read_text())
        if existing!=cp and (existing['t']!=trial['record']['t'] or not np.array_equal(existing['state'],trial['record']['y'])):
            raise ValueError('pending commit does not extend the saved checkpoint')
    atomic(checkpoint_path,cp)
    if commit['accepted']:
        rhs_path=pending_path.parent/'source_rhs.json'
        if rhs_path.exists():
            rhs=json.loads(rhs_path.read_text())
            expected=hashlib.sha256(np.asarray(trial['record']['y']).tobytes()).hexdigest()
            if rhs['state_sha256']!=expected:raise ValueError('committed source RHS identity mismatch')
            rhs_path.unlink()
    sample_archive=pending_path.parent/commit['sample_archive']
    if sample_path.exists():sample_path.replace(sample_archive)
    pending_path.replace(pending_path.parent/commit['archive'])


def initial(n,amplitude):
    _,skeleton=setup('wave',n,0.);a=skeleton[:-4].reshape(-1,20)
    cells,_=model.geometry(a[:,:2],a[:,2],0.);masses=[];factors=[]
    for i,poly in enumerate(cells):
        raw,points,weights,rho,velocity,pressure=physical_raw(poly,0.,amplitude)
        mass=raw[0];c=raw[1:3]/mass;mean=raw[7:9]/mass
        S=raw[3:7].reshape(2,2)-np.outer(raw[1:3],raw[1:3])/mass
        C=raw[9:13].reshape(2,2)-np.outer(raw[1:3],raw[7:9])/mass
        B=np.linalg.solve(S,C).T;residual=velocity-mean-(points-c)@B.T
        _,R=np.linalg.qr(np.sqrt(weights*rho)[:,None]*residual,mode='reduced')
        factors.append(R.T);masses.append(mass);a[i,3:]=raw[1:]
    mass=np.asarray(masses);state=rebuild(skeleton,mass,factors)
    adjustment=float(np.max(abs(state-skeleton)))
    return mass,state,factors,adjustment


def score(t,state,mass,amplitude,factors=None):
    a=np.asarray(state)[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],t)
    states=[model.decode(poly,np.r_[m,row[3:]]) for m,row,poly in zip(mass,a,cells)]
    quadratic=os.environ.get('LAGFORCE_QUADRATIC_VELOCITY','0')=='1'
    if quadratic:
        if factors is None:raise ValueError('quadratic velocity scoring requires the primary factors')
        from quadratic_velocity_reconstruction import prepare,correction
        prepare(states,cells,edges,mass,factors)
    centers,grad=kinetic.pressure_gradients(states,edges)
    from quadratic_pressure_reconstruction import enabled as quadratic_pressure_enabled
    pressure_polynomials=None
    if quadratic_pressure_enabled():
        from quadratic_pressure_reconstruction import prepare as prepare_pressure,value as pressure_value
        pressure_polynomials,_=prepare_pressure(states,cells)
    errors=dict(pressure_cell=0.,velocity_cell=0.,density_cell=0.,pressure_field=0.,velocity_field=0.,density_field=0.,transverse_velocity=0.,energy_cell=0.)
    exact_mass=0.;reference_quad=0.;minrho=float('inf')
    for i,(m,row,poly,s) in enumerate(zip(mass,a,cells,states)):
        raw,points,weights,rho,velocity,pressure=physical_raw(poly,t,amplitude,12)
        fine,*_=physical_raw(poly,t,amplitude,16)
        reference_quad+=float(np.max(abs(fine-raw)))/amplitude
        area=weights.sum();exact_mass+=raw[0]
        errors['pressure_cell']+=abs(area*s['P']-weights@pressure)/amplitude
        errors['velocity_cell']+=area*np.linalg.norm(s['v']-raw[7:9]/raw[0])/amplitude
        errors['density_cell']+=abs(m-raw[0])/amplitude
        errors['energy_cell']+=abs(row[19]-raw[17])/amplitude
        numerical_velocity=s['v']+(points-s['c'])@s['B'].T
        if quadratic:numerical_velocity+=correction(s,points)
        numerical_pressure=s['P']+(points-centers[i])@grad[i]
        if pressure_polynomials is not None:numerical_pressure=pressure_value(pressure_polynomials[i],points)
        numerical_density=s['density'].value(points);minrho=min(minrho,float(np.min(numerical_density)))
        errors['velocity_field']+=weights@np.linalg.norm(numerical_velocity-velocity,axis=1)/amplitude
        errors['pressure_field']+=weights@abs(numerical_pressure-pressure)/amplitude
        errors['density_field']+=weights@abs(numerical_density-rho)/amplitude
        errors['transverse_velocity']+=weights@abs(numerical_velocity[:,1])/amplitude
    return dict(t=t,amplitude=amplitude,errors_over_amplitude={k:float(v) for k,v in errors.items()},
        reference_quadrature_difference=reference_quad,reference_global_mass_relative=float(abs(exact_mass-mass.sum())/mass.sum()),
        density_min=minrho,reference_valid=bool(reference_quad<1e-5 and abs(exact_mass-mass.sum())/mass.sum()<1e-10))


def run_with_executor(executor=None):
    started=time.monotonic();data=Path('simple_wave_config.json').read_bytes();cfg=json.loads(data)
    pressure_solver=os.environ.get('LAGFORCE_PRESSURE_RIEMANN','hllc_davis')
    if pressure_solver not in ('hllc_davis','exact_euler_star') or cfg.get('pressure_riemann_solver','hllc_davis')!=pressure_solver:
        raise ValueError('pressure Riemann solver configuration mismatch')
    backend='vector_sparse' if os.environ.get('LAGFORCE_VECTOR_KINETIC_ROOT','0')=='1' else 'legacy_dense'
    if cfg.get('mass_root_backend','legacy_dense')!=backend:raise ValueError('mass root backend configuration mismatch')
    velocity_closure='affine'
    if os.environ.get('LAGFORCE_QUADRATIC_VELOCITY','0')=='1':
        from quadratic_velocity_reconstruction import closure_name
        velocity_closure=closure_name()
    if cfg.get('velocity_reconstruction','affine')!=velocity_closure:raise ValueError('velocity reconstruction configuration mismatch')
    from quadratic_velocity_reconstruction import fit_mode
    velocity_fit=fit_mode()
    if cfg.get('velocity_fit','radial_fixed_mean')!=velocity_fit:raise ValueError('velocity fit configuration mismatch')
    if velocity_fit!='radial_fixed_mean' and velocity_closure=='affine':raise ValueError('velocity fit requires quadratic reconstruction')
    from quadratic_velocity_reconstruction import stencil_mode
    velocity_stencil=stencil_mode()
    if cfg.get('velocity_stencil','area')!=velocity_stencil:raise ValueError('velocity stencil configuration mismatch')
    if velocity_stencil!='area' and not velocity_fit.startswith('tensor_free_'):raise ValueError('cell-width velocity stencil requires tensor fit')
    if cfg.get('event_aligned_retries',False) and cfg.get('geometry_path','quadratic_midpoint')!='quadratic_midpoint':
        raise ValueError('event-aligned retries require quadratic midpoint geometry')
    from quadratic_velocity_reconstruction import budget_backend_name
    if cfg.get('quadratic_budget_backend','scalar')!=budget_backend_name():raise ValueError('quadratic budget backend configuration mismatch')
    from quadratic_pressure_reconstruction import closure_name as pressure_closure_name
    pressure_closure=pressure_closure_name()
    if cfg.get('pressure_reconstruction','linear')!=pressure_closure:raise ValueError('pressure reconstruction configuration mismatch')
    pressure_reflection=cfg.get('pressure_reflection_y','')
    if pressure_reflection!=os.environ.get('LAGFORCE_PRESSURE_REFLECT_Y',''):raise ValueError('pressure reflection configuration mismatch')
    if pressure_reflection and pressure_closure=='linear':raise ValueError('pressure reflection requires polynomial reconstruction')
    from pressure_fit_stencil import mode as pressure_fit_mode
    pressure_fit=pressure_fit_mode()
    if cfg.get('pressure_fit','radial_fixed_mean')!=pressure_fit:raise ValueError('pressure fit configuration mismatch')
    if pressure_fit!='radial_fixed_mean' and pressure_closure=='linear':raise ValueError('pressure fit requires polynomial reconstruction')
    from polynomial_pressure_reconstruction import evaluation_backend_name
    if cfg.get('pressure_evaluation_backend','monomial')!=evaluation_backend_name():raise ValueError('pressure evaluation backend configuration mismatch')
    from smooth_mass_path import evaluation_backend_name as path_evaluation_backend
    if cfg.get('mass_path_evaluation_backend','scalar')!=path_evaluation_backend():raise ValueError('mass path evaluation backend configuration mismatch')
    from vector_moment_flux import backend_name as moment_flux_backend
    if cfg.get('moment_flux_evaluation_backend','scalar')!=moment_flux_backend():raise ValueError('moment flux evaluation backend configuration mismatch')
    active=os.environ.get('LAGFORCE_ACTIVE_POWER_CELLS','0')=='1'
    geometry_backend='seeded_active' if os.environ.get('LAGFORCE_SEEDED_POWER_CELLS','0')=='1' else ('active' if active else 'exhaustive')
    if cfg.get('power_cell_evaluation_backend','active' if active else 'exhaustive')!=geometry_backend:raise ValueError('power cell evaluation backend configuration mismatch')
    if cfg.get('geometry_path','quadratic_midpoint') not in ('quadratic_midpoint','material_mass_lift'):raise ValueError('unknown geometry path')
    if cfg.get('geometry_step_controller','combined') not in ('combined','physical_path'):raise ValueError('unknown geometry step controller')
    if cfg.get('geometry_step_controller')=='physical_path' and cfg.get('geometry_path')!='material_mass_lift':raise ValueError('physical path controller needs geometry ODE')
    initial_function=initial;score_function=score
    problem=cfg.get('problem','nonlinear_simple_wave')
    if problem=='steady_vortex':
        from steady_vortex_evolution import initial as initial_function,score as score_function
    elif problem=='sod_gamma5_3':
        from sod_evolution import initial as initial_function,score as score_function
    elif problem!='nonlinear_simple_wave':raise ValueError('unknown physical benchmark')
    from kinetic_speed_smoothing import parameters
    smoothing,smoothing_length=parameters()
    if cfg.get('kinetic_smoothing_coefficient',0.)!=smoothing or cfg.get('kinetic_smoothing_length',1.)!=smoothing_length:raise ValueError('kinetic smoothing configuration mismatch')
    identity=hashlib.sha256(data).hexdigest();end=float(cfg['end']);max_dt=float(cfg['max_dt']);amp=cfg['amplitude']
    install_domain();input_path=Path('initial_state.json')
    if input_path.exists():
        original=json.loads(input_path.read_text())
        if original['config_sha256']!=identity:raise ValueError('simple-wave source identity mismatch')
        mass=np.asarray(original['mass']);y0=np.asarray(original['state']);L0=[np.asarray(v) for v in original['factors']]
    else:
        if cfg.get('preserved_initial_state'):
            source,mass,y0,L0=preserved_initial(cfg)
            original=dict(t=source['t'],mass=mass,state=y0,factors=L0,config_sha256=identity,
                provenance=dict(source=cfg['preserved_initial_state'],sha256=cfg['preserved_initial_state_sha256']))
        else:
            mass,y0,L0,adjustment=initial_function(cfg['n'],amp)
            original=dict(t=0.,mass=mass,state=y0,factors=L0,config_sha256=identity,roundoff_Pi_completion=adjustment)
        atomic(input_path,original);atomic('initial_score.json',score_function(original['t'],y0,mass,amp,L0))
    from physical_gate_accounting import budget_duration
    execution_duration=end-original['t']
    duration=budget_duration(cfg,original['t'])
    checkpoint=Path('checkpoint.json');pending_path=Path('pending.json');sample_path=Path('pending_samples.json')
    finish_pending_commit(pending_path,checkpoint,sample_path,identity)
    if checkpoint.exists():
        cp=json.loads(checkpoint.read_text())
        if cp['config_sha256']!=identity or not np.array_equal(cp['mass'],mass):raise ValueError('simple-wave continuation identity mismatch')
        t=cp['t'];y=np.asarray(cp['state']);L=[np.asarray(v) for v in cp['factors']];dt=cp['next_dt'];steps=cp['steps'];rejected=cp['rejected']
        keys=('absolute','signed','quad_delta','quad_error_bound','absolute_quad_error_bound') if cfg.get('adaptive_path_quadrature') else ('absolute','signed','quad_delta')
        accum={k:np.asarray(cp[k]) for k in keys};peak=cp['peak']
    else:
        t=original['t'];y=y0.copy();L=L0;dt=max_dt;steps=rejected=0
        keys=('absolute','signed','quad_delta','quad_error_bound','absolute_quad_error_bound') if cfg.get('adaptive_path_quadrature') else ('absolute','signed','quad_delta')
        accum={k:np.zeros(len(mass)) for k in keys}
        peak=dict(rate=0.,generator=0.,ledger=0.,stage_mass=0.,min_pressure=1e300,min_factor_R=0.)
    def commit_trial(accepted,archive,sample_archive):
        cp=dict(config_sha256=identity,t=t,next_dt=dt,steps=steps,rejected=rejected,
            mass=mass,state=y,factors=L,peak=peak,**accum)
        trial['commit']=dict(checkpoint=cp,accepted=accepted,archive=archive,sample_archive=sample_archive)
        atomic(pending_path,trial)
        finish_pending_commit(pending_path,checkpoint,sample_path,identity)
    while t<end-1e-15:
        dt=min(dt,end-t)
        if dt<execution_duration*1e-8:raise RuntimeError('simple-wave time controller exhausted')
        if pending_path.exists():
            trial=json.loads(pending_path.read_text())
            if trial['record']['t']!=t or trial['record']['dt']!=dt:raise ValueError('simple-wave pending mismatch')
        else:
            base_path=Path('source_rhs.json');digest=hashlib.sha256(y.tobytes()).hexdigest()
            if base_path.exists():
                base=json.loads(base_path.read_text())
                if base['state_sha256']!=digest or base['t']!=t:raise ValueError('saved source RHS mismatch')
                if base['diag'].get('velocity_reconstruction','affine')!=velocity_closure:raise ValueError('saved source RHS velocity closure mismatch')
                if base['diag'].get('velocity_fit','radial_fixed_mean')!=velocity_fit:raise ValueError('saved source RHS velocity fit mismatch')
                if base['diag'].get('velocity_stencil','area')!=velocity_stencil:raise ValueError('saved source RHS velocity stencil mismatch')
                if base['diag'].get('pressure_reconstruction','linear')!=pressure_closure:raise ValueError('saved source RHS pressure closure mismatch')
                if base['diag'].get('pressure_reflection_y','')!=pressure_reflection:raise ValueError('saved source RHS pressure reflection mismatch')
                if base['diag'].get('pressure_fit','radial_fixed_mean')!=pressure_fit:raise ValueError('saved source RHS pressure fit mismatch')
                if base['diag'].get('pressure_riemann_solver','hllc_davis')!=pressure_solver:raise ValueError('saved source RHS pressure Riemann solver mismatch')
                if base['diag'].get('kinetic_smoothing_coefficient',0.)!=smoothing or base['diag'].get('kinetic_smoothing_length',1.)!=smoothing_length:raise ValueError('saved RHS smoothing mismatch')
                loaded=load_saved_rhs(base);d0,diag0=loaded['dy'],loaded['diag']
            else:
                if t==original['t'] and cfg.get('preserved_source_rhs'):
                    base=checked_json(cfg['preserved_source_rhs'],cfg['preserved_source_rhs_sha256'])
                    if base['state_sha256']!=digest or base['t']!=t:raise ValueError('preserved source RHS mismatch')
                    if base['diag'].get('velocity_reconstruction','affine')!=velocity_closure:raise ValueError('preserved source RHS velocity closure mismatch')
                    if base['diag'].get('velocity_fit','radial_fixed_mean')!=velocity_fit:raise ValueError('preserved source RHS velocity fit mismatch')
                    if base['diag'].get('velocity_stencil','area')!=velocity_stencil:raise ValueError('preserved source RHS velocity stencil mismatch')
                    if base['diag'].get('pressure_reconstruction','linear')!=pressure_closure:raise ValueError('preserved source RHS pressure closure mismatch')
                    if base['diag'].get('pressure_reflection_y','')!=pressure_reflection:raise ValueError('preserved source RHS pressure reflection mismatch')
                    if base['diag'].get('pressure_fit','radial_fixed_mean')!=pressure_fit:raise ValueError('preserved source RHS pressure fit mismatch')
                    if base['diag'].get('pressure_riemann_solver','hllc_davis')!=pressure_solver:raise ValueError('preserved source RHS pressure Riemann solver mismatch')
                    if 'mass_sha256' in base and base['mass_sha256']!=hashlib.sha256(mass.tobytes()).hexdigest():
                        raise ValueError('preserved source RHS mass mismatch')
                    if base['diag'].get('kinetic_smoothing_coefficient',0.)!=smoothing or base['diag'].get('kinetic_smoothing_length',1.)!=smoothing_length:raise ValueError('preserved RHS smoothing mismatch')
                    loaded=load_saved_rhs(base);d0,diag0=loaded['dy'],loaded['diag']
                    atomic(base_path,base)
                else:
                    d0,diag0=kinetic.rhs(t,y,mass,covariance_factors=L,face_order=4)
                    atomic(base_path,dict(state_sha256=digest,t=t,dy=d0,diag=diag0))
            halfL=advance(L,diag0,dt/2);half=rebuild(y+dt*d0/2,mass,halfL)
            dm,diag=kinetic.rhs(t+dt/2,half,mass,covariance_factors=halfL,face_order=4)
            newL=advance(L,diag,dt);new=rebuild(y+dt*dm,mass,newL)
            keys=('cov_linearization_factors','cov_A','cov_K','cov_Qfactor')
            record=dict(t=t,dt=dt,y=y,L=L,d0=d0,dm=dm,diag={k:diag[k] for k in keys})
            if cfg.get('dense_moment_arithmetic')=='longdouble':record['dense_moment_arithmetic']='longdouble'
            trial=dict(record=record,new=new,newL=newL,stage_mass=max(diag0['mass_flux'],diag['mass_flux']))
            trial['pressure_reconstruction_diagnostics']=[pressure_diagnostics(diag0),pressure_diagnostics(diag)]
            if cfg.get('adaptive_path_quadrature'):trial['sample_lattice']=cfg.get('path_sample_lattice',128)
            if cfg.get('geometry_path','quadratic_midpoint')=='material_mass_lift':
                trial['geometry_pending']=True;atomic(pending_path,trial)
                trial=complete_geometry_trial(trial,mass,cfg)
            atomic(pending_path,trial)
        if trial.get('geometry_pending'):
            trial=complete_geometry_trial(trial,mass,cfg);atomic(pending_path,trial)
        record=restore(trial['record']);samples=json.loads(sample_path.read_text()) if sample_path.exists() else {}
        adaptive=bool(cfg.get('adaptive_path_quadrature',False))
        lattice=cfg.get('path_sample_lattice',128)
        if adaptive:
            from smooth_path_quadrature import ensure_sample_lattice
            trial,samples=ensure_sample_lattice(trial,pending_path,sample_path,lattice)
        counts=tuple(v for v in (9,17,33,65,129,257,513,1025) if v-1<=lattice) if adaptive else (9,)
        for count in counts:
            indices=node_indices(count,lattice) if adaptive else list(range(9));divisor=lattice if adaptive else 8
            missing=[j for j in indices if str(j) not in samples]
            if executor is None:results=((j,sample_record(record,mass,j/divisor)) for j in missing)
            else:
                from predictive_evolution_batch import record_digest
                digest=record_digest(trial['record'])
                futures={executor.submit(sample_record_file,str(pending_path.resolve()),digest,mass,j/divisor):j for j in missing}
                results=((futures[f],f.result()) for f in as_completed(futures))
            unsaved=0
            try:
                for j,value in results:
                    samples[str(j)]=value;unsaved+=1
                    if unsaved>=cfg.get('path_save_batch',1):atomic(sample_path,samples);unsaved=0
            finally:
                if unsaved:atomic(sample_path,samples)
            values=[samples[str(j)] for j in indices]
            checked=assess_values(values,mass,dt,duration,cfg,check_absolute_quadrature=adaptive)
            if not requires_refinement(checked,cfg,adaptive,'geometry_lift' in record):break
            print('simple wave',cfg['n'],'refining',count,'quad ratio',checked['quadrature_ratio'],flush=True)
        signed=checked['signed'];absolute=checked['absolute'];quad=checked['quad_delta']
        rate=checked['rate'];gen=checked['generator'];pressure=checked['pressure'];gram=checked['gram'];ratio=checked['ratio']
        if not np.isfinite(ratio) or pressure<=0 or gram<-1e-20:raise RuntimeError('simple-wave physical path state failed')
        if ratio>1:
            if cfg.get('fixed_dt',False):raise RuntimeError('fixed time mesh failed path gate: '+str(ratio))
            next_dt=dt/2
            if cfg.get('event_aligned_retries',False):
                from power_flip_time import rejected_step_proposal
                import moment_ale_evolve as event_geometry
                proposal=rejected_step_proposal(record,lambda x,w,time:event_geometry.geometry(x,w,time)[1],
                                               minimum_dt=execution_duration*1e-8)
                next_dt=proposal['dt'];trial['event_retry_proposal']=proposal;atomic(pending_path,trial)
            rejected+=1;dt=next_dt
            commit_trial(False,'rejection_'+str(rejected)+'.json','rejection_'+str(rejected)+'_samples.json')
            print('simple wave',cfg['n'],'rejected',ratio,'dt',dt,flush=True);continue
        new=np.asarray(trial['new']);newL=[np.asarray(v) for v in trial['newL']]
        corrected=new.copy();corrected[-4:]-=y0[-4:]
        d=diagnostic(t+dt,corrected,mass,y0)
        if d['ledger_error']>=1e-10 or d['min_pressure']<=0 or trial['stage_mass']>=1e-10:
            raise RuntimeError('simple-wave accepted state gate failed: '+str(d))
        for k,v in [('signed',signed),('absolute',absolute),('quad_delta',quad)]:accum[k]+=v
        if adaptive:
            accum['quad_error_bound']+=abs(quad);accum['absolute_quad_error_bound']+=abs(checked['absolute_quad_delta'])
        peak['rate']=max(peak['rate'],rate);peak['generator']=max(peak['generator'],gen)
        peak['ledger']=max(peak['ledger'],d['ledger_error']);peak['stage_mass']=max(peak['stage_mass'],trial['stage_mass'])
        peak['min_pressure']=min(peak['min_pressure'],pressure,d['min_pressure']);peak['min_factor_R']=min(peak['min_factor_R'],gram)
        update_pressure_peak(peak,trial.get('pressure_reconstruction_diagnostics',[]))
        steps+=1;t=min(end,t+dt);y,L=new,newL
        atomic('step_'+str(steps)+'_report.json',dict(t=t,dt=dt,ratio=ratio,rate=rate,generator=gen,
            absolute=absolute,signed=signed,quad_delta=quad,diagnostic=d,path_nodes=count,
            absolute_quad_delta=checked['absolute_quad_delta'],quadrature_ratio=checked['quadrature_ratio'],pressure_reconstruction_diagnostics=trial.get('pressure_reconstruction_diagnostics',[])))
        dt=choose_next_dt(cfg,dt,max_dt,checked)
        commit_trial(True,'step_'+str(steps)+'_trial.json','step_'+str(steps)+'_samples.json')
        print('simple wave',cfg['n'],'accepted',steps,'t',t,'ratio',ratio,flush=True)
    final=score_function(t,y,mass,amp,L);errors=final['errors_over_amplitude']
    passed=final['reference_valid'] and all(errors[k]<cfg.get('physical_goals',{}).get(k,cfg['physical_goal']) for k in ('pressure_cell','velocity_cell','density_cell','pressure_field','velocity_field','density_field'))
    from physical_gate_accounting import numerical_gates
    final_numerical_gates=numerical_gates(accum,mass,peak,cfg);numerical=all(final_numerical_gates.values())
    report=dict(completed=True,problem=problem,n=cfg['n'],amplitude=amp,t=t,source_time=original['t'],steps=steps,rejected=rejected,
        validation_end=float(cfg.get('validation_end',end)),path_budget_duration=duration,
        covers_validation_target=bool(t>=float(cfg.get('validation_end',end))-1e-15),
        numerical_passed=bool(numerical),numerical_gate_details=final_numerical_gates,adaptive_path_quadrature=bool(cfg.get('adaptive_path_quadrature',False)),
        physical_accuracy_passed=bool(passed),final_score=final,peak=peak,
        integrated_absolute_mass_defect=float(np.max(accum['absolute']/mass)),quadrature_mass_difference=float(np.max(abs(accum['quad_delta'])/mass)),
        config=cfg,seconds=time.monotonic()-started,
        path_workers=int(os.environ.get('LAGFORCE_SMOOTH_PATH_WORKERS','1')),
        note=problem+' exact Euler reference used only at initialization/scoring. Continuous positive density from numerical six moments; current Sylvester positive midpoint factors. Fixed original mass, material generators. New physical-error allocation for this benchmark; old5e-11 contact budget unchanged. Sampled dense-path quadrature evidence, not a continuous-time bound.')
    if cfg.get('adaptive_path_quadrature'):
        report['signed_quadrature_absolute_sum']=float(np.max(accum['quad_error_bound']/mass))
        report['absolute_quadrature_absolute_sum']=float(np.max(accum['absolute_quad_error_bound']/mass))
    atomic('simple_wave_report.json',report)
    print('simple wave complete',cfg['n'],report['numerical_passed'],report['physical_accuracy_passed'],flush=True)
    raise SystemExit(0 if report['numerical_passed'] else 2)


def run():
    workers=int(os.environ.get('LAGFORCE_SMOOTH_PATH_WORKERS','1'))
    if workers<1:raise ValueError('positive path worker count required')
    if workers==1:return run_with_executor()
    with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize_worker) as executor:
        return run_with_executor(executor)


if __name__=='__main__':run()
