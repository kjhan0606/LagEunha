"""Checkpointed short evolution, automatically chained after static gates.
Numerical-only density, current positive midpoint factors, 17/65 path checks.
"""
import json,time,hashlib,os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import multiprocessing as mp
import numpy as np
from adaptive_coupled_mass import atomic,simpson,budget_ratio,restore,plain
from adaptive_moment_batch import setup,diagnostic
from factor_moment_evolve import rebuild,advance_factors_sylvester as advance
from coupled_factor_path_audit import dense_state
from shear_contact_check import affine_errors
import kinetic_factor_rhs as kinetic
from predictive_boundary_density import BoundaryDensity
from predictive_mass_path import evaluate_path


_path_provider=None
_path_record_key=None
_path_record_data=None


def record_digest(record):
    return hashlib.sha256(json.dumps(plain(record),sort_keys=True,separators=(',',':')).encode()).hexdigest()


def path_sample_file(record_path,expected,mass,index):
    """Load each immutable trial once per worker; tasks send only its digest."""
    global _path_record_key,_path_record_data
    key=(str(Path(record_path).resolve()),expected)
    if key!=_path_record_key:
        record=json.loads(Path(record_path).read_text())['record']
        if record_digest(record)!=expected:raise RuntimeError('path worker record digest mismatch')
        _path_record_data=restore(record);_path_record_key=key
    return path_sample(_path_record_data,mass,index)


def submit_path_samples(executor,record,mass,indices,divisor=1.,record_path=None):
    if record_path is None:
        return {executor.submit(path_sample,record,mass,index/divisor):index for index in indices}
    expected=record_digest(record);path=str(Path(record_path).resolve())
    return {executor.submit(path_sample_file,path,expected,mass,index/divisor):index for index in indices}


def physical_samples_valid(result):
    return result['rate']<=5e-9 and result['generator']<=5e-9 and result['pressure']>0 and result['factor_R']>=-1e-20


def rejected_trial_dt(dt,max_dt,state,enabled=False):
    """Bounded alternative search when repeated halving enters a noise floor."""
    state=dict(state);state['rejections']=state.get('rejections',0)+1
    if state.get('recovering',False):
        index=state['candidate_index']+1
        if index<3:
            state['candidate_index']=index
            return (.8,.6,.4)[index]*max_dt,state
        state['recovering']=False
        return state['return_dt']/2,state
    if enabled and not state.get('used',False) and state['rejections']>=4 and dt<max_dt/16:
        state.update(used=True,recovering=True,candidate_index=0,return_dt=float(dt))
        return .8*max_dt,state
    return dt/2,state


def initialize_case(cfg):
    """Rebuild the original IC to validate masses and install its domain."""
    kind=cfg.get('kind','shear_contact')
    if kind=='saved_curved_contact':
        from saved_curved_contact import initialize_saved
        return initialize_saved(cfg)
    if kind=='common_planar_contact':
        from common_interface_contact import initialize
        mass,y,_,_=initialize(cfg['n'],cfg['amplitude'],cfg['angle'])
        return mass,y
    if kind!='shear_contact':raise ValueError('unsupported evolution IC: '+kind)
    return setup(kind,cfg['n'],cfg['amplitude'])


def density_options(cfg,vector_geometry):
    options=dict(mode='neighbor',segments=64,radial=1e-7)
    options.update(cfg.get('density_provider',{}))
    options['vector_geometry']=vector_geometry
    options['smooth_phase_floor']=bool(cfg.get('smooth_phase_floor',False))
    if options['segments']<1 or int(options['segments'])!=options['segments']:
        raise ValueError('face segments must be a positive integer')
    return options


def initialize_path_worker(n,ratio,vector_geometry,smooth_phase_floor=False,options=None):
    global _path_provider
    # Both supported ICs use the same affine moving domain. Workers need
    # its boundary law, not a new mass array or any exact interior field.
    from adaptive_moment_batch import install_affine_shear_domain
    install_affine_shear_domain()
    _path_provider=BoundaryDensity(**(options if options is not None else
        dict(mode='neighbor',vector_geometry=vector_geometry,smooth_phase_floor=smooth_phase_floor)))


def path_sample(record,mass,index,provider=None):
    r=restore(record);N=len(mass);s=index/64
    state,L,dy=dense_state(r['y'],r['L'],r['d0'],r['dm'],r['diag'],mass,r['dt'],s)
    rates=dy[:-4].reshape(N,20)
    f,d=evaluate_path(r['t']+s*r['dt'],state,L,mass,(rates[:,:2],rates[:,2]),provider if provider is not None else _path_provider)
    if not np.all(np.isfinite(f)):raise RuntimeError('nonfinite path mass rate')
    return dict(mass_rate=f.tolist(),generator=d['generator_defect'],pressure=d['min_pressure'],factor_R=d['min_factor_R'])


def assess(record,mass,duration,provider,nodes,samples,cache_path,label,executor=None,record_path=None):
    r=restore(record);N=len(mass)
    missing=[i for i in range(0,65,64//(nodes-1)) if str(i) not in samples]
    if executor is None:results=((index,path_sample(record,mass,index,provider)) for index in missing)
    else:
        futures=submit_path_samples(executor,record,mass,missing,record_path=record_path)
        results=((futures[future],future.result()) for future in as_completed(futures))
    for index,value in results:
        samples[str(index)]=value
        if len(samples)%8==0 or len(samples)==65:atomic(cache_path,samples)
        if len(samples)%16==0 or len(samples)==65:print(label,'path nodes',len(samples),'of65',flush=True)
    atomic(cache_path,samples)
    values=[samples[str(k)] for k in range(0,65,64//(nodes-1))];flux=np.array([x['mass_rate'] for x in values])
    signed=simpson(flux,r['dt']);coarse=simpson(flux[::2],r['dt']);absolute=simpson(abs(flux),r['dt'])
    result=dict(rate=duration*float(np.max(abs(flux)/mass)),generator=max(x['generator'] for x in values),
        pressure=min(x['pressure'] for x in values),factor_R=min(x['factor_R'] for x in values),
        signed=signed,absolute=absolute,quad_delta=signed-coarse)
    result['ratio']=budget_ratio(result['rate'],float(np.max(absolute/mass)),float(np.max(abs(result['quad_delta'])/mass)),result['generator'],r['dt'],duration,5e-11)
    return result


def advance_interval(task,initial_rhs=None):
    name,interval,duration,initial_steps=task
    root=Path('evolution')/name/interval;root.mkdir(parents=True,exist_ok=True)
    result_path=root/'result.json'
    if result_path.exists():
        result=json.loads(result_path.read_text())
        if result['passed']:return result
        raise RuntimeError('failed interval needs diagnosed revision: '+str(result_path))
    original=json.loads((Path('inputs')/name/'checkpoint.json').read_text());cfg=original['config']
    rebuilt,_=initialize_case(cfg);mass=np.asarray(original['mass'])
    np.testing.assert_allclose(mass,rebuilt,rtol=1e-13,atol=1e-15)
    source=original if interval=='short' else json.loads((Path('evolution')/name/'short'/'checkpoint.json').read_text())
    np.testing.assert_array_equal(mass,np.asarray(source['mass']))
    y0=np.asarray(source['state']);L0=[np.asarray(f) for f in source['factors']];start_t=source['t'];target=start_t+duration
    control_horizon=float(cfg.get('path_control_horizon',duration));smooth_phase_floor=bool(cfg.get('smooth_phase_floor',False))
    face_order=int(cfg.get('face_order',4))
    if face_order not in (2,4):raise ValueError('validated face orders are2 and4')
    if control_horizon<duration:raise ValueError('control horizon shorter than evolution interval')
    vector_geometry=os.environ.get('LAGFORCE_VECTOR_GEOMETRY','0')=='1'
    options=density_options(cfg,vector_geometry)
    config=dict(case=name,interval=interval,duration=duration,start=start_t,end=target,initial_steps=initial_steps,
        mass_budget=5e-11,segments=options['segments'],radial=options['radial'],mode=options['mode'],source_state_sha256=hashlib.sha256(y0.tobytes()).hexdigest())
    if os.environ.get('LAGFORCE_ANALYTIC_CAP','0')=='1':config['cap_inverse']='analytic'
    if os.environ.get('LAGFORCE_CENTERED_MOMENT_FIT','0')=='1':config['moment_fit_frame']='cell_centered'
    if os.environ.get('LAGFORCE_ACTIVE_POWER_CELLS','0')=='1':config['power_cell_backend']='active_all_constraints'
    if 'density_provider' in cfg:config['density_provider']=options
    if 'face_order' in cfg:config['face_order']=face_order
    if cfg.get('kind')=='common_planar_contact':
        config['initial_condition']={k:cfg[k] for k in ('kind','n','amplitude','angle')}
    if cfg.get('kind')=='saved_curved_contact':
        config['initial_condition']={k:cfg[k] for k in ('kind','n','initial_state_sha256','reference_sha256')}
    if smooth_phase_floor:config['smooth_phase_floor']=True
    for key in ('path_control_horizon','initial_dt','refined_path_audit','nonmonotone_step_search','adaptive_path_quadrature','max_path_quadrature_nodes'):
        if key in cfg:config[key]=cfg[key]
    cp_path=root/'checkpoint.json'
    if cp_path.exists():
        cp=json.loads(cp_path.read_text())
        if cp['config']!=config:raise RuntimeError('interval configuration mismatch')
        y=np.asarray(cp['state']);L=[np.asarray(f) for f in cp['factors']];t=cp['t'];dt=cp['next_dt'];steps=cp['steps'];rejected=cp['rejected']
        accum={k:np.asarray(cp[k]) for k in ('signed','absolute','quad_delta')}
        peak=cp['peak']
        search_state=cp.get('search_state',{})
    else:
        y=y0.copy();L=L0;t=start_t;dt=min(duration/initial_steps,float(cfg.get('initial_dt',duration/initial_steps)));steps=rejected=0
        accum={k:np.zeros(len(mass)) for k in ('signed','absolute','quad_delta')}
        peak=dict(rate=0.,generator=0.,ledger=0.,pressure_error=0.,stage_mass=0.,min_pressure=1e300,min_factor_R=0.)
        search_state={}
    refined_audit=bool(cfg.get('refined_path_audit',False))
    nonmonotone_search=bool(cfg.get('nonmonotone_step_search',False))
    adaptive_quadrature=bool(cfg.get('adaptive_path_quadrature',False))
    max_quadrature_nodes=int(cfg.get('max_path_quadrature_nodes',1025))
    if max_quadrature_nodes not in (1025,2049,4097):raise ValueError('unsupported maximum path quadrature nodes')
    if max_quadrature_nodes!=1025 and not adaptive_quadrature:raise ValueError('higher node cap requires adaptive quadrature')
    if nonmonotone_search and not refined_audit:raise ValueError('alternative dt search requires refined path audit')
    if adaptive_quadrature and not refined_audit:raise ValueError('adaptive quadrature requires refined path audit')
    provider=BoundaryDensity(**options);started=time.monotonic();base_cache=None
    if initial_rhs is not None:
        if t!=start_t or steps or rejected or cp_path.exists():
            raise ValueError('initial RHS reuse requires the untouched source state')
        base_cache=initial_rhs
    path_workers=int(os.environ.get('LAGFORCE_PATH_WORKERS','1'))
    record_file_cache=path_workers>1 and os.environ.get('LAGFORCE_PATH_RECORD_CACHE','0')=='1'
    path_executor=ProcessPoolExecutor(max_workers=path_workers,mp_context=mp.get_context('spawn'),initializer=initialize_path_worker,
        initargs=(cfg['n'],cfg['amplitude'],vector_geometry,smooth_phase_floor,options)) if path_workers>1 else None
    pending_path=root/'pending.json';nodes_path=root/'pending_nodes.json';refined_path=root/'pending_refined_nodes.json'
    try:
        while t<target-1e-15:
            dt=min(dt,target-t)
            if dt<duration*1e-8:raise RuntimeError('path controller exhausted step size')
            label=name+'/'+interval+'/step'+str(steps+1)
            if pending_path.exists():
                pending=json.loads(pending_path.read_text())
                if pending['t']<t and abs(pending['t']+pending['dt']-t)<1e-14 and (root/('step_'+str(steps)+'.json')).exists():
                    pending_path.replace(root/('committed_pending_'+str(steps)+'.json'))
                    if nodes_path.exists():nodes_path.replace(root/('committed_nodes_'+str(steps)+'.json'))
                    if refined_path.exists():refined_path.replace(root/('refined_nodes_step_'+str(steps)+'.json'))
                    continue
                if pending['t']!=t or pending['dt']!=dt:raise RuntimeError('pending step mismatch')
                record=pending['record'];new=np.asarray(pending['new']);newL=[np.asarray(x) for x in pending['newL']]
                stage_mass=pending['stage_mass'];samples=json.loads(nodes_path.read_text()) if nodes_path.exists() else {}
            else:
                args=dict(face_density_provider=provider,face_segments=provider.segments,face_order=face_order,exact_bin_transport=True)
                if base_cache is None:
                    # A failure at the accepted state is structural, not repaired by smaller dt.
                    base_cache=kinetic.rhs(t,y,mass,covariance_factors=L,**args)
                d0,diag0=base_cache
                halfL=advance(L,diag0,dt/2);half=rebuild(y+dt*d0/2,mass,halfL)
                dm,diag=kinetic.rhs(t+dt/2,half,mass,covariance_factors=halfL,**args)
                newL=advance(L,diag,dt);new=rebuild(y+dt*dm,mass,newL)
                keys=('cov_linearization_factors','cov_A','cov_K','cov_Qfactor')
                record=dict(t=t,dt=dt,y=y,L=L,d0=d0,dm=dm,diag={k:diag[k] for k in keys})
                stage_mass=max(diag0['mass_flux'],diag['mass_flux'])
                pending=dict(t=t,dt=dt,record=record,new=new,newL=newL,stage_mass=stage_mass)
                atomic(pending_path,pending);samples={}
            record_path=pending_path if record_file_cache else None
            control=assess(record,mass,control_horizon,provider,17,samples,nodes_path,label,path_executor,record_path=record_path)
            checked=control
            if control['ratio']<=1 or adaptive_quadrature and physical_samples_valid(control):
                checked=assess(record,mass,control_horizon,provider,65,samples,nodes_path,label,path_executor,record_path=record_path)
            if refined_audit and (checked['ratio']<=1 or adaptive_quadrature and physical_samples_valid(checked)):
                from predictive_path_refinement import assess_refined,assess_adaptive
                evaluate=assess_adaptive if adaptive_quadrature else assess_refined
                kwargs=dict(max_nodes=max_quadrature_nodes) if adaptive_quadrature else {}
                checked=evaluate(record,mass,control_horizon,provider,samples,refined_path,label,path_executor,record_path=record_path,**kwargs)
            if not np.isfinite(checked['ratio']) or checked['pressure']<=0 or checked['factor_R']<-1e-20:
                raise RuntimeError('dense path finite/pressure/factor gate failed')
            if adaptive_quadrature and checked.get('refined_nodes',0)>=max_quadrature_nodes and checked.get('quadrature_ratio',0)>1:
                raise RuntimeError('path quadrature unresolved at'+str(max_quadrature_nodes)+'; preserve trial for integration-accuracy diagnosis')
            if checked['ratio']>1:
                rejected+=1
                atomic(root/('rejection_'+str(rejected)+'.json'),dict(t=t,dt=dt,ratio=checked['ratio'],record=record,control=control,audit=checked,samples=samples))
                if refined_path.exists():refined_path.replace(root/('rejected_refined_nodes_'+str(rejected)+'.json'))
                pending_path.unlink();nodes_path.unlink(missing_ok=True)
                old_dt=dt;dt,search_state=rejected_trial_dt(dt,min(duration/initial_steps,target-t),search_state,nonmonotone_search)
                if dt>old_dt:
                    atomic(root/('alternative_dt_'+str(rejected)+'.json'),dict(t=t,old_dt=old_dt,new_dt=dt,search_state=search_state,
                        note='A different trial length is proposed. It must pass17/65 plus257/focused65 and all unchanged state/mass gates.'))
                atomic(cp_path,dict(config=config,t=t,next_dt=dt,steps=steps,rejected=rejected,mass=mass,state=y,factors=L,peak=peak,search_state=search_state,**accum))
                print(label,'rejected','ratio',checked['ratio'],'new dt',dt,flush=True);continue
            corrected=new.copy();corrected[-4:]-=y0[-4:]
            d=diagnostic(t+dt,corrected,mass,y0);affine=affine_errors(new,mass)
            if not np.all(np.isfinite(new)) or d['min_pressure']<=0 or d['ledger_error']>=1e-10 or max(affine.values())>=1e-8 or d['pressure_departure']>=1e-8 or stage_mass>=1e-10:
                raise RuntimeError('evolution state/pressure/affine/ledger gate failed: '+str(d))
            for key in accum:accum[key]+=checked[key]
            peak['rate']=max(peak['rate'],checked['rate']);peak['generator']=max(peak['generator'],checked['generator'])
            peak['ledger']=max(peak['ledger'],d['ledger_error']);peak['pressure_error']=max(peak['pressure_error'],d['pressure_departure'])
            peak['stage_mass']=max(peak['stage_mass'],stage_mass);peak['min_pressure']=min(peak['min_pressure'],checked['pressure'],d['min_pressure'])
            peak['min_factor_R']=min(peak['min_factor_R'],checked['factor_R'])
            steps+=1;t=min(target,t+dt);y,L=new,newL
            atomic(root/('step_'+str(steps)+'.json'),dict(t=t,dt=dt,record=record,control=control,audit=checked,affine=affine,ledger=d['ledger_error'],samples=samples))
            factor=float(np.clip(.9/np.sqrt(max(checked['ratio'],1e-12)),.5,1.5))
            dt=min(duration/initial_steps,dt*factor)
            search_state={}
            atomic(cp_path,dict(config=config,t=t,next_dt=dt,steps=steps,rejected=rejected,mass=mass,state=y,factors=L,peak=peak,search_state=search_state,**accum))
            if refined_path.exists():refined_path.replace(root/('refined_nodes_step_'+str(steps)+'.json'))
            pending_path.unlink();nodes_path.unlink(missing_ok=True);base_cache=None
            print(label,'accepted t',t,'ratio',checked['ratio'],flush=True)
        absolute=float(max(accum['absolute']/mass));quadrature=float(max(abs(accum['quad_delta'])/mass))
        state_gate=peak['min_pressure']>0 and peak['min_factor_R']>=-1e-20 and peak['ledger']<1e-10 and peak['pressure_error']<1e-8
        path_gate=peak['rate']<1e-8 and absolute<1e-10 and quadrature<1e-11 and peak['generator']<1e-8
        result=dict(case=name,interval=interval,status='completed',passed=bool(state_gate and path_gate),steps=steps,rejected=rejected,
            vector_geometry=vector_geometry,
            path_workers=path_workers,
            path_record_file_cache=record_file_cache,max_path_quadrature_nodes=max_quadrature_nodes,
            cap_inverse=config.get('cap_inverse','legacy_bisection'),
            moment_fit_frame=config.get('moment_fit_frame','global_raw'),
            power_cell_backend=config.get('power_cell_backend','exhaustive'),
            face_order=face_order,
            refined_path_audit=refined_audit,nonmonotone_step_search=nonmonotone_search,
            adaptive_path_quadrature=adaptive_quadrature,
            smooth_phase_floor=smooth_phase_floor,source_time=start_t,path_control_horizon=control_horizon,
            density_provider=options,
            t=t,duration=duration,integrated_absolute_mass_defect=absolute,integrated_signed_mass_defect=float(max(abs(accum['signed'])/mass)),
            quadrature_mass_difference=quadrature,peak=peak,seconds=time.monotonic()-started,
            note='Continues preserved checkpoint; numerical-only density. Positive midpoint/factors;17 control/65 audit nodes, with257+focused65 when refined_path_audit is enabled. Prescribed path rates audited without weight re-solve; sampled checks, not a continuous proof.')
    except Exception as exc:
        result=dict(case=name,interval=interval,status='failed',passed=False,t=t,steps=steps,error=repr(exc),seconds=time.monotonic()-started)
    finally:
        if path_executor is not None:path_executor.shutdown()
    atomic(result_path,result);print(name,interval,result['status'],result['passed'],flush=True);return result


def run_basic_sequence():
    static=json.loads(Path('static_report.json').read_text())
    if not static['passed']:raise RuntimeError('static coupling gates failed')
    all_results=[]
    # Two workers remain inside this foreground job; one BLAS/OpenMP thread each.
    for interval,duration,steps in (('short',1e-5,1),('extended',1e-4,4)):
        tasks=[(name,interval,duration,steps) for name in ('n3_r10','n3_r100','n4_r100')]
        with ProcessPoolExecutor(max_workers=2,mp_context=mp.get_context('spawn')) as executor:
            results=list(executor.map(advance_interval,tasks))
        all_results.extend(results);atomic('evolution_progress.json',dict(results=all_results))
        if not all(x['passed'] for x in results):break
    report=dict(completed=len(all_results)==6,passed=len(all_results)==6 and all(x['passed'] for x in all_results),results=all_results)
    atomic('evolution_report.json',report)
    raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run_basic_sequence()
