"""Approved full disk h=.005/.0025 comparison, conditional coarse continuation.

No steering, fixed masses, Re=infinity, true circular wall, constrained Psi.
Stops and preserves any numerical/physical gate failure. No automatic retry,
timestep reduction, reset, pressure repair, or replacement of old trajectories.
"""
import argparse,gc,hashlib,json,os,shutil,subprocess,sys,time,traceback,fcntl
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
from types import SimpleNamespace
import numpy as np

base=Path('/gpfs/kjhan/LagForce');fast=base/'Codex_full_disk_fast_modules_v1_20261019'
old=base/'Codex_full_disk_arc_roots_modules_v7_20261019'
sys.path.insert(0,str(fast))
import full_disk_stage_memory as stage_api
import full_disk_density_memory as density_api
import full_disk_velocity_memory as velocity_api
from full_disk_factor_memory import update
from moment_compiled_backend import Backend

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
p.add_argument('--mode',choices=('bundle','coarse','fine','continue'),default='bundle')
p.add_argument('--restart-root',type=Path)
p.add_argument('--stop-file',type=Path)
p.add_argument('--gpu-lock',type=Path)
args=p.parse_args();run=args.output.resolve();run.mkdir(exist_ok=False)
for name in ('work','steps','checkpoints','factors'): (run/name).mkdir()
started=time.perf_counter();providers=1;steps_accepted=0;current_label='preflight'
IC=base/'Codex_full_curved_kepler_IC_20261019/initial.npz'
initial={k:v for k,v in np.load(IC).items()};n=len(initial['ids'])
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
state_hash=lambda state:hashlib.sha256(np.ascontiguousarray(state['state']).view(np.uint8)).hexdigest()
def write(path,data):Path(path).write_text(json.dumps(data,indent=2))
def emit(data):
    row=dict(elapsed_seconds=time.perf_counter()-started,**data)
    with (run/'events.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
    print(json.dumps(row),flush=True)
def command(script,arguments,log):
    with Path(log).open('w') as out:
        subprocess.run([sys.executable,str(script),*map(str,arguments)],check=True,stdout=out,stderr=subprocess.STDOUT)
def load_npz(path):
    with np.load(path) as data:return {k:data[k] for k in data.files}
class SerializedGPU:
    """Independent CPU paths share one Ada without overlapping GPU buffers."""
    def __init__(self,backend,lock):self.backend=backend;self.lock=lock
    @property
    def path(self):return self.backend.path
    @property
    def kind(self):return self.backend.kind
    def invoke(self,method,*a):
        with self.lock.open('a') as handle:
            fcntl.flock(handle.fileno(),fcntl.LOCK_EX)
            return getattr(self.backend,method)(*a)
    def factor(self,*a):return self.invoke('factor',*a)
    def stage(self,*a):return self.invoke('stage',*a)
def save_state(path,state):
    path=Path(path);assert not path.exists(),str(path)
    np.savez(path,**state)
def latest(branch,state,step):
    temp=run/f'latest_{branch}_next.npz';save_state(temp,state)
    target=run/f'latest_{branch}.npz'
    if target.exists():os.replace(target,run/f'previous_{branch}.npz')
    os.replace(temp,target)
    write(run/f'latest_{branch}.json',dict(step=step,time=float(state['time']),checkpoint=str(target),sha256=sha(target),
        providers_started=providers,accepted_steps_total=steps_accepted,scope='Restart state; not a claim of checkpoint physical validation'))
def factor_update(base_state,linear,saved,h,time_value,label,shadow=False):
    try:return update(base_state,linear,saved,h,time_value,backend,shadow_audit=shadow)
    except Exception:
        failure=run/('failure_factor_'+label);failure.mkdir(exist_ok=False)
        save_state(failure/'base.npz',base_state);save_state(failure/'linear.npz',linear)
        witness=dict(saved);q=witness['incoming_Q']
        complete=q.nbytes<=1024**3
        if not complete:witness['incoming_Q']=q[:int(witness['incoming_offsets'][64])]
        np.savez(failure/'saved_stage.npz',**witness)
        write(failure/'manifest.json',dict(h=h,time=time_value,incoming_Q_complete=complete,
            incoming_Q_bytes=q.nbytes,limits='Complete raw rates and states; incoming columns capped at 1GiB, otherwise first64 cells'))
        raise

raw0=np.column_stack([initial['mass'],initial['state'][:-4].reshape(n,20)[:,3:]]).astype(np.longdouble)
r0=np.linalg.norm(initial['state'][:-4].reshape(n,20)[:,:2],axis=1)
physics=json.loads((IC.parent/'report.json').read_text())['physics'];P0=physics['P0'];gamma=physics['gamma']
def cheap_checks(state,mesh,step):
    assert np.array_equal(state['ids'],initial['ids']) and np.array_equal(state['mass'],initial['mass'])
    a=state['state'][:-4].reshape(n,20);raw=np.column_stack([state['mass'],a[:,3:]]).astype(np.longdouble)
    thermal=raw[:,17]-.5*(raw[:,13]+raw[:,16]);area=mesh['moments'][:,0]
    assert np.isfinite(raw).all() and np.isfinite(state['Psi']).all() and (thermal>0).all() and (area>0).all()
    pe=(gamma-1)*thermal/(P0*area)-1;mass=state['mass'].astype(np.longdouble)
    pressure={}
    for label,mask in dict(all=np.ones(n,bool),inner=r0<.4).items():
        pressure[label]=dict(RMS=float(np.sqrt(np.average(pe[mask]**2,weights=mass[mask]))),maximum=float(abs(pe[mask]).max()))
    change=raw-raw0
    vector=np.r_[change[:,7:9].sum(0),change[:,17].sum(),(change[:,10]-change[:,11]).sum()]
    ledger=state['state'][-4:]-initial['state'][-4:];residual=vector-ledger
    eps=np.finfo(float).eps;ops=128+4*n+4*step;g=ops*eps/(1-ops*eps)
    scale=np.r_[abs(raw[:,7:9]).sum(0)+abs(raw0[:,7:9]).sum(0),abs(raw[:,17]).sum()+abs(raw0[:,17]).sum(),
        abs(raw[:,10:12]).sum()+abs(raw0[:,10:12]).sum()]+abs(ledger)
    bound=g*scale;assert np.all(abs(residual)<=bound)
    combined=change[:,17].sum()+(state['Psi'].astype(np.longdouble)-initial['Psi'].astype(np.longdouble)).sum()
    energy_scale=abs(raw[:,17]).sum()+abs(raw0[:,17]).sum()+abs(state['Psi']).sum()+abs(initial['Psi']).sum()
    energy_bound=128*eps*(2*step+1)*energy_scale
    assert abs(combined)<=energy_bound
    return dict(pressure=pressure,minimum_thermal=float(thermal.min()),
        P_E_L_ledger_residual=np.asarray(residual,float).tolist(),P_E_L_ledger_bound=np.asarray(bound,float).tolist(),
        combined_energy_change=float(combined),combined_energy_roundoff_bound=float(energy_bound),mass_bitwise_initial=True)

def prepare_state(state,label,checkpoint=False):
    global current_label
    current_label=label;folder=run/('checkpoints' if checkpoint else 'work')/label;folder.mkdir()
    path=folder/'accepted.npz';save_state(path,state)
    geom=folder/'geometry'
    emit(dict(status='geometry_started',label=label,time=float(state['time'])))
    command(old/'prepare_full_curved_stage_geometry_workers.py',
        ['--state',path,'--workers',8,'--output',geom],folder/'geometry.log')
    native=json.loads((geom/'full_exact_supports.json').read_text());mesh=load_npz(geom/'mesh.npz')
    key=sha(path)+':'+sha(geom/'full_exact_supports.json')
    emit(dict(status='density_started',label=label,time=float(state['time']),independent_q48=checkpoint))
    fit=density_api.prepare(state,native,mesh,key,workers=8,independent_q48=checkpoint)
    write(folder/'density_prepare.json',fit.report)
    checks=cheap_checks(state,mesh,max(1,int(round(float(state['time'])/.0025))))
    write(folder/'state_checks.json',checks)
    result=SimpleNamespace(state=state,path=path,root=folder,geom=geom,native=native,mesh=mesh,key=key,
        fits=fit.rows,checkpoint=checkpoint,checks=checks)
    if checkpoint:checkpoint_audit(result)
    return result

def checkpoint_audit(bundle):
    folder=bundle.root;density=folder/'density';density.mkdir()
    with (density/'cells.jsonl').open('w') as f:
        for row in bundle.fits:f.write(json.dumps(row)+'\n')
    # Existing independent physical comparison consumes this measured q48
    # report; no second Newton fit or density reconstruction is performed.
    rows=bundle.fits
    write(density/'report.json',dict(status='passed_full_stage_Psi_density_preparation',state_sha256=sha(bundle.path),
        signed_q48_Psi_error=sum(row['Psi_error'] for row in rows),absolute_q48_Psi_error=sum(abs(row['Psi_error']) for row in rows),
        maximum_q48_Psi_error=max(abs(row['Psi_error']) for row in rows),cells=n,
        time=float(bundle.state['time']),source='Full same-state memory fits with independent q48'))
    emit(dict(status='checkpoint_reference_started',label=folder.name,time=float(bundle.state['time'])))
    command(old/'project_steady_reference_on_full_curved_state_workers.py',
        ['--state',bundle.path,'--geometry',bundle.geom,'--workers',8,'--output',folder/'reference'],folder/'reference.log')
    command(fast/'compare_full_curved_checkpoint_physics.py',
        ['--state',bundle.path,'--geometry',bundle.geom,'--density',density,'--reference',folder/'reference','--output',folder/'physics'],folder/'physics.log')
    report=json.loads((folder/'physics/report.json').read_text());gates={}
    for label,row in report['regions'].items():
        if not row['cells']:continue
        pressure=row['relative_pressure_error'];velocity=row['fields']['velocity']['relative_RMS']
        phase=row.get('velocity_direction_error_degrees')
        gates[label]=dict(pressure_RMS=pressure['RMS'],pressure_maximum=pressure['maximum'],velocity_relative_RMS=velocity,
            phase_RMS_degrees=phase['RMS'] if phase else None,
            passed=pressure['RMS']<.01 and pressure['maximum']<.1 and velocity is not None and velocity<.01 and (not phase or phase['RMS']<1.))
    passed=all(row['passed'] for row in gates.values())
    write(folder/'physical_gate.json',dict(status='passed_checkpoint_physical_monitor' if passed else 'failed_checkpoint_physical_monitor',
        time=float(bundle.state['time']),regions=gates,criteria_source='Astra 2026-10-19 operational monitoring; not final long-time success tolerance'))
    emit(dict(status='checkpoint_complete',label=folder.name,time=float(bundle.state['time']),physical_gate_passed=passed,pressure={k:v['pressure_RMS'] for k,v in gates.items()}))
    assert passed,'Checkpoint physical monitoring failed; state and diagnostics preserved'

def release(bundle):
    if bundle is not None and not bundle.checkpoint:
        assert bundle.root.parent==run/'work'
        # Only disposable work created by this invocation; old outputs and
        # restart checkpoints never enter this directory.
        shutil.rmtree(bundle.root)

def provider(bundle,label):
    global providers,current_label
    assert providers<439;current_label=label
    emit(dict(status='velocity_gravity_started',label=label,time=float(bundle.state['time'])))
    velocity=velocity_api.prepare(bundle.state,bundle.native,bundle.mesh,bundle.fits,bundle.key,workers=8)
    providers+=1
    emit(dict(status='actual_provider_started',label=label,time=float(bundle.state['time']),unique_provider_index=providers))
    result=stage_api.evaluate(bundle.state,bundle.native,bundle.mesh,bundle.fits,velocity.prepared,velocity.gravity_nodes,
        bundle.key,backend,run/'steps'/label,callback=lambda row:emit(dict(label=label,provider_event=row)),provider_budget=439)
    write(run/'steps'/label/'preparation.json',dict(velocity_gravity=velocity.report,geometry=json.loads((bundle.geom/'report.json').read_text()),
        input_state_sha256=sha(bundle.path),geometry_support_sha256=sha(bundle.geom/'full_exact_supports.json')))
    del velocity.gravity_nodes;velocity.owner.close();del velocity;gc.collect()
    emit(dict(status='actual_provider_complete',label=label,time=float(bundle.state['time']),unique_provider_index=providers,seconds=result.report['total_seconds']))
    return result.saved

checkpoint_steps={'coarse':{2,5,10,20,50,100,200},'fine':{4,10,20}}
def evolve(branch,h,start_state,start_step,end_step,bundle=None):
    global steps_accepted,current_label
    state=start_state
    for step in range(start_step+1,end_step+1):
        if args.stop_file is not None and args.stop_file.exists():
            raise RuntimeError('Sibling path failed; stopped at the next step boundary with restart state retained')
        label=f'{branch}_{step:04d}';current_label=label
        if step==1:
            assert float(state['time'])==0.;first=initial_stage
        else:
            assert bundle is not None
            first=provider(bundle,label+'_start')
        midpoint,half_report=factor_update(state,state,first,h/2,(step-.5)*h,label+'_midpoint')
        write(run/'factors'/(label+'_midpoint.json'),half_report)
        if step!=1:del first
        release(bundle);bundle=None
        half_bundle=prepare_state(midpoint,label+'_midpoint',False)
        half_stage=provider(half_bundle,label+'_midpoint')
        candidate,full_report=factor_update(state,midpoint,half_stage,h,step*h,label+'_full',shadow=step in checkpoint_steps[branch])
        write(run/'factors'/(label+'_full.json'),full_report)
        del half_stage,midpoint;release(half_bundle);del half_bundle;gc.collect()
        # Endpoint geometry/fit admissibility is required before numerical
        # acceptance; the same preparation is reused by the next stage0.
        bundle=prepare_state(candidate,label,step in checkpoint_steps[branch])
        state=candidate;steps_accepted+=1;latest(branch,state,step)
        emit(dict(status='step_accepted',branch=branch,step=step,time=float(state['time']),h=h,
            total_accepted_steps=steps_accepted,providers_started=providers,checks=bundle.checks))
        write(run/'progress.json',dict(status='running',branch=branch,step=step,time=float(state['time']),
            accepted_steps=steps_accepted,providers_started=providers,latest_label=label))
    return state,bundle

def compare_paths():
    results=[]
    for cs,fs in ((2,4),(5,10),(10,20)):
        ca=run/'coarse/checkpoints'/f'coarse_{cs:04d}'/'physics';fb=run/'fine/checkpoints'/f'fine_{fs:04d}'/'physics'
        a=load_npz(ca/'fields.npz');b=load_npz(fb/'fields.npz');assert np.array_equal(a['ids'],b['ids'])
        assert np.array_equal(a['actual_mass'],b['actual_mass']);m=a['actual_mass'].astype(np.longdouble)
        regions={};passed=True
        for label in ('all','inner','cavity','disk_core','inner_contact','outer_contact','dynamic_buffer','wall'):
            mask=a[label]|b[label]
            if not mask.any():continue
            w=m[mask];rms=lambda x:float(np.sqrt(np.sum(w*np.asarray(x,dtype=np.longdouble)**2)/w.sum()))
            pe=b['relative_pressure_error'][mask]-a['relative_pressure_error'][mask]
            ve=(b['actual_velocity']-b['reference_velocity'])-(a['actual_velocity']-a['reference_velocity'])
            vr=np.sqrt(np.sum(a['reference_velocity'][mask]**2,axis=1))
            velocity=rms(np.linalg.norm(ve[mask],axis=1))/rms(vr)
            valid=mask&a['phase_defined']&b['phase_defined']
            phase=None
            if valid.any():
                d=b['velocity_direction_error_degrees'][valid]-a['velocity_direction_error_degrees'][valid]
                d=(d+180)%360-180;phase=float(np.sqrt(np.average(d*d,weights=m[valid])))
            ok=rms(pe)<.002 and velocity<.002 and (phase is None or phase<.2);passed &= ok
            regions[label]=dict(pressure_error_difference_RMS=rms(pe),velocity_error_difference_relative_RMS=velocity,
                phase_error_difference_RMS_degrees=phase,passed=ok)
        results.append(dict(time=cs*.005,regions=regions,passed=passed))
    write(run/'two_dt_physical_gate.json',dict(status='passed_two_dt_operational_gate' if all(row['passed'] for row in results) else 'failed_two_dt_operational_gate',checkpoints=results))
    return all(row['passed'] for row in results)

try:
    assert os.environ['OPENBLAS_NUM_THREADS']==os.environ['OMP_NUM_THREADS']=='1' and os.environ['OMP_PROC_BIND']=='false'
    requirements={'factor_memory_check':'passed_factor_memory_against_accepted','fast_provider_comparison':'passed_provider_replay_bitwise',
        'density_memory_check':'passed_production_density_fits_bitwise','velocity_memory_check':'passed_in_memory_velocity_bitwise'}
    for label,status in requirements.items():assert json.loads((base/f'Codex_full_disk_{label}_20261019/report.json').read_text())['status']==status
    free=os.statvfs(run).f_bavail*os.statvfs(run).f_frsize;estimated_output=8*1024**3
    assert free>2*estimated_output
    factor_seconds=json.loads((base/'Codex_full_disk_factor_memory_check_20261019/report.json').read_text())['factor']['seconds']
    # Components measured on the identical moving state. Provider preparation
    # previously included file load/hash: excluded here, with a separate 0–42s
    # allowance for unmeasured in-memory state restoration and changing cells.
    stage_seconds=21.957+66.1046+40.554+(227.8546-41.648)
    forecasts={}
    for target,stages,steps in [('conditional_t1',438,220),('same_h_t6_projection',2438,1220),('same_h_t10_projection',4038,2020)]:
        forecasts[target]=dict(core_hours=(stages*stage_seconds+2*steps*factor_seconds)/3600,
            with_restoration_allowance_hours=(stages*(stage_seconds+42)+2*steps*factor_seconds)/3600)
    source_files=[fast/name for name in ('full_disk_stage_memory.py','full_disk_factor_memory.py','full_disk_density_memory.py','full_disk_velocity_memory.py','curved_face_nodes_fast.py','curved_contact_priors_fast.py','same_stage_face_cache.py','closed_root_transport_arrays.py','face_record_batch_fast.py','native_exact_prior_bisection.cpp','native_exact_prior_bisection.py','libexact_prior.so','compare_full_curved_checkpoint_physics.py')]+[Path(__file__)]
    manifest=dict(status='starting_approved_continuous_bundle',mode=args.mode,initial_checkpoint=str(IC),initial_sha256=sha(IC),
        cells=n,Re='infinity',physical_viscosity=0.,steering=False,fixed_mass=True,geometry='power cells intersected with true radius2.5 disk',
        steps=dict(coarse_h=.005,fine_h=.0025,comparison_end=.05,conditional_end=1.,maximum_accepted_steps=220,maximum_unique_providers=439),
        checkpoint_times=dict(coarse=[.01,.025,.05,.1,.25,.5,1.],fine=[.01,.025,.05]),
        output_preflight=dict(initial_checkpoint_bytes=IC.stat().st_size,checkpoint_count=10,estimated_total_bytes=estimated_output,free_bytes=free,
            policy='Compact restart states/checkpoints and summaries; anonymous gravity buffers and disposable successful work; failures retained'),
        cost_preflight=dict(measured_component_stage_seconds=stage_seconds,measured_factor_update_seconds=factor_seconds,
            forecasts=forecasts,limits='Component timing forecast, not measured whole-production timing. Checkpoint audits extra; state/node count and outer iterations may change.'),
        source_sha256={str(q):sha(q) for q in source_files+[base/'Codex_compiled_midpoint_build_20261018/liblagmoment.so',old/'prepare_full_curved_stage_geometry_workers.py',old/'project_steady_reference_on_full_curved_state_workers.py']},python=sys.version,worker_environment={key:os.environ[key] for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','OMP_PROC_BIND')},pid=os.getpid(),host=os.uname().nodename)
    write(run/'manifest.json',manifest);emit(dict(status='preflight_passed',output_estimate_GiB=8,cost_forecasts=forecasts))
    if args.mode=='bundle':
        # Fixed independent IC/dt paths. The coordinator waits for process
        # completions, never polls Slurm or SSH in a loop.
        available=int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:')))*1024
        assert len(os.sched_getaffinity(0))>=20 and 2*(8*16*1024**3+24*1024**3)<available/2
        write(run/'concurrency.json',dict(independent_paths=2,CPU_workers_per_path=8,worker_BLAS_OMP_threads=1,
            GPU_policy='One shared Ada with GPU calls serialized by an advisory lock; CPU preparation runs concurrently',available_host_memory=available))
        stop_file=run/'STOP';gpu_lock=run/'GPU.lock'
        def child(mode):
            commandline=[sys.executable,str(Path(__file__).resolve()),'--mode',mode,'--output',str(run/mode),
                '--stop-file',str(stop_file),'--gpu-lock',str(gpu_lock)]
            if mode=='continue':commandline+=['--restart-root',str(run/'coarse')]
            with (run/(mode+'.log')).open('w') as out:
                status=subprocess.run(commandline,stdout=out,stderr=subprocess.STDOUT).returncode
            return mode,status
        failures=[]
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures=[pool.submit(child,mode) for mode in ('coarse','fine')]
            for future in as_completed(futures):
                mode,status=future.result();emit(dict(status='comparison_child_terminal',branch=mode,exit_code=status))
                if status:failures.append(mode);stop_file.touch(exist_ok=True)
        if failures:raise RuntimeError(dict(status='independent_comparison_path_failed',branches=failures))
        if not compare_paths():
            write(run/'result.json',dict(status='stopped_two_dt_operational_gate',accepted_steps=30,providers_started=59,long_time_stability_verified=False))
        else:
            emit(dict(status='conditional_coarse_continuation_started',from_time=.05,to_time=1.,h=.005))
            mode,status=child('continue')
            if status:raise RuntimeError('Conditional continuation failed; child restart states retained')
            write(run/'result.json',dict(status='completed_conditional_t1_bundle',accepted_steps=220,providers_started=439,
                coarse_time=1.,fine_time=.05,long_time_stability_verified=False,limits='t1 intermediate validation; historical t6–10 target remains.'))
    elif args.mode in ('coarse','fine'):
        backend=SerializedGPU(Backend(base/'Codex_compiled_midpoint_build_20261018/liblagmoment.so','cuda',1),args.gpu_lock or run/'GPU.lock')
        initial_stage=load_npz(base/'Codex_full_disk_initial_stage0_arc_roots_retry_20261019/stage0.npz')
        assert str(initial_stage['initial_state_sha256'])==state_hash(initial)
        initial.setdefault('weight_rate_seed',np.zeros(n))
        h,end=(.005,10) if args.mode=='coarse' else (.0025,20)
        emit(dict(status='comparison_evolution_started',branch=args.mode,h=h,initial_RHS_reused=True))
        state,bundle=evolve(args.mode,h,initial,0,end)
        write(run/'result.json',dict(status='completed_comparison_branch',branch=args.mode,time=.05,
            accepted_steps=steps_accepted,providers_started=providers,long_time_stability_verified=False))
    else:
        assert args.restart_root is not None
        backend=SerializedGPU(Backend(base/'Codex_compiled_midpoint_build_20261018/liblagmoment.so','cuda',1),args.gpu_lock or run/'GPU.lock')
        providers=59;steps_accepted=30
        folder=args.restart_root/'checkpoints/coarse_0010';state=load_npz(folder/'accepted.npz')
        assert float(state['time'])==.05
        assert json.loads((folder/'physical_gate.json').read_text())['status']=='passed_checkpoint_physical_monitor'
        geom=folder/'geometry';mesh=load_npz(geom/'mesh.npz');native=json.loads((geom/'full_exact_supports.json').read_text())
        fits=[json.loads(line) for line in (folder/'density/cells.jsonl').read_text().splitlines()]
        bundle=SimpleNamespace(state=state,path=folder/'accepted.npz',root=folder,geom=geom,native=native,mesh=mesh,
            key=sha(folder/'accepted.npz')+':'+sha(geom/'full_exact_supports.json'),fits=fits,checkpoint=True)
        coarse,bundle=evolve('coarse',.005,state,10,200,bundle)
        write(run/'result.json',dict(status='completed_conditional_t1_bundle',accepted_steps=steps_accepted,providers_started=providers,
            coarse_time=1.,fine_time=.05,long_time_stability_verified=False,limits='t1 intermediate validation; historical t6–10 target remains.'))
except Exception as exc:
    children={}
    if args.mode=='bundle':
        for name in ('coarse','fine','continue'):
            for filename in ('result.json','progress.json'):
                path=run/name/filename
                if path.exists():children[name]=json.loads(path.read_text());break
    write(run/'result.json',dict(status='failed_continuous_bundle',current_label=current_label,accepted_steps=steps_accepted,
        providers_started=providers,error=repr(exc),traceback=traceback.format_exc(),long_time_stability_verified=False,
        retry=False,timestep_reduced=False,old_outputs_preserved=True,child_states=children,
        counts_scope='Coordinator local counts only; child_states authoritative for physical progress' if args.mode=='bundle' else 'This physical branch (continuation includes preceding30 steps/59 providers)'))
    emit(dict(status='bundle_failed',label=current_label,error=repr(exc)))
    raise
