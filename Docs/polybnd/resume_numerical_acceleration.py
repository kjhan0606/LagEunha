"""Preserve an approved trajectory and its pending work across equivalent kernels.

Only a proven computational backend changes. Physical model, numerical budgets,
accepted arrays, counters, accumulated defects and pending samples are retained.
Historical archives remain in the immutable source directory and are referenced.
"""
import fcntl,hashlib,json,os,shutil
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from simple_wave_evolution import finish_pending_commit

RUNTIME=('initial_state.json','initial_score.json','checkpoint.json','source_rhs.json','pending.json','pending_samples.json')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reidentify(data,identity):
    result=dict(data);result['config_sha256']=identity
    return result


def validate_arrays(initial,cp,cfg):
    mass=np.asarray(initial['mass']);state=np.asarray(cp['state']);n=cfg['n']**2
    if mass.shape!=(n,) or state.shape!=(20*n+4,) or len(cp['factors'])!=n:raise ValueError('acceleration source shape mismatch')
    if not np.array_equal(cp['mass'],mass) or np.any(mass<=0):raise ValueError('acceleration source fixed mass mismatch')
    if not np.all(np.isfinite(state)) or not initial['t']<=cp['t']<cfg['end']:raise ValueError('invalid acceleration source time/state')


def prepare(root=Path('.')):
    root=Path(root);ready=root/'numerical_resume_ready.json';config_path=root/'simple_wave_config.json'
    if ready.exists():
        marker=json.loads(ready.read_text());cp=json.loads((root/'checkpoint.json').read_text())
        if marker['config_sha256']!=sha(config_path) or cp['config_sha256']!=marker['config_sha256']:raise ValueError('prepared acceleration identity mismatch')
        return dict(marker,already_prepared=True)
    control=json.loads((root/'numerical_resume_config.json').read_text());source=Path(control['source']).resolve()
    if source==root.resolve() or source.stat().st_uid!=os.getuid() or 'LagEunha' not in (source/'OWNER').read_text():raise ValueError('unowned or identical acceleration source')
    if not any(v.startswith('exit=') for v in (source/'execution.txt').read_text().splitlines()):raise ValueError('acceleration source has not terminated')
    with (source/'.run.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise ValueError('acceleration source is still locked/running')
        for name,key in [('simple_wave_config.json','source_config_sha256'),('initial_state.json','source_initial_sha256'),('checkpoint.json','source_checkpoint_sha256')]:
            if sha(source/name)!=control[key]:raise ValueError('acceleration source hash mismatch: '+name)
        cfg=json.loads((source/'simple_wave_config.json').read_text());old_id=control['source_config_sha256']
        original=root/'resume_original';original.mkdir(exist_ok=True)
        shutil.copy2(source/'simple_wave_config.json',original/'source_config.json')
        copied={}
        for name in RUNTIME:
            path=source/name
            if path.exists():
                shutil.copy2(path,original/name);shutil.copy2(path,root/name);copied[name]=sha(path)
        initial=json.loads((root/'initial_state.json').read_text());cp=json.loads((root/'checkpoint.json').read_text())
        if initial['config_sha256']!=old_id or cp['config_sha256']!=old_id:raise ValueError('source configuration identity mismatch')
        validate_arrays(initial,cp,cfg)
        pending=root/'pending.json';samples=root/'pending_samples.json'
        trial=json.loads(pending.read_text()) if pending.exists() else None
        committed=bool(trial is not None and 'commit' in trial)
        if trial is not None and trial.get('sample_lattice_upgrade'):
            name=trial['sample_lattice_upgrade']['archive']
            if Path(name).name!=name:raise ValueError('invalid pending lattice archive')
            shutil.copy2(source/name,root/name)
        finish_pending_commit(pending,root/'checkpoint.json',samples,old_id)
        cp=json.loads((root/'checkpoint.json').read_text());validate_arrays(initial,cp,cfg)
        rhs_path=root/'source_rhs.json'
        if rhs_path.exists():
            rhs=json.loads(rhs_path.read_text());state_hash=hashlib.sha256(np.asarray(cp['state']).tobytes()).hexdigest()
            if rhs['state_sha256']!=state_hash or rhs['t']!=cp['t']:raise ValueError('cached source RHS does not match recovered source')
        if pending.exists():
            trial=json.loads(pending.read_text());record=trial['record']
            if record['t']!=cp['t'] or record['dt']!=cp['next_dt'] or not np.array_equal(record['y'],cp['state']):raise ValueError('pending work does not extend the recovered source')
        for key in ('preserved_initial_state','preserved_source_rhs','preserved_source_config'):
            if key not in cfg:continue
            name=Path(cfg[key])
            if name.is_absolute():continue
            if '..' in name.parts:raise ValueError('relative preserved input escapes source')
            target=root/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source/name,target)
            if sha(target)!=cfg[key+'_sha256']:raise ValueError('preserved acceleration input hash mismatch')
        proof=Path(control['equivalence_report'])
        if sha(proof)!=control['equivalence_report_sha256']:raise ValueError('acceleration equivalence report mismatch')
        verification=json.loads(proof.read_text())
        if not verification['completed'] or verification['rhs_maximum_absolute_difference']>=1e-12 or verification['independent_mass_flux_maximum_absolute_difference']>=1e-12:raise ValueError('unverified acceleration kernel equivalence')
        if verification.get('physical_covariance_rate_maximum_absolute_difference',0.)>=1e-12:raise ValueError('unverified physical covariance rate equivalence')
        cfg.update(quadratic_budget_backend='batched',numerical_acceleration_source=str(source),numerical_acceleration_equivalence_report=str(proof),
            resolve_mass_quadrature_before_reject=cfg.get('resolve_mass_quadrature_before_reject',False))
        if control.get('pressure_evaluation_backend')=='horner':
            if cfg.get('pressure_reconstruction') not in ('polynomial_cell_mean_degree3','polynomial_cell_mean_degree4'):raise ValueError('Horner acceleration requires a higher polynomial pressure')
            cfg['pressure_evaluation_backend']='horner'
        for key,allowed in [('mass_path_evaluation_backend','cell_batched'),('moment_flux_evaluation_backend','batched'),('power_cell_evaluation_backend','seeded_active')]:
            if key in control:
                if control[key]!=allowed:raise ValueError('unsupported acceleration backend: '+key)
                cfg[key]=allowed
        cfg['note']=cfg.get('note','')+' Equivalent numerical-kernel continuation: accepted primary arrays, budgets, counters and pending fluid/path work retained; original history stays in the source directory.'
        atomic(config_path,cfg);identity=sha(config_path)
        atomic(root/'initial_state.json',reidentify(initial,identity));atomic(root/'checkpoint.json',reidentify(cp,identity))
        saved=json.loads((root/'checkpoint.json').read_text())
        for key in cp:
            if key!='config_sha256' and saved[key]!=cp[key]:raise ValueError('approved source changed during acceleration handoff: '+key)
        if pending.exists() and sha(pending)!=copied['pending.json']:raise ValueError('uncommitted fluid stages changed during acceleration handoff')
        if samples.exists() and sha(samples)!=copied['pending_samples.json']:raise ValueError('cached path samples changed during acceleration handoff')
        marker=dict(prepared=True,source=str(source),source_time=cp['t'],steps=cp['steps'],rejected=cp['rejected'],all_accepted_arrays_and_accumulators_unchanged=True,
            source_pending_commit_completed=committed,pending_fluid_stages_reused=pending.exists(),pending_sample_cache_reused=samples.exists(),source_files_sha256=copied,
            config_sha256=identity,physical_model_and_budgets_unchanged=True,numerical_backend_changes={k:cfg[k] for k in ['quadratic_budget_backend','pressure_evaluation_backend','mass_path_evaluation_backend','moment_flux_evaluation_backend','power_cell_evaluation_backend'] if k in cfg},note='Equivalent computational kernels only: same geometry constraints, pressure polynomial, positive samples and shared moment products. Native mass-quadrature rejection policy remains unchanged by default. Original archives are referenced, not recopied.')
        atomic(ready,marker);return marker


if __name__=='__main__':print(prepare(),flush=True)
