"""Preserve an approved prefix and compare a declared numerical error allocation.

The physical RHS/velocity closure/smoothing are unchanged. The geometry update
returns to the same material-stage RK2 polynomial. Saved unapproved fluid stages
can be reused exactly; the older geometry and audit samples are retained.
"""
import json,hashlib,shutil
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from simple_wave_evolution import finish_pending_commit


def native_trial(trial):
    record=trial['record'];a=np.asarray(record['y'])[:-4].reshape(-1,20);dm=np.asarray(record['dm'])[:-4].reshape(-1,20)
    new=np.asarray(trial['new']);old_moments=new[:-4].reshape(-1,20)[:,3:].copy()
    new[:-4].reshape(-1,20)[:,:3]=a[:,:3]+record['dt']*dm[:,:3]
    np.testing.assert_array_equal(new[:-4].reshape(-1,20)[:,3:],old_moments)
    record.pop('geometry_lift',None);trial.pop('geometry_pending',None);trial.update(record=record,new=new)
    return trial


def run():
    control=json.loads(Path('native_resume_config.json').read_text());source=Path(control['source'])
    if Path('native_resume_ready.json').exists():print('saved native continuation already prepared',flush=True);return
    if 'LagEunha' not in (source/'OWNER').read_text():raise ValueError('wrong source project')
    if not (source/'execution.txt').read_text().splitlines()[-1].startswith('exit='):raise ValueError('source still running')
    if hashlib.sha256((source/'initial_state.json').read_bytes()).hexdigest()!=control['source_initial_sha256']:raise ValueError('source initial hash mismatch')
    checkpoint_source=source/'checkpoint.json'
    actual=hashlib.sha256(checkpoint_source.read_bytes()).hexdigest() if checkpoint_source.exists() else None
    if actual!=control.get('source_checkpoint_sha256'):raise ValueError('source checkpoint hash mismatch')
    for p in source.glob('*.json'):
        if p.name not in ('native_resume_config.json','native_resume_ready.json'):shutil.copy2(p,p.name)
    Path('preserved').mkdir(exist_ok=True)
    for p in (source/'preserved').glob('*.json'):shutil.copy2(p,Path('preserved')/p.name)
    old_data=Path('simple_wave_config.json').read_bytes();old_identity=hashlib.sha256(old_data).hexdigest();cfg=json.loads(old_data)
    for key,value in control['expected'].items():
        if cfg.get(key,'nonlinear_simple_wave' if key=='problem' else None)!=value:raise ValueError('source benchmark mismatch: '+key)
    finish_pending_commit(Path('pending.json'),Path('checkpoint.json'),Path('pending_samples.json'),old_identity)
    original=json.loads(Path('initial_state.json').read_text());mass=np.asarray(original['mass'])
    if original['config_sha256']!=old_identity:raise ValueError('source initial identity mismatch')
    if Path('checkpoint.json').exists():before=json.loads(Path('checkpoint.json').read_text())
    else:
        trial=json.loads(Path('pending.json').read_text()) if Path('pending.json').exists() else None
        if trial and (trial['record']['t']!=original['t'] or not np.array_equal(trial['record']['y'],original['state'])):raise ValueError('missing checkpoint after an approved step')
        zeros=np.zeros(len(mass));before=dict(config_sha256=old_identity,t=original['t'],next_dt=trial['record']['dt'] if trial else cfg['max_dt'],steps=0,rejected=0,
            mass=mass,state=original['state'],factors=original['factors'],peak=dict(rate=0.,generator=0.,ledger=0.,stage_mass=0.,min_pressure=1e300,min_factor_R=0.),
            absolute=zeros,signed=zeros,quad_delta=zeros,quad_error_bound=zeros,absolute_quad_error_bound=zeros)
        atomic('checkpoint.json',before);before=json.loads(Path('checkpoint.json').read_text())
    if before['config_sha256']!=old_identity:raise ValueError('source checkpoint identity mismatch')
    np.testing.assert_array_equal(before['mass'],mass)
    for name in ['simple_wave_config.json','initial_state.json','checkpoint.json','pending.json','pending_samples.json']:
        p=Path(name)
        if p.exists():shutil.copy2(p,Path('preserved')/('before_native_budget_'+name))
    old_budgets={k:cfg[k] for k in ['path_mass_budget','path_quadrature_budget','path_rate_budget','generator_path_budget']}
    cfg.update(control['budgets']);cfg.update(geometry_path='quadratic_midpoint',geometry_step_controller='combined',dense_moment_arithmetic='longdouble',path_sample_lattice=1024,path_save_batch=8)
    cfg['native_budget_source']=str(source);cfg['native_budget_transition_time']=before['t'];cfg['native_budget_original']=old_budgets
    cfg['note']='Same physical RHS and primary moments; material RK2 geometry with explicit numerical-error allocation study. All original accepted arrays retained; physical exact-field goals unchanged. Source strict-budget run is not reclassified.'
    for key in ['geometry_rtol','geometry_atol','geometry_mass_root_tolerance','geometry_precision_policy','geometry_precision_divisor']:cfg.pop(key,None)
    Path('simple_wave_config.json').write_text(json.dumps(cfg,indent=2)+'\n');identity=hashlib.sha256(Path('simple_wave_config.json').read_bytes()).hexdigest()
    for name in ['initial_state.json','checkpoint.json']:
        data=json.loads(Path(name).read_text());data['config_sha256']=identity;atomic(name,data)
    pending_reused=False
    if Path('pending.json').exists():
        trial=json.loads(Path('pending.json').read_text());trial=native_trial(trial);trial['record']['dense_moment_arithmetic']='longdouble';atomic('pending.json',trial);pending_reused=True
        if Path('pending_samples.json').exists():Path('pending_samples.json').replace('preserved_prior_geometry_samples.json')
    after=json.loads(Path('checkpoint.json').read_text());after['config_sha256']=old_identity
    if after!=before:raise RuntimeError('approved physical state changed in numerical-budget transition')
    result=dict(prepared=True,source=str(source),source_time=before['t'],steps=before['steps'],rejected=before['rejected'],all_accepted_arrays_unchanged=True,
        pending_fluid_stages_reused=pending_reused,original_budgets=old_budgets,new_budgets=control['budgets'],physical_goals_unchanged=True,config_sha256=identity)
    atomic('native_resume_ready.json',result);print(result,flush=True)


if __name__=='__main__':run()
