"""Recover a stopped owned run and change only its geometry integrator.

Approved dynamical arrays and accumulated audits remain unchanged. Any saved
unapproved RK2 fluid stages are reused; their legacy geometry sample cache is
retained separately because it describes a different path.
"""
import json,hashlib,shutil,os
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic,restore
from saved_curved_contact import checked_json
from simple_wave_evolution import finish_pending_commit,install_domain
from geometry_mass_lift import integrate,lifted_state


def run():
    config=json.loads(Path('geometry_handoff_config.json').read_text());source=Path(config['source'])
    if Path('geometry_handoff_ready.json').exists():
        print('saved handoff already prepared; preserve all new approvals',flush=True);return
    validation=json.loads(Path(config['time_validation']).read_text())
    orders=[r['observed_order'] for r in validation['fields'].values()]
    if not validation['all_numerical_gates_passed'] or any(p is None or not 1.8<p<2.2 for p in orders):
        raise ValueError('geometry integrator time validation failed')
    if 'LagEunha' not in (source/'OWNER').read_text():raise ValueError('source project mismatch')
    if not (source/'execution.txt').read_text().splitlines()[-1].startswith('exit='):raise ValueError('source still running')
    checked_json(source/'checkpoint.json',config['checkpoint_sha256'])
    for p in source.glob('*.json'):
        if p.name not in ('geometry_handoff_config.json','geometry_handoff_ready.json'):shutil.copy2(p,p.name)
    Path('preserved').mkdir(exist_ok=True)
    for p in (source/'preserved').glob('*.json'):shutil.copy2(p,Path('preserved')/p.name)
    old_config=Path('simple_wave_config.json').read_bytes();old_identity=hashlib.sha256(old_config).hexdigest();cfg=json.loads(old_config)
    finish_pending_commit(Path('pending.json'),Path('checkpoint.json'),Path('pending_samples.json'),old_identity)
    before=json.loads(Path('checkpoint.json').read_text());original=json.loads(Path('initial_state.json').read_text())
    if before['config_sha256']!=old_identity or original['config_sha256']!=old_identity:raise ValueError('source configuration identity mismatch')
    np.testing.assert_array_equal(before['mass'],original['mass'])
    for name in ['simple_wave_config.json','checkpoint.json','initial_state.json','pending.json','pending_samples.json']:
        p=Path(name)
        if p.exists():shutil.copy2(p,Path('preserved')/('before_geometry_handoff_'+name))
    cfg.update(geometry_path='material_mass_lift',max_dt=config['max_dt'],geometry_step_controller='physical_path')
    for key in ('path_sample_lattice','path_save_batch','geometry_rtol','geometry_atol','dense_moment_arithmetic','geometry_mass_root_tolerance'):
        if key in config:cfg[key]=config[key]
    cfg['geometry_handoff_source']=str(source)
    cfg['geometry_handoff_note']='All accepted states/raw moments/factors/masses/budget accumulators preserved. Saved unapproved midpoint fluid stages reused. Independent validated geometry integration replaces only the old interpolated site/weight path.'
    Path('simple_wave_config.json').write_text(json.dumps(cfg,indent=2)+'\n');identity=hashlib.sha256(Path('simple_wave_config.json').read_bytes()).hexdigest()
    for name in ['initial_state.json','checkpoint.json']:
        data=json.loads(Path(name).read_text());data['config_sha256']=identity
        if name=='checkpoint.json' and 'next_dt_override' in config:data['next_dt']=config['next_dt_override']
        atomic(name,data)
    if 'next_dt_override' in config:
        if Path('pending.json').exists():Path('pending.json').replace('preserved_pending_before_dt_probe.json')
        if Path('pending_samples.json').exists():Path('pending_samples.json').replace('preserved_pending_samples_before_dt_probe.json')
    install_domain();pending_reused=False
    if Path('pending.json').exists():
        trial=json.loads(Path('pending.json').read_text());record=restore(trial['record'])
        if cfg.get('dense_moment_arithmetic')=='longdouble':record['dense_moment_arithmetic']='longdouble'
        if 'geometry_lift' not in record or config.get('refine_geometry',False):
            record.pop('geometry_lift',None)
            record['geometry_lift']=integrate(record,np.asarray(before['mass']),rtol=cfg.get('geometry_rtol',1e-12),atol=cfg.get('geometry_atol',1e-14),root_tolerance=cfg.get('geometry_mass_root_tolerance',2e-13))
            endpoint,_,_=lifted_state(record,np.asarray(before['mass']),1.)
            candidate=np.asarray(trial['new']);saved_moments=candidate[:-4].reshape(-1,20)[:,3:].copy()
            candidate[:-4].reshape(-1,20)[:,:3]=endpoint[:-4].reshape(-1,20)[:,:3]
            np.testing.assert_array_equal(candidate[:-4].reshape(-1,20)[:,3:],saved_moments)
            trial.update(record=record,new=candidate);trial.pop('geometry_pending',None);atomic('pending.json',trial)
            if Path('pending_samples.json').exists():
                archive='preserved_legacy_pending_samples_'+hashlib.sha256(Path('pending_samples.json').read_bytes()).hexdigest()+'.json'
                Path('pending_samples.json').replace(archive)
        pending_reused=True
    after=json.loads(Path('checkpoint.json').read_text());after['config_sha256']=old_identity
    if 'next_dt_override' in config:after['next_dt']=before['next_dt']
    if after!=before:raise RuntimeError('approved checkpoint changed during handoff')
    report=dict(prepared=True,source=str(source),source_time=before['t'],accepted_steps=before['steps'],rejected=before['rejected'],
        checkpoint_sha256_before=config['checkpoint_sha256'],fixed_mass_and_all_accepted_arrays_unchanged=True,pending_fluid_stages_reused=pending_reused,
        config_sha256=identity,time_orders=orders,geometry_path='material_mass_lift',next_dt_before=before['next_dt'],next_dt_override=config.get('next_dt_override'))
    atomic('geometry_handoff_ready.json',report);print(report,flush=True)


if __name__=='__main__':run()
