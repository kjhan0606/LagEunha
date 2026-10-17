"""Read-only endpoint sensitivity to error budgets or timestep caps."""
import hashlib,json,os
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from analyze_simple_wave_time import sampled_fields
from simple_wave_evolution import install_domain


def scales_and_goals(config):
    problem=config.get('problem','nonlinear_simple_wave')
    if problem=='nonlinear_simple_wave':
        scales=dict(density=config['amplitude'],velocity=config['amplitude'],pressure=config['amplitude'])
    elif problem=='steady_vortex':
        from steady_vortex_exact import parameters
        scales=dict(density=1.,velocity=config['amplitude'],pressure=parameters(config['amplitude'])[1])
    else:raise ValueError('unknown physical budget benchmark')
    goals={k:config.get('physical_goals',{}).get(k+'_field',config['physical_goal']) for k in scales}
    return scales,goals


def compatible(configs,initials):
    keys=('n','problem','amplitude','end','physical_goal','physical_goals','mass_root_backend','velocity_reconstruction','pressure_reconstruction',
          'kinetic_smoothing_coefficient','kinetic_smoothing_length','geometry_path','dense_moment_arithmetic',
          'native_budget_source','native_budget_transition_time')
    for key in keys:
        if configs[0].get(key)!=configs[1].get(key):raise ValueError('budget comparison model/source mismatch: '+key)
    for key in ('t','mass','state','factors'):
        if not np.array_equal(initials[0][key],initials[1][key]):raise ValueError('budget comparison initial array mismatch: '+key)


def field_difference(left,right,scales):
    result={}
    for key,scale in scales.items():
        delta=left[key]-right[key]
        norm=np.linalg.norm(delta,axis=1) if delta.ndim==2 else abs(delta)
        result[key]=float(np.mean(norm)/scale)
    return result


def endpoint_matches(t,end):
    return bool(abs(t-end)<=64*np.finfo(float).eps*max(abs(t),abs(end)))


def sampled_pair(cps,configs,points):
    from physical_field_accuracy_audit import environment
    from analyze_simple_wave_time import check_environment
    fields=[]
    for cp,config in zip(cps,configs):
        os.environ.update(environment(config));check_environment(config)
        fields.append(sampled_fields(cp,points))
    return fields


def run():
    install_domain();cfg=json.loads(Path('physical_budget_comparison_config.json').read_text());roots=list(map(Path,cfg['roots']))
    if len(roots)!=2:raise ValueError('two error-budget endpoints required')
    configs=[];reports=[];cps=[];initials=[];hashes=[]
    for root in roots:
        config=json.loads((root/'simple_wave_config.json').read_text());report=json.loads((root/'simple_wave_report.json').read_text());cp=json.loads((root/'checkpoint.json').read_text());initial=json.loads((root/'initial_state.json').read_text())
        identity=hashlib.sha256((root/'simple_wave_config.json').read_bytes()).hexdigest()
        if cp['config_sha256']!=identity or initial['config_sha256']!=identity:raise ValueError('budget comparison source identity mismatch')
        if not report['completed'] or not report['numerical_passed'] or not endpoint_matches(cp['t'],config['end']) or report['config']!=config:raise ValueError('budget endpoint incomplete or numerical gate failed')
        np.testing.assert_array_equal(cp['mass'],initial['mass'])
        configs.append(config);reports.append(report);cps.append(cp);initials.append(initial);hashes.append(hashlib.sha256((root/'checkpoint.json').read_bytes()).hexdigest())
    compatible(configs,initials)
    from quadratic_velocity_reconstruction import closure_name
    actual=closure_name() if os.environ.get('LAGFORCE_QUADRATIC_VELOCITY','0')=='1' else 'affine'
    from kinetic_speed_smoothing import parameters
    for c in configs:
        if c.get('velocity_reconstruction','affine')!=actual:raise ValueError('comparison velocity environment mismatch')
        from quadratic_pressure_reconstruction import closure_name as pressure_closure_name
        if c.get('pressure_reconstruction','linear')!=pressure_closure_name():raise ValueError('comparison pressure environment mismatch')
        if parameters()!=(c.get('kinetic_smoothing_coefficient',0.),c.get('kinetic_smoothing_length',1.)):raise ValueError('comparison smoothing environment mismatch')
    scales,goals=scales_and_goals(configs[0]);rows=[]
    for n in cfg.get('probe_grids',[128,256]):
        # Interior irrational offsets avoid symmetry faces without random probes.
        x,y=np.meshgrid((np.arange(n)+.3819660112501051)/n,(np.arange(n)+.6180339887498949)/n)
        points=np.column_stack([x.ravel(),y.ravel()]);fields=sampled_pair(cps,configs,points)
        delta=field_difference(*fields,scales)
        rows.append(dict(grid=n,points=len(points),differences_over_physical_scale=delta,fraction_of_physical_goal={k:delta[k]/goals[k] for k in goals}))
    screen=cfg.get('fraction_of_goal_limit',.05)
    result=dict(completed=True,comparison_kind=cfg.get('comparison_kind','path_error_budget'),roots=list(map(str,roots)),checkpoint_sha256=hashes,n=configs[0]['n'],t=cps[0]['t'],same_model_and_initial_arrays=True,
        fixed_mass_arrays_exact=True,all_numerical_gates_passed=True,probe_comparisons=rows,physical_goals=goals,fraction_of_goal_limit=screen,
        sampled_budget_sensitivity_passed=bool(max(v for row in rows for v in row['fraction_of_physical_goal'].values())<screen),
        physical_accuracy_passed=[r['physical_accuracy_passed'] for r in reports],
        note='Saved completed fields only; no dynamical stage repeated. Same material model and initial arrays, explicit numerical budgets or timestep caps differ as identified in the source configurations. Each endpoint uses its own recorded geometry/pressure evaluation backend. Two Eulerian probe grids estimate L1 differences; this is not a rigorous spatial quadrature enclosure or proof that either solution meets its exact-field accuracy goal.')
    atomic('physical_budget_comparison_report.json',result);print(result,flush=True)


if __name__=='__main__':run()
