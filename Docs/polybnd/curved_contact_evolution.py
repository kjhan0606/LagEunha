"""Continue the physical curved contact from a preserved accepted state."""
import hashlib,json,math,time
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from preserved_experiment_input import install_input
from predictive_evolution_batch import advance_interval,initialize_case
from curved_contact_first_step import load_saved_rhs
from saved_curved_contact import checked_json,score


def run():
    started=time.monotonic();entry=json.loads(Path('curved_evolution_config.json').read_text())
    source=checked_json(entry['source'],entry['source_sha256'])
    initial=checked_json(entry['initial_state'],entry['initial_state_sha256'])
    cfg=dict(kind='saved_curved_contact',n=initial['n'],amplitude=100.,smooth_phase_floor=True,
        initial_state=entry['initial_state'],initial_state_sha256=entry['initial_state_sha256'],
        reference=entry['reference'],reference_sha256=entry['reference_sha256'],
        density_provider=dict(mode='moment',segments=64,radial=1e-7,core_level=3,align_core=True,
            vector_geometry=True,construct_prior_witness=True,lp_retry=True,phase_fraction_floor=1e-7),
        face_order=4,refined_path_audit=True,adaptive_path_quadrature=True,
        max_path_quadrature_nodes=4097,path_control_horizon=entry['control_horizon'])
    cfg['density_provider'].update(entry.get('density_provider',{}))
    mass,_=initialize_case(cfg);np.testing.assert_array_equal(mass,source['mass'])
    if source['n']!=initial['n'] or not 0<=source['t']<entry['target']:
        raise ValueError('curved continuation source identity/time mismatch')
    name=entry['case'];duration=entry['target']-source['t']
    steps=max(1,math.ceil(duration/entry['max_dt']-1e-12))
    cp=dict(t=source['t'],mass=mass,state=source['state'],factors=source['factors'],config=cfg,
        provenance=dict(source=entry['source'],source_sha256=entry['source_sha256'],approved_prefix_reused=True))
    install_input(Path('inputs')/name/'checkpoint.json',cp)
    rhs=None
    result_path=Path('evolution')/name/'short'/'checkpoint.json'
    if entry.get('source_rhs') and not result_path.exists():
        if source!=initial:raise ValueError('initial RHS may only accompany its exact saved initial state')
        rhs_data=checked_json(entry['source_rhs'],entry['source_rhs_sha256'])
        gate=checked_json(entry['source_gate'],entry['source_gate_sha256'])
        if not gate['passed'] or gate['n']!=initial['n']:raise ValueError('saved initial RHS gate failed')
        saved=load_saved_rhs(rhs_data);rhs=(saved['dy'],saved['diag'])
    result=advance_interval((name,'short',duration,steps),initial_rhs=rhs)
    report=dict(completed=True,passed=bool(result['passed']),result=result,source_time=source['t'],
        source_sha256=entry['source_sha256'],source_rhs_reused=rhs is not None,
        note='Original curved mass and accepted prefix retained. Exact advected128-gon is used only for scoring; numerical RHS/weights receive moments only.')
    if result['passed']:
        final=json.loads(result_path.read_text());report['final_score']=score(final,cfg)
        report['passed']=bool(report['passed'] and report['final_score']['reference_valid'])
    report['seconds']=time.monotonic()-started
    atomic('curved_evolution_report.json',report)
    print('curved continuation',name,report['passed'],flush=True)
    raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()
