"""Join an equivalent validated prefix to the tested continuous-floor event.
No fluid state is recomputed, projected, or changed in this finalization.
"""
import json,hashlib
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import setup,diagnostic
from shear_contact_check import affine_errors


def read(path):return json.loads(Path(path).read_text())


def run():
    dense=read('predecessors/prefix_equivalence_report.json');stages=read('predecessors/prefix_stage_report.json')
    event_result=read('predecessors/event_result.json');prefix=read('predecessors/prefix_checkpoint.json');event=read('predecessors/event_checkpoint.json')
    initial_event=read('predecessors/event_initial_checkpoint.json')
    if not dense['passed'] or not stages['passed'] or not event_result['passed']:raise RuntimeError('predecessor gate failed')
    if not initial_event['config']['smooth_phase_floor'] or event_result['path_control_horizon']!=1e-4:raise RuntimeError('event physics/control mismatch')
    for key in ('state','mass','factors'):np.testing.assert_array_equal(prefix[key],initial_event[key])
    if prefix['t']!=initial_event['t'] or event['t']!=prefix['config']['end']:raise RuntimeError('checkpoint interval mismatch')
    start_sha=hashlib.sha256(np.asarray(prefix['state']).tobytes()).hexdigest()
    if event['config']['source_state_sha256']!=start_sha:raise RuntimeError('event source state mismatch')
    expected={(row['case'],row['interval'],row['step']) for row in dense['results']}
    actual={(row['case'],row['interval'],row['step']) for row in stages['results']}
    if actual!=expected:raise RuntimeError('prefix stage/dense manifests differ')
    if sum(c=='n4_r100' and i=='extended' for c,i,_ in expected)!=prefix['steps']:raise RuntimeError('wrong prefix length')
    results=[]
    for name in ('n3_r10','n3_r100','n4_r100'):
        source=read(Path('inputs')/name/'checkpoint.json');mass=np.asarray(source['mass'])
        (Path('endpoints')/name).mkdir(parents=True,exist_ok=True)
        for interval in ('short','extended'):
            if (name,interval)==('n4_r100','extended'):continue
            cp=read(Path('cases')/name/interval/'checkpoint.json');r=read(Path('cases')/name/interval/'result.json')
            np.testing.assert_array_equal(cp['mass'],mass)
            if not r['passed'] or cp['t']!=r['t']:raise RuntimeError('incomplete inherited case')
            if cp['steps']!=sum(c==name and i==interval for c,i,_ in expected):raise RuntimeError('inherited interval not covered by equivalence check')
            atomic(Path('endpoints')/name/(interval+'.json'),cp);results.append(dict(r,prefix_equivalence_verified=True))
    name='n4_r100';source=read(Path('inputs')/name/'checkpoint.json');mass=np.asarray(source['mass']);cfg=source['config']
    short=read(Path('endpoints')/name/'short.json')
    np.testing.assert_array_equal(prefix['mass'],mass);np.testing.assert_array_equal(event['mass'],mass)
    if prefix['config']['source_state_sha256']!=hashlib.sha256(np.asarray(short['state']).tobytes()).hexdigest():raise RuntimeError('wrong extended-prefix origin')
    setup('shear_contact',cfg['n'],cfg['amplitude'])
    y=np.asarray(event['state']);relative=y.copy();relative[-4:]-=np.asarray(short['state'])[-4:]
    d=diagnostic(event['t'],relative,mass,np.asarray(short['state']));affine=affine_errors(y,mass)
    accum={key:np.asarray(prefix[key])+event[key] for key in ('signed','absolute','quad_delta')}
    peak={key:max(prefix['peak'][key],event['peak'][key]) for key in ('rate','generator','pressure_error','stage_mass')}
    peak.update(ledger=max(prefix['peak']['ledger']+event['peak']['ledger'],d['ledger_error']),
        min_pressure=min(prefix['peak']['min_pressure'],event['peak']['min_pressure']),min_factor_R=min(prefix['peak']['min_factor_R'],event['peak']['min_factor_R']))
    absolute=float(np.max(accum['absolute']/mass));quad=float(np.max(abs(accum['quad_delta'])/mass))
    passed=bool(peak['rate']<1e-8 and absolute<1e-10 and quad<1e-11 and peak['generator']<1e-8 and peak['ledger']<1e-10 and peak['min_pressure']>0 and peak['min_factor_R']>=-1e-20 and d['pressure_departure']<1e-8 and max(affine.values())<1e-8)
    joined=dict(event,config={**prefix['config'],'smooth_phase_floor':True,'completion_method':'equivalent_prefix_plus_smoothed_event'},
        steps=prefix['steps']+event['steps'],rejected=prefix['rejected']+event['rejected'],peak=peak,**accum)
    atomic(Path('endpoints')/name/'extended.json',joined)
    result=dict(case=name,interval='extended',status='completed',passed=passed,t=event['t'],duration=prefix['config']['duration'],
        steps=joined['steps'],rejected=joined['rejected'],reused_prefix_steps=prefix['steps'],new_event_steps=event['steps'],
        integrated_absolute_mass_defect=absolute,integrated_signed_mass_defect=float(np.max(abs(accum['signed'])/mass)),quadrature_mass_difference=quad,peak=peak,
        final_ledger_error=d['ledger_error'],affine_errors=affine,completion_method='equivalent_prefix_plus_smoothed_event')
    results.append(result)
    report=dict(completed=len(results)==6,passed=len(results)==6 and all(r['passed'] for r in results),results=results,
        prefix_steps_checked=dense['steps'],prefix_dense_samples=dense['samples'],prefix_half_stages=stages['stages'],
        note='All inherited midpoint stages and65-node dense samples are outside the smoothing transition. Reused accepted states are unchanged. N4 continued from its actual step126 checkpoint with the separately validated smooth floor. This is sampled/stage equivalence, not an unsampled existence proof. Baseline AQ is preserved separately.')
    atomic('evolution_report.json',report);print(json.dumps(report),flush=True)
    raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()
