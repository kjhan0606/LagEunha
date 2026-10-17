"""Compare dt,dt/2,dt/4 from one saved state; reuse the completed coarse run."""
import json,hashlib,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import numpy as np
from adaptive_coupled_mass import atomic
from predictive_evolution_batch import advance_interval,initialize_case,density_options
from common_interface_contact import initialize
from common_contact_evolution import reference_moments
from adaptive_moment_batch import diagnostic


def digest(value):return hashlib.sha256(np.asarray(value,dtype=float).tobytes()).hexdigest()


def difference(a,b,mass):
    aa=np.asarray(a['state'])[:-4].reshape(-1,20);bb=np.asarray(b['state'])[:-4].reshape(-1,20)
    delta=aa-bb;w=delta[:,2]-delta[:,2].mean()
    return dict(generator=float(np.max(abs(delta[:,:2]))),gauge_fixed_weight=float(np.max(abs(w))),
        raw_moments_over_mass=float(np.max(abs(delta[:,3:])/mass[:,None])),ledger=float(np.max(abs(np.asarray(a['state'])[-4:]-np.asarray(b['state'])[-4:]))))


def uniform_grid(checkpoint,step_sizes,duration,nominal_steps):
    """A rejection-free adaptive run can still shrink its accepted step sizes."""
    return bool(checkpoint.get('rejected',0)==0 and checkpoint['steps']==nominal_steps
        and len(step_sizes)==nominal_steps
        and np.allclose(step_sizes,duration/nominal_steps,rtol=1e-12,atol=0.))


def observed_orders(differences,floor,uniform):
    return {key:(float(np.log2(differences[0][key]/differences[1][key]))
        if uniform and min(differences[0][key],differences[1][key])>floor else None)
        for key in differences[0]}


def end_time_agrees(actual,target,duration):
    tolerance=64*np.finfo(float).eps*max(abs(target),abs(duration),np.finfo(float).tiny)
    return bool(np.isfinite(actual) and abs(actual-target)<=tolerance)


def endpoint(cp,source,cfg):
    import moment_ale_evolve as model
    mass=np.asarray(cp['mass']);t=cp['t'];y=np.asarray(cp['state']);a=y[:-4].reshape(-1,20)
    _,_,refs,phases=initialize(cfg['n'],cfg['amplitude'],cfg['angle'])
    cells,_=model.geometry(a[:,:2],a[:,2],t);exact,coverage,_=reference_moments(cells,refs,phases,t)
    raw=np.column_stack([mass,a[:,3:9]]);corrected=y.copy();corrected[-4:]-=np.asarray(source['state'])[-4:]
    d=diagnostic(t,corrected,mass,np.asarray(source['state']))
    return dict(cell_mass_reference_L1=float(np.sum(abs(raw[:,0]-exact[:,0]))/mass.sum()),
        q_reference_L1=float(np.sum(abs(raw[:,1:3]-exact[:,1:3]))/mass.sum()),
        Q_reference_L1=float(np.sum(abs(raw[:,3:7]-exact[:,3:7]))/mass.sum()),
        reference_coverage=coverage,pressure_error=d['pressure_departure'],ledger_error=d['ledger_error'])


def analyze_saved(config,base=Path('.'),seconds=0.):
    """Read completed integration artifacts without re-running any step."""
    base=Path(base)
    for path,expected in config['source_sha256'].items():
        if hashlib.sha256((base/path).read_bytes()).hexdigest()!=expected:
            raise RuntimeError('time source changed: '+path)
    captured=json.loads((base/config['source_checkpoint']).read_text())
    initial=json.loads((base/config['initial_config_checkpoint']).read_text())
    cfg=initial['config'];original=dict(captured,config=cfg);mass,_=initialize_case(cfg)
    np.testing.assert_array_equal(initial['mass'],mass)
    np.testing.assert_array_equal(original['mass'],mass)
    coarse_path=base/config['coarse_checkpoint'];coarse=json.loads(coarse_path.read_text())
    roots=[coarse_path.parent]+[base/'evolution'/('steps'+str(n))/'short' for n in config['fine_steps']]
    result_paths=[base/config['coarse_result']]+[root/'result.json' for root in roots[1:]]
    results=[json.loads(path.read_text()) for path in result_paths]
    if not all(result['passed'] for result in results):raise RuntimeError('time evolution did not pass')
    checkpoints=[coarse]+[json.loads((root/'checkpoint.json').read_text()) for root in roots[1:]]
    if len(checkpoints)!=3:raise ValueError('three dt levels required')
    duration=float(config['duration']);target=original['t']+duration
    for cp in checkpoints:
        np.testing.assert_array_equal(cp['mass'],mass)
        if not end_time_agrees(cp['t'],target,duration):raise RuntimeError('end times differ')
        if cp['config']['source_state_sha256']!=digest(original['state']):raise RuntimeError('source state differs')
        if cp['config']['density_provider']!=density_options(cfg,True):raise RuntimeError('model differs')
        if cp['config'].get('face_order',4)!=coarse['config'].get('face_order',4):raise RuntimeError('face quadrature differs')
    differences=[difference(checkpoints[i],checkpoints[i+1],mass) for i in range(2)]
    scales=max(1.,max(float(np.max(abs(np.asarray(cp['state'])[:-4].reshape(-1,20)[:,3:])/mass[:,None])) for cp in checkpoints))
    floor=512*np.finfo(float).eps*scales
    nominal=[coarse['config']['initial_steps']]+config['fine_steps']
    sizes=[[json.loads((root/('step_'+str(i)+'.json')).read_text())['dt']
        for i in range(1,cp['steps']+1)] for root,cp in zip(roots,checkpoints)]
    uniform=all(uniform_grid(cp,grid,duration,n) for cp,grid,n in zip(checkpoints,sizes,nominal))
    halving=nominal==[nominal[0],2*nominal[0],4*nominal[0]]
    order=observed_orders(differences,floor,uniform and halving)
    rows=[endpoint(cp,original,cfg) for cp in checkpoints]
    return dict(completed=True,passed=bool(all(r['reference_coverage']<1e-11 and r['ledger_error']<1e-10 and r['pressure_error']<1e-8 for r in rows)),
        source_time=original['t'],end_time=target,actual_end_times=[cp['t'] for cp in checkpoints],coarse_reused=True,
        nominal_steps=nominal,actual_steps=[cp['steps'] for cp in checkpoints],uniform_without_rejection=uniform,
        accepted_dt_ranges=[[min(grid),max(grid)] for grid in sizes],differences=differences,
        conservative_roundoff_screen=floor,observed_orders_above_screen=order,endpoint_reference_errors=rows,seconds=seconds,
        note='Saved integrations only, same source/model/grid/face bins. Actual time grids verified. Endpoint agreement permits only64*eps*time-scale roundoff. Order withheld below floating-point screen or for nonuniform/nonhalving grids. Sampled path gates are not a continuous bound.')


def run():
    started=time.monotonic();config=json.loads(Path('time_config.json').read_text())
    for path,expected in config['source_sha256'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=expected:raise RuntimeError('time source changed: '+path)
    captured=json.loads(Path(config['source_checkpoint']).read_text())
    initial=json.loads(Path(config['initial_config_checkpoint']).read_text())
    cfg=initial['config'];original=dict(captured,config=cfg);mass,_=initialize_case(cfg)
    # Interval checkpoint config is provenance, not the original IC config.
    # Preserve both original files verbatim and assemble only explicit inputs.
    Path('captured_source_checkpoint.json').write_bytes(Path(config['source_checkpoint']).read_bytes())
    Path('captured_initial_config_checkpoint.json').write_bytes(Path(config['initial_config_checkpoint']).read_bytes())
    np.testing.assert_array_equal(initial['mass'],mass)
    np.testing.assert_array_equal(original['mass'],mass)
    coarse_path=Path(config['coarse_checkpoint']);coarse=json.loads(coarse_path.read_text())
    coarse_result=json.loads(Path(config['coarse_result']).read_text())
    if not coarse_result['passed']:raise RuntimeError('coarse run failed')
    if coarse['config']['source_state_sha256']!=digest(original['state']):raise RuntimeError('coarse start differs')
    if coarse['config']['density_provider']!=density_options(cfg,True):raise RuntimeError('coarse density model differs')
    duration=float(config['duration']);tasks=[]
    for steps in config['fine_steps']:
        name='steps'+str(steps);path=Path('inputs')/name/'checkpoint.json'
        if path.exists():raise RuntimeError('refuse duplicate time experiment input')
        path.parent.mkdir(parents=True)
        atomic(path,dict(original,provenance=dict(source=config['source_checkpoint'],source_sha256=hashlib.sha256(Path(config['source_checkpoint']).read_bytes()).hexdigest(),purpose='Time refinement from accepted saved state; fixed density model, mass, grid and face bins.')))
        tasks.append((name,'short',duration,steps))
    with ProcessPoolExecutor(max_workers=2,mp_context=mp.get_context('spawn')) as executor:
        results=list(executor.map(advance_interval,tasks))
    atomic('time_evolution_report.json',dict(passed=all(r['passed'] for r in results),results=results))
    if not all(r['passed'] for r in results):raise SystemExit(2)
    report=analyze_saved(config,seconds=time.monotonic()-started)
    atomic('time_refinement_report.json',report);print('time refinement',report['passed'],report['differences'],report['observed_orders_above_screen'],flush=True)
    raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--analyze-source',type=Path);args=parser.parse_args()
    if args.analyze_source is None:run()
    else:
        config=json.loads((args.analyze_source/'time_config.json').read_text())
        report=analyze_saved(config,args.analyze_source);report['analysis_source']=str(args.analyze_source)
        atomic('time_refinement_report.json',report)
        print(json.dumps(report),flush=True)
        raise SystemExit(0 if report['passed'] else 2)
