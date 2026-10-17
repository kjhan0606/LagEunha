"""Check phase-scale sensitivity on saved steps without replaying evolution."""
import json,hashlib
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import numpy as np
from shapely.geometry import Polygon
from adaptive_coupled_mass import atomic,restore
from coupled_factor_path_audit import dense_state
from adaptive_moment_batch import setup
import moment_ale_evolve as model


def initialize():setup('shear_contact',3,100.)


def classify(t,y,mass,scales):
    a=y[:-4].reshape(-1,20);cells,_=model.geometry(a[:,:2],a[:,2],t)
    means=np.array([m/Polygon(v).area for m,v in zip(mass,cells)])
    low=float(means.min());high=float(means.max());gap=np.minimum(means-low,high-means)
    roundoff=32*np.finfo(float).eps*max(1.,abs(low),abs(high))
    minimum=max(roundoff,min(scales)*(high-low));maximum=max(roundoff,max(scales)*(high-low))
    flat=gap<=minimum;mixed=gap>=4*maximum;unresolved=~(flat|mixed)
    return dict(passed=not np.any(unresolved),unresolved_cells=np.flatnonzero(unresolved),
        gap_fraction=gap/(high-low) if high>low else np.zeros(len(gap)),flat_cells=int(np.sum(flat)))


def inspect(entry):
    payload=Path(entry['file']).read_bytes()
    if hashlib.sha256(payload).hexdigest()!=entry['sha256']:raise RuntimeError('saved step changed')
    r=restore(json.loads(payload)['record']);mass=np.asarray(entry['mass']);failures=[];flat_counts=[]
    for i in range(65):
        state,_,_=dense_state(r['y'],r['L'],r['d0'],r['dm'],r['diag'],mass,r['dt'],i/64)
        result=classify(r['t']+r['dt']*i/64,state,mass,entry['scales']);flat_counts.append(result['flat_cells'])
        if not result['passed']:failures.append(dict(kind='dense',index=i,**result))
    # Factor rebuilding changes Pi only, so actual midpoint m/q/Q/x/w are
    # exactly those of y + dt*d0/2. No RHS or moment-fit/LP is repeated here.
    half=r['y']+.5*r['dt']*r['d0'];result=classify(r['t']+.5*r['dt'],half,mass,entry['scales'])
    if not result['passed']:failures.append(dict(kind='actual_midpoint',**result))
    return dict(case=entry['case'],interval=entry['interval'],step=entry['step'],passed=not failures,
        dense_samples=65,actual_midpoint_stages=1,flat_cell_count_range=[min(flat_counts),max(flat_counts)],failures=failures)


def run():
    cfg=json.loads(Path('scale_audit_config.json').read_text());root=Path(cfg['source'])
    if not json.loads((root/'evolution_report.json').read_text())['passed']:raise RuntimeError('source evolution incomplete/failed')
    manifest=[]
    for case in cfg['cases']:
        original=json.loads((root/'inputs'/case/'checkpoint.json').read_text());mass=original['mass']
        primary=original['config']['density_provider']['phase_fraction_floor']
        if primary not in cfg['scales']:raise RuntimeError('primary scale omitted')
        for interval in ('short','extended'):
            path=root/'evolution'/case/interval;cp=json.loads((path/'checkpoint.json').read_text())
            if not json.loads((path/'result.json').read_text())['passed']:raise RuntimeError('failed interval')
            np.testing.assert_array_equal(cp['mass'],mass)
            for step in range(1,cp['steps']+1):
                p=path/('step_'+str(step)+'.json')
                manifest.append(dict(case=case,interval=interval,step=step,file=str(p),mass=mass,scales=cfg['scales'],sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    atomic('scale_manifest.json',manifest)
    with ProcessPoolExecutor(max_workers=2,mp_context=mp.get_context('spawn'),initializer=initialize) as executor:
        results=list(executor.map(inspect,manifest))
    passed=all(r['passed'] for r in results)
    report=dict(passed=passed,steps=len(results),dense_samples=sum(r['dense_samples'] for r in results),
        actual_midpoint_stages=sum(r['actual_midpoint_stages'] for r in results),scales=cfg['scales'],results=results,
        note='All sampled cells either flat at every selected scale or outside every blend band, with the same inferred moment normal. Therefore stage/path reconstruction equations coincide at these states for the tested positive scales. Actual midpoint m/q/Q checked separately. No RHS/LP/evolution replay and no continuous-path or epsilon->0 convergence proof.')
    atomic('scale_equivalence_report.json',report);print('scale equivalence',passed,'steps',len(results),flush=True)
    raise SystemExit(0 if passed else 2)


if __name__=='__main__':run()
