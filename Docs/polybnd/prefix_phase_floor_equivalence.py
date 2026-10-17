"""Check whether preserved accepted prefix samples activate the new blend.
Only x,w geometry is needed; fluid moments/factors and evolved state unchanged.
"""
import json,time,hashlib
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import numpy as np
from shapely.geometry import Polygon
from adaptive_moment_batch import setup
from adaptive_coupled_mass import atomic
import moment_ale_evolve as model


def initialize():setup('shear_contact',3,10.)


def inspect(entry):
    path=Path(entry['file']);data=path.read_bytes()
    if hashlib.sha256(data).hexdigest()!=entry['sha256']:raise RuntimeError('prefix input changed')
    saved=json.loads(data);r=saved['record'];y=np.asarray(r['y']);d0=np.asarray(r['d0']);dm=np.asarray(r['dm']);mass=np.asarray(entry['mass'])
    active=[];minimum=1e300
    for node in range(65):
        s=node/64;a=(y+r['dt']*(s*d0+s*s*(dm-d0)))[:-4].reshape(-1,20)
        polys,_=model.geometry(a[:,:2],a[:,2],r['t']+s*r['dt'])
        means=np.array([m/Polygon(v).area for m,v in zip(mass,polys)]);low=float(means.min());high=float(means.max())
        threshold=32*np.finfo(float).eps*max(1.,abs(low),abs(high));gap=np.minimum(means-low,high-means)
        changing=(gap>threshold)&(gap<4*threshold);positive=gap[gap>threshold]
        if len(positive):minimum=min(minimum,float(positive.min()/threshold))
        if np.any(changing):active.append(dict(node=node,cells=np.flatnonzero(changing).tolist(),ratios=(gap[changing]/threshold).tolist()))
    return dict(case=entry['case'],interval=entry['interval'],step=entry['step'],passed=not active,nodes=65,min_nonflat_gap_over_threshold=minimum,active=active)


def run():
    manifest=json.loads(Path('prefix_manifest.json').read_text());started=time.monotonic();rows=[]
    with ProcessPoolExecutor(max_workers=2,mp_context=mp.get_context('spawn'),initializer=initialize) as executor:
        for row in executor.map(inspect,manifest):
            rows.append(row)
            if len(rows)%20==0 or len(rows)==len(manifest):
                atomic('prefix_progress.json',dict(done=len(rows),total=len(manifest),results=rows));print('prefix',len(rows),'of',len(manifest),'active',sum(not r['passed'] for r in rows),flush=True)
    report=dict(passed=all(r['passed'] for r in rows),steps=len(rows),samples=sum(r['nodes'] for r in rows),results=rows,seconds=time.monotonic()-started,
        note='If no sample activates tau<gap<4tau, old and optional smooth prior formulas coincide at all preserved65-node geometry samples. This preserves previous fluid/path validation without repeating LP/flux evolution. Sampled equivalence, not an unsampled continuous proof.')
    atomic('prefix_equivalence_report.json',report)
    raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()
