"""Actual explicit-midpoint stage is distinct from the quadratic dense midpoint."""
import json,hashlib
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import numpy as np
from shapely.geometry import Polygon
from prefix_phase_floor_equivalence import initialize
from adaptive_coupled_mass import atomic
import moment_ale_evolve as model


def inspect(entry):
    data=Path(entry['file']).read_bytes()
    if hashlib.sha256(data).hexdigest()!=entry['sha256']:raise RuntimeError('prefix input changed')
    r=json.loads(data)['record'];mass=np.asarray(entry['mass'])
    a=(np.asarray(r['y'])+.5*r['dt']*np.asarray(r['d0']))[:-4].reshape(-1,20)
    polys,_=model.geometry(a[:,:2],a[:,2],r['t']+.5*r['dt'])
    means=np.array([m/Polygon(v).area for m,v in zip(mass,polys)]);low=float(means.min());high=float(means.max())
    threshold=32*np.finfo(float).eps*max(1.,abs(low),abs(high));gap=np.minimum(means-low,high-means)
    active=(gap>threshold)&(gap<4*threshold);positive=gap[gap>threshold]
    return dict(case=entry['case'],interval=entry['interval'],step=entry['step'],passed=not np.any(active),
        active_cells=np.flatnonzero(active),min_nonflat_gap_over_threshold=float(positive.min()/threshold) if len(positive) else None)


def run():
    dense=json.loads(Path('prefix_equivalence_report.json').read_text());manifest=json.loads(Path('prefix_manifest.json').read_text())
    with ProcessPoolExecutor(max_workers=2,mp_context=mp.get_context('spawn'),initializer=initialize) as executor:rows=list(executor.map(inspect,manifest))
    report=dict(passed=bool(dense['passed'] and all(r['passed'] for r in rows)),stages=len(rows),dense_samples=dense['samples'],results=rows,
        note='Actual RK midpoint geometries checked in addition to the65-node quadratic dense path; no fluid step or LP repeated.')
    atomic('prefix_stage_report.json',report);print('combined prefix equivalence',report['passed'],'stages',len(rows),'dense samples',dense['samples'],flush=True)
    raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()
