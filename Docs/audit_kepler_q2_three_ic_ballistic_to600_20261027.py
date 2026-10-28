"""Exact-orbit fixed-mass W0 Voronoi-area reference through t=600."""

import argparse
import json
from pathlib import Path
import time

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from audit_kepler_sparse_native_ballistic_20261027 import disk_geometry
from audit_kepler_q2_three_ic_native_t65_20261027 import SOURCES
from prepare_kepler_sparse_area_restart_20261027 import digest


TIMES = (6.5, 10., 30., 100., 600.)


def metrics(initial, area, take):
    v0 = initial['volume'][take].astype(float)
    m0 = initial['mass'][take].astype(float)
    selected_area = area[take]
    change = v0/selected_area-1.
    return {
        'cells': int(np.sum(take)),
        'mean_abs_density_change_mass_weighted': float(
            np.sum(np.abs(change)*m0, dtype=np.longdouble)
            / np.sum(m0, dtype=np.longdouble)),
        'p95_abs_density_change': float(np.percentile(np.abs(change), 95)),
        'max_abs_density_change': float(np.max(np.abs(change))),
        'min_area_ratio': float(np.min(selected_area/v0)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    started = time.monotonic()
    cases = {}
    for label, source in SOURCES.items():
        audit = json.loads(source.read_text())
        initial = Path(audit['initial'])
        if (audit['time'] != 6.5 or not all(audit['checks'].values())
                or digest(initial) != audit['initial_sha256']):
            raise ValueError(f'{label}: invalid accepted initial state')
        p0, _ = ordered(initial, 55351)
        r0 = np.hypot(p0['x'].astype(float)-3., p0['y'].astype(float)-3.)
        selected = np.flatnonzero((r0 >= .5) & (r0 < 2.))
        if len(selected) != 49038:
            raise ValueError(f'{label}: unexpected disk mask')
        by_time = {}
        for t in TIMES:
            _, area, safety = disk_geometry(p0, t, selected)
            if not np.all(np.isfinite(area) & (area > 0)):
                raise ValueError(f'{label} t={t}: invalid exact-orbit areas')
            disk_initial = p0[selected]
            disk_r0 = r0[selected]
            by_time[str(t)] = {
                'disk': metrics(disk_initial, area, np.ones(len(selected), dtype=bool)),
                'inner': metrics(disk_initial, area, disk_r0 < .75),
                'minimum_periodic_image_safety_ratio': safety,
            }
        expected = audit['disk_ballistic_density_abs_mass_weighted']
        actual = by_time['6.5']['disk']['mean_abs_density_change_mass_weighted']
        if abs(actual-expected) > 1e-8:
            raise ValueError(f'{label}: t6.5 geometric anchor changed: {actual} vs {expected}')
        cases[label] = {'initial': str(initial),
                        'initial_sha256': audit['initial_sha256'],
                        'by_time': by_time}
    result = {'times': TIMES, 'cases': cases,
              'elapsed_seconds': time.monotonic()-started,
              'scope': 'Exact softened circular motion of initial generators with original fixed masses and rebuilt W0 Voronoi areas. This is a kinematic/geometric reference, not a strict lower bound on all numerical trajectories or a native t=600 hydro stability proof.'}
    with args.output.open('x') as out:
        json.dump(result, out, indent=2)
        out.write('\n')
    print(json.dumps({label:{t:values['disk']['mean_abs_density_change_mass_weighted']
                             for t,values in case['by_time'].items()}
                      for label,case in cases.items()},indent=2))


if __name__=='__main__':
    main()
