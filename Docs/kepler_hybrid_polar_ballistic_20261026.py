"""Exact-orbit Voronoi control for the full 384² paper-profile hybrid ICs.

One-off geometry analysis; no hydro evolution or native GPU performance claim.
"""

from pathlib import Path
import json
import struct
import time

import numpy as np
from scipy.spatial import Voronoi

N = 384
L = 6.
H = L / N
EPS = .01175
TIMES = (.825, 2.221441469079183, 120., 600.)
GLASS = {'polar_q1': Path('/tmp/kepler_hybrid_polar_q1_384.bin'),
         'polar_q2': Path('/tmp/kepler_hybrid_polar_q2_384.bin')}


def read_glass(path):
    with path.open('rb') as stream:
        n, lx, ly = struct.unpack('<idd', stream.read(20))
        if (n, lx, ly) != (N*N, L, L):
            raise ValueError('Unexpected glass header')
        x = np.fromfile(stream, dtype='<f8', count=n)
        y = np.fromfile(stream, dtype='<f8', count=n)
    if len(x) != n or len(y) != n:
        raise ValueError('Truncated glass positions')
    return np.column_stack((x, y))


def polygon_areas(xy, selected):
    vor = Voronoi(xy, qhull_options='Qbb Qc Qz')
    area = np.empty(len(selected))
    for k, i in enumerate(selected):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f'Unbounded annulus cell {i}')
        polygon = vor.vertices[region]
        area[k] = .5 * abs(np.dot(polygon[:, 0], np.roll(polygon[:, 1], -1)) -
                        np.dot(polygon[:, 1], np.roll(polygon[:, 0], -1)))
    return area


def main():
    started = time.monotonic()
    ids = np.arange(N*N)
    cartesian = np.column_stack(((ids % N + .5)*H, (ids // N + .5)*H))
    layouts = {'cartesian': cartesian}
    layouts.update({name: read_glass(path) for name, path in GLASS.items()})
    cases = []
    for name, xy0 in layouts.items():
        dx0, dy0 = xy0[:, 0]-3., xy0[:, 1]-3.
        radius = np.hypot(dx0, dy0)
        selected = np.flatnonzero((radius >= .5) & (radius < 2.))
        a0 = polygon_areas(xy0, selected)
        mass = 1.01*a0  # paper profile has constant Sigma=1.01 in the disk
        outputs = []
        for t in TIMES:
            if name == 'cartesian' and t > 2.3:
                continue
            angle = t * (radius*radius+EPS*EPS)**(-.75)
            xy = np.column_stack((3. + dx0*np.cos(angle)-dy0*np.sin(angle),
                                  3. + dx0*np.sin(angle)+dy0*np.cos(angle)))
            at = polygon_areas(xy, selected)
            error = abs(a0/at - 1.)
            bins = []
            for lo, hi in ((.5, .75), (.75, 1.), (1., 2.), (.5, 2.)):
                take = (radius[selected] >= lo) & (radius[selected] < hi)
                bins.append(dict(rlo=lo, rhi=hi, cells=int(np.sum(take)),
                                 density_change_mass_weighted=float(np.sum(mass[take]*error[take])/np.sum(mass[take])),
                                 density_change_p95=float(np.quantile(error[take], .95)),
                                 density_change_max=float(np.max(error[take]))))
            outputs.append(dict(time=t, bins=bins))
        cases.append(dict(layout=name, total_cells=len(xy0), disk_cells=len(selected),
                          initial_disk_cell_area_min=float(np.min(a0)),
                          initial_disk_cell_area_max=float(np.max(a0)),
                          results=outputs))
        print(f'finished {name}', flush=True)
    result = dict(description='Exact softened Kepler orbits followed by independent unweighted Voronoi cell areas; fixed initial cell masses',
                  box=L, grid=N, epsilon=EPS, cases=cases,
                  wall_seconds=time.monotonic()-started,
                  limitation='Geometry-only; no hydrodynamics, central stability, periodic wrapping of outer background, or GPU throughput. The hybrid boundary at r=2.35 is outside the reported disk annulus.')
    output = Path(__file__).with_name('KEPLER_HYBRID_POLAR_BALLISTIC_2026-10-26.json')
    with output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(output)


if __name__ == '__main__':
    main()
