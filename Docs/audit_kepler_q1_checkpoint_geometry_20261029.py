"""Independent saved-state geometry and circular-orbit density attribution.

Read-only postprocessing on a compute host; never a hydrodynamic restart.
"""
import argparse
import hashlib
import json
from pathlib import Path
import socket
import time
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from make_kepler_bounded_polygon_ic_20261027 import bounded_areas
from audit_kepler_native_radial_20261026 import mean


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent(xy):
    a, _, n4 = bounded_areas(xy, 4.)
    b, _, n5 = bounded_areas(xy, 5.)
    domain = .5 * 1024 * 2.75**2 * np.sin(2*np.pi/1024)
    agreement = float(np.max(abs(a/b-1)))
    partition = abs(float(np.sum(b, dtype=np.longdouble))-domain)
    if agreement >= 1e-8 or partition >= 1e-7:
        raise ValueError('Independent guard-radius/area partition failed')
    return b, {'guard_radius_max_relative_difference': agreement,
               'domain_area_absolute_error': partition, 'clipped_cells': [n4, n5]}


def audit(initial, endpoint):
    if socket.gethostname().split('.')[0].lower() in ('grammar', 'syntax'):
        raise RuntimeError('Run geometry analysis on a compute host')
    started = time.monotonic()
    p0, m0 = ordered_global(initial); ids = p0['u4if'] & ID_MASK
    p, mt = ordered_global(endpoint, ids)
    if m0['t'] != 0 or mt['t'] <= 0 or len(p) != 90352:
        raise ValueError('Requires the bounded Q1 IC and a complete evolved state')
    if not np.array_equal(p0['mass'], p['mass']) or any(np.any(q['w2'] != 0) for q in (p0, p)):
        raise ValueError('This diagnostic requires fixed mass and W0')
    for q in (p0, p):
        if not all(np.isfinite(q[k]).all() for k in ['x', 'y', 'vx', 'vy', 'mass', 'volume', 'den']):
            raise ValueError('Nonfinite native state')
        if not all(np.all(q[k] > 0) for k in ['mass', 'volume', 'den']):
            raise ValueError('Nonpositive native state')
    x = p0['x']-3.; y = p0['y']-3.; r0 = np.hypot(x, y)
    phase = mt['t'] * (r0*r0+.01175**2)**(-.75)
    circular = np.column_stack((3.+x*np.cos(phase)-y*np.sin(phase),
                                3.+x*np.sin(phase)+y*np.cos(phase)))
    native_xy = np.column_stack((p['x'], p['y']))
    native_area, native_check = independent(native_xy)
    circular_area, circular_check = independent(circular)
    volume = p['volume'].astype(float); volume0 = p0['volume'].astype(float)
    geometry = volume/native_area-1
    measures = {'native_density_change': p['den'].astype(float)/p0['den'].astype(float)-1,
                'circular_geometry_density_change': volume0/circular_area-1,
                'native_density_relative_to_circular_geometry': circular_area/volume-1,
                'native_volume_relative_to_independent_geometry': geometry}
    bands = [(0., .1), (.1, .5), (.5, .75), (.75, 1.), (1., 2.), (.5, 2.), (2., 2.65)]
    regions = []
    for lo, hi in bands:
        mask = (r0 >= lo) & (r0 < hi); mass = p0['mass'][mask].astype(float)
        if not np.any(mask): raise ValueError('Empty reference region')
        values = {}
        for name, data in measures.items():
            a = data[mask]
            values[name] = {'signed_mass_mean': mean(a, mass),
                            'absolute_mass_mean': mean(abs(a), mass),
                            'absolute_cell_p99': float(np.percentile(abs(a), 99)),
                            'absolute_max': float(np.max(abs(a)))}
        regions.append({'rlo': lo, 'rhi': hi, 'cells': int(mask.sum()), 'measures': values})
    maximum = float(np.max(abs(geometry)))
    worst = np.argsort(abs(geometry))[-10:][::-1]
    return {'passed': maximum < 1e-6,
            'checks': {'independent_guards_and_partitions': True,
                       'native_volume_matches_independent_below_1e_minus6': maximum < 1e-6},
            'time': mt['t'], 'sites': len(p), 'native_geometry': native_check,
            'circular_reference_geometry': circular_check,
            'native_volume_max_relative_error': maximum, 'regions': regions,
            'worst_native_geometry_IDs': [int(v) for v in ids[worst]],
            'source_hashes': {str(q): digest(q) for q in (initial, endpoint, Path(__file__),
                Path(__file__).with_name('make_kepler_bounded_polygon_ic_20261027.py'))},
            'analysis_seconds': time.monotonic()-started,
            'scope': 'All-site independent SciPy Voronoi areas clipped to the same fixed polygon, with two guard radii. Circular reference rotates the same original sites at their Plummer circular frequencies and keeps their masses. It isolates geometric density breathing and dynamical departures; it is not gas evolution, entropy attribution, an exact finite-cell fluid solution, boundary convergence, or a long-time promotion.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--initial', type=Path, required=True)
    parser.add_argument('--endpoint', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): raise FileExistsError('Refusing to overwrite independent evidence')
    result = audit(args.initial, args.endpoint)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2); stream.write('\n')
    print(json.dumps(result, indent=2))
    if not result['passed']: raise SystemExit(1)
