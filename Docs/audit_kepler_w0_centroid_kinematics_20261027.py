"""Finite-difference area-centroid velocity on saved native W0 Kepler states."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import Voronoi

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest

N = 55351


def centroids(points, active):
    vor = Voronoi(points)
    result = np.empty((len(active), 2))
    for k, i in enumerate(active):
        region = vor.regions[vor.point_region[i]]
        if len(region) < 3 or -1 in region:
            raise ValueError(f'Unbounded active cell {i}')
        # Local coordinates avoid cancellation for small cells.
        v = vor.vertices[region] - points[i]
        nxt = np.roll(v, -1, axis=0)
        cross = v[:, 0]*nxt[:, 1] - nxt[:, 0]*v[:, 1]
        twice_area = cross.sum()
        if abs(twice_area) < 1e-14:
            raise ValueError(f'Degenerate active cell {i}')
        result[k] = points[i] + np.sum((v+nxt)*cross[:, None], axis=0)/(3*twice_area)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('receipt', type=Path)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    grid = np.array([(x, y) for x in range(-2, 3) for y in range(-2, 3)], dtype=float)
    common_velocity = np.array([.3, -.7])
    common_fd = (centroids(grid+1e-5*common_velocity, [12])
                 -centroids(grid-1e-5*common_velocity, [12]))/(2e-5)
    if np.max(np.abs(common_fd-common_velocity)) > 1e-9:
        raise ValueError('Uniform-translation centroid test failed')
    receipt = json.loads(args.receipt.read_text())
    if not all(receipt['checks'].values()):
        raise ValueError('Native parent audit did not pass')
    initial, endpoint = (Path(receipt[key]) for key in ('initial', 'endpoint'))
    for key, path in (('initial', initial), ('endpoint', endpoint)):
        if digest(path) != receipt[key+'_sha256']:
            raise ValueError(f'Changed {key} state')
    p0, _ = ordered(initial, N)
    p, meta = ordered(endpoint, N)
    if abs(float(meta['t'])-receipt['time']) > 1e-10 or np.any(p['w2'] != 0):
        raise ValueError('Expected the accepted W0 endpoint at the receipt time')
    xy0 = np.column_stack((p0['x'].astype(float)-3., p0['y'].astype(float)-3.))
    xy = np.column_stack((p['x'].astype(float)-3., p['y'].astype(float)-3.))
    velocity = np.column_stack((p['vx'].astype(float), p['vy'].astype(float)))
    r0 = np.linalg.norm(xy0, axis=1)
    active = np.flatnonzero((r0 >= .5) & (r0 < 2.))
    mass = p['mass'][active].astype(float)
    eps = 1e-5
    derivative = []
    for h in (eps, eps/2):
        derivative.append((centroids(xy+h*velocity, active)
                           -centroids(xy-h*velocity, active))/(2*h))
    delta = derivative[1]-velocity[active]
    error = np.linalg.norm(derivative[1]-derivative[0], axis=1)
    mismatch = np.linalg.norm(delta, axis=1)
    speed = np.linalg.norm(velocity[active], axis=1)
    weighted = lambda value: float(np.average(value, weights=mass))
    rms = lambda value: float(np.sqrt(np.average(value*value, weights=mass)))
    bands = {}
    for low, high in ((.5, .75), (.75, 1.), (1., 1.5), (1.5, 2.)):
        selected = (r0[active] >= low) & (r0[active] < high)
        bands[f'{low:g}-{high:g}'] = {
            'count': int(selected.sum()),
            'mass_weighted_relative_rms': float(np.sqrt(np.average(
                (mismatch[selected]/speed[selected])**2, weights=mass[selected]))),
        }
    checks = {'uniform_translation': True, 'finite': bool(np.all(np.isfinite(delta))),
              'finite_difference_converged': bool(rms(error) < 1e-3*max(rms(mismatch), 1e-14)),
              'nonzero_centroid_velocity_mismatch': bool(rms(mismatch) > 1e-6)}
    result = {'receipt': str(args.receipt.resolve()), 'receipt_sha256': digest(args.receipt),
              'endpoint_sha256': receipt['endpoint_sha256'], 'time': receipt['time'],
              'disk_sites': len(active), 'epsilon': eps,
              'mass_weighted_rms_speed': rms(speed),
              'mass_weighted_rms_centroid_minus_generator_speed': rms(mismatch),
              'mass_weighted_relative_rms': rms(mismatch/speed),
              'mass_weighted_mean_relative': weighted(mismatch/speed),
              'finite_difference_rms_error_estimate': rms(error),
              'radial_bands_initial_radius': bands, 'checks': checks,
              'scope': 'Instantaneous W0 Voronoi geometry using saved generator velocities; '
                       'zero face mass flux and uniform-density mass-centroid identification '
                       'are a hypothetical closure, not asserted properties of the native solver. '
                       'This is not an evolution or proof of long-time physical error.'}
    with args.output.open('x') as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write('\n')
    print(json.dumps(result, indent=2, allow_nan=False))
    if not all(checks.values()):
        raise SystemExit('Centroid kinematic diagnostic checks failed')


if __name__ == '__main__':
    main()
