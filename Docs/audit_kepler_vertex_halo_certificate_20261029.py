"""Read-only periodic W=0 geometry audit; never evolves the saved state."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
from scipy.spatial import Voronoi

sys.path.insert(0, '/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
from audit_snapshot import read_snapshot


def influence_bounds(owner, vertices):
    """W=0 potential cutters lie in the union of these vertex disks.

    A generator q cannot cut at a vertex v when |v-q| >= |v-owner|.
    A linear Voronoi half-plane satisfied at all vertices contains the
    complete convex polygon. Return an enclosing rectangle of the disks.
    """
    radii = np.linalg.norm(vertices - owner, axis=1)
    return (np.min(vertices - radii[:, None], axis=0),
            np.max(vertices + radii[:, None], axis=0), float(np.max(radii)))


def rectangle_margin(bounds, lower, upper):
    lo, hi = bounds[:2]
    return float(np.min(np.concatenate((lo - lower, upper - hi))))


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('snapshot', type=Path)
    parser.add_argument('--expected-sha256', required=True)
    parser.add_argument('--cellsize', type=float, default=.0625)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    started = time.perf_counter()
    assert digest(args.snapshot) == args.expected_sha256
    p, meta = read_snapshot(args.snapshot)
    xy = np.column_stack((p['x'], p['y'])).astype(np.float64)
    length = np.array([meta['Lx'], meta['Ly']], dtype=float)
    assert np.array_equal(length, [6., 6.]) and len(xy) == 384**2
    assert meta['psize'] == 576 and np.isfinite(xy).all()
    assert np.all(xy >= 0) and np.all(xy < length)
    assert np.all(p['w2'] == 0) and np.all(p['volume'] > 0)
    assert np.isfinite(args.cellsize) and args.cellsize == .0625
    ids = (p['u4if'].astype(np.uint64) & np.uint64((1 << 56) - 1)).astype(np.int64)
    assert np.array_equal(np.sort(ids), np.arange(len(xy)))
    assert len(np.unique(xy, axis=0)) == len(xy), 'Coincident generators need separate treatment'
    shifts = [(0, 0)] + [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1) if (i, j) != (0, 0)]
    all_xy = np.concatenate([xy + np.array(shift) * length for shift in shifts])
    vor = Voronoi(all_xy)  # default deterministic Qhull options; no coordinate jitter
    built = time.perf_counter()
    n = len(xy)
    areas = np.empty(n)
    old_extent = np.empty((n, 4))
    new_extent = np.empty((n, 4))
    radii = np.empty(n)
    periodic_margin = np.empty(n)
    vertices_by_cell = []
    for i in range(n):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f'Unbounded central periodic cell {int(ids[i])}')
        vertices = vor.vertices[region]
        assert np.isfinite(vertices).all() and len(vertices) >= 3
        local = vertices - xy[i]
        areas[i] = .5 * abs(np.sum(local[:, 0] * np.roll(local[:, 1], -1)
                                  - local[:, 1] * np.roll(local[:, 0], -1)))
        lo, hi, radius = influence_bounds(xy[i], vertices)
        radii[i] = radius
        old_extent[i] = [xy[i, 0] - 2*radius, xy[i, 1] - 2*radius,
                         xy[i, 0] + 2*radius, xy[i, 1] + 2*radius]
        new_extent[i] = [*lo, *hi]
        periodic_margin[i] = rectangle_margin((lo, hi), -length, 2*length)
        vertices_by_cell.append(region)
    # Certify that the nine-tile construction did not omit an influencing
    # generator outside that tiled domain. Require a positive numeric margin.
    slack = 1024 * np.finfo(float).eps * max(length)
    assert np.all(periodic_margin > slack), 'Nine images do not certify the full periodic tessellation'
    assert np.all(areas > 0) and np.isfinite(areas).all()
    relative = abs(areas / p['volume'].astype(float) - 1.)
    rows = []
    for halo in (0, 4, 8, 12, 16, 20, 24, 28, 32, 40):
        lower, upper = -halo*args.cellsize, length + halo*args.cellsize
        old_margin = np.min(np.column_stack((old_extent[:, :2] - lower,
                                             upper - old_extent[:, 2:])), axis=1)
        new_margin = np.min(np.column_stack((new_extent[:, :2] - lower,
                                             upper - new_extent[:, 2:])), axis=1)
        assert np.all(new_margin >= old_margin - slack)
        rows.append({'halo_cells': halo, 'old_unproven': int(np.sum(old_margin <= slack)),
                     'vertex_unproven': int(np.sum(new_margin <= slack)),
                     'old_min_margin': float(np.min(old_margin)),
                     'vertex_min_margin': float(np.min(new_margin))})
    lower28, upper28 = -28*args.cellsize, length + 28*args.cellsize
    old_margin28 = np.min(np.column_stack((old_extent[:, :2] - lower28,
                                           upper28 - old_extent[:, 2:])), axis=1)
    new_margin28 = np.min(np.column_stack((new_extent[:, :2] - lower28,
                                           upper28 - new_extent[:, 2:])), axis=1)
    details = []
    for i in np.flatnonzero(old_margin28 <= slack):
        details.append({'id': int(ids[i]), 'owner': xy[i].tolist(),
                        'radius_from_center_3_3': float(np.linalg.norm(xy[i]-3.)),
                        'circumradius': float(radii[i]),
                        'old_halo28_margin': float(old_margin28[i]),
                        'vertex_halo28_margin': float(new_margin28[i]),
                        'native_volume': float(p['volume'][i]),
                        'periodic_reference_area': float(areas[i]),
                        'area_relative_difference': float(relative[i]),
                        'vertices': vor.vertices[vertices_by_cell[i]].tolist()})
    checks = {'periodic_tiling_certified': bool(np.all(periodic_margin > slack)),
              'positive_reference_areas': bool(np.all(areas > 0)),
              'periodic_area_tiles_36': bool(abs(np.sum(areas)-36.) < 1e-6),
              'reference_matches_saved_volume_max_rel_below_1e_5': bool(np.max(relative) < 1e-5),
              'vertex_certificate_dominates_old_certificate': True}
    result = {'snapshot': str(args.snapshot), 'input_sha256': args.expected_sha256,
              'time': float(meta['t']), 'cells': n, 'W0_only': True,
              'reference': 'Qhull Voronoi of nine periodic image tiles, independently certified by vertex disks',
              'strict_numeric_margin': slack, 'periodic_reference_min_margin': float(np.min(periodic_margin)),
              'area_sum': float(np.sum(areas)),
              'area_native_relative_p99': float(np.quantile(relative, .99)),
              'area_native_relative_max': float(np.max(relative)),
              'halo_comparison': rows, 'old_halo28_unproven_cells': details,
              'reference_construction_seconds': built-started,
              'total_seconds': time.perf_counter()-started,
              'max_rss_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'checks': checks, 'passed': bool(all(checks.values())),
              'scope': 'Read-only accepted-state geometry reference; not native CUDA parity, not a failed RK4-stage reconstruction, not t>6 stability',
              'production_changed': False}
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'old_halo28_unproven_cells'}, indent=2))
    if not result['passed']:
        raise SystemExit('Geometry reference gate failed; result retained')


if __name__ == '__main__':
    main()
