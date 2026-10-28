"""Audit the paper-profile t=.001 output using generator IDs to reconstruct t=0.

The first native output is at t=.001. Initial masses and Cartesian positions
are exactly recoverable from low 56 index bits in the packed particle ID.
"""

from pathlib import Path
import hashlib
import json
import sys

import numpy as np

SOURCE = Path('/gpfs/kjhan/LagForce/Kepler_HopkinsIC384_H200_20261026/khout.000009.dat')
ROOT = Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
sys.path.insert(0, str(ROOT))
from audit_snapshot import read_snapshot


def main():
    p, meta = read_snapshot(SOURCE)
    if (meta['nx'], meta['ny'], meta['Lx'], meta['Ly']) != (384, 384, 6., 6.):
        raise ValueError('Unexpected native IC grid')
    index = (p['u4if'].astype(np.uint64) & np.uint64((1 << 56) - 1)).astype(np.int64)
    if not np.array_equal(np.sort(index), np.arange(384 * 384)):
        raise ValueError('Packed particle IDs do not cover the initial grid')
    h = 6. / 384.
    x0 = (index % 384 + .5) * h
    y0 = (index // 384 + .5) * h
    r0 = np.hypot(x0 - 3., y0 - 3.)
    shape = np.where(r0 < .5, (r0 / .5) ** 3,
                     np.where(r0 <= 2., 1., (1. + (r0 - 2.) / .1) ** -3))
    rho0 = .01 + shape
    expected_mass = rho0 * h * h
    mass_relative_error = abs(p['mass'].astype(float) / expected_mass - 1.)
    angle = float(meta['t']) * (r0 * r0 + .01175 ** 2) ** (-.75)
    dx0, dy0 = x0 - 3., y0 - 3.
    xb = 3. + dx0 * np.cos(angle) - dy0 * np.sin(angle)
    yb = 3. + dx0 * np.sin(angle) + dy0 * np.cos(angle)
    position_error = np.hypot(p['x'] - xb, p['y'] - yb)
    ring = (r0 >= .5) & (r0 <= 2.)
    finite = all(np.isfinite(p[key]).all() for key in
                 ('x', 'y', 'vx', 'vy', 'mass', 'volume', 'ie', 'den', 'pressure', 'w2'))
    pressure_ring = p['pressure'][ring].astype(float)
    density_error_ring = abs(p['den'][ring].astype(float) / rho0[ring] - 1.)
    checks = dict(finite=finite,
                  exact_initial_mass=bool(mass_relative_error.max() < 1e-12),
                  positive_mass=bool(np.all(p['mass'] > 0)),
                  positive_volume=bool(np.all(p['volume'] > 0)),
                  positive_U=bool(np.all(p['ie'] > 0)),
                  W0=bool(np.all(p['w2'] == 0)),
                  exact_orbit_ring_to_1e_minus_10=bool(position_error[ring].max() < 1e-10),
                  ring_pressure_within_1_percent=bool(np.max(abs(pressure_ring / 1e-6 - 1.)) < .01))
    digest = hashlib.sha256()
    with SOURCE.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(block)
    result = dict(source=str(SOURCE), sha256=digest.hexdigest(), native_time=float(meta['t']),
                  cells=len(p), effective_resolution_across_width4=256,
                  initial_mass_relative_error_max=float(mass_relative_error.max()),
                  ring_cells=int(ring.sum()),
                  ring_position_error_ballistic_p50_p99_max=np.quantile(position_error[ring], [.5, .99, 1.]).tolist(),
                  ring_density_error_initial_p50_p99_max=np.quantile(density_error_ring, [.5, .99, 1.]).tolist(),
                  ring_pressure_p01_p50_p99=np.quantile(pressure_ring, [.01, .5, .99]).tolist(),
                  all_pressure_min_max=[float(np.min(p['pressure'])), float(np.max(p['pressure']))],
                  checks=checks,
                  limitation='First saved native output is t=.001, not t=0; reconstructed initial mass and exact orbit are verified by IDs. This is a preflight only, not long-time stability.')
    output = Path(__file__).with_name('HOPKINS_PREFLIGHT_RECONSTRUCT_T001_2026-10-26.json')
    with output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))
    if not all(checks.values()):
        raise SystemExit('Preflight audit gate failed')


if __name__ == '__main__':
    main()
