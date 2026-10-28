"""Audit exact t=0 and t=.001 native Hopkins-profile polar IC snapshots."""

from pathlib import Path
import hashlib
import json
import struct
import sys

import numpy as np

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
sys.path.insert(0, str(ROOT))
from audit_snapshot import read_snapshot


def read_glass(path):
    with path.open('rb') as stream:
        n, lx, ly = struct.unpack('<idd', stream.read(20))
        if (n, lx, ly) != (384*384, 6., 6.):
            raise ValueError('Unexpected glass header')
        x = np.fromfile(stream, dtype='<f8', count=n)
        y = np.fromfile(stream, dtype='<f8', count=n)
    return np.column_stack((x, y))


def low_index(p):
    return (p['u4if'].astype(np.uint64) & np.uint64((1 << 56) - 1)).astype(np.int64)


def main():
    run = Path.cwd()
    initial_file = run / 'khout.000000.dat'
    evolved = sorted(run.glob('khout.*.dat'))
    if len(evolved) != 2 or evolved[0] != initial_file:
        raise ValueError('Expected only exact initial and one endpoint snapshot')
    p0, meta0 = read_snapshot(initial_file)
    p1, meta1 = read_snapshot(evolved[1])
    if meta0['t'] != 0. or abs(meta1['t'] - .001) > 1e-12:
        raise ValueError('Unexpected snapshot times')
    i0, i1 = low_index(p0), low_index(p1)
    expected_ids = np.arange(384*384)
    if not (np.array_equal(np.sort(i0), expected_ids) and
            np.array_equal(np.sort(i1), expected_ids)):
        raise ValueError('Snapshot IDs do not cover the grid')
    p0, p1 = p0[np.argsort(i0)], p1[np.argsort(i1)]
    xy = read_glass(run / 'glass.bin')
    position_initial_max = float(np.max(np.hypot(p0['x']-xy[:, 0], p0['y']-xy[:, 1])))
    r0 = np.hypot(xy[:, 0]-3., xy[:, 1]-3.)
    profile = np.where(r0 < .5, (r0/.5)**3,
                       np.where(r0 <= 2., 1., (1. + (r0-2.)/.1)**-3))
    rho_target = .01+profile
    mass_expected = rho_target*p0['volume'].astype(float)
    mass_initial_relative_error_max = float(np.max(abs(p0['mass']/mass_expected-1.)))
    disk = (r0 >= .5) & (r0 < 2.)
    angle = float(meta1['t'])*(r0*r0+.01175**2)**(-.75)
    dx0, dy0 = xy[:, 0]-3., xy[:, 1]-3.
    xb = 3.+dx0*np.cos(angle)-dy0*np.sin(angle)
    yb = 3.+dx0*np.sin(angle)+dy0*np.cos(angle)
    orbit_error = np.hypot(p1['x']-xb, p1['y']-yb)
    density_change = abs(p1['den'][disk].astype(float)/p0['den'][disk].astype(float)-1.)
    finite = all(np.isfinite(p[key]).all() for p in (p0, p1) for key in
                 ('x','y','vx','vy','mass','volume','ie','den','pressure','w2'))
    checks = dict(finite=bool(finite),
                  initial_position_matches_glass=bool(position_initial_max < 1e-12),
                  initial_mass_matches_target_density_times_cell_volume_at_snapshot_precision=bool(mass_initial_relative_error_max < 1e-6),
                  mass_exact_across_step=bool(np.array_equal(p0['mass'], p1['mass'])),
                  volume_positive=bool(np.all(p0['volume']>0) and np.all(p1['volume']>0)),
                  U_positive=bool(np.all(p0['ie']>0) and np.all(p1['ie']>0)),
                  W0=bool(np.all(p0['w2']==0) and np.all(p1['w2']==0)),
                  disk_orbit_error_below_1e_minus_6=bool(np.max(orbit_error[disk]) < 1e-6),
                  disk_pressure_near_1e_minus_6_initial=bool(np.max(abs(p0['pressure'][disk]/1e-6-1.)) < 1e-6))
    digest = {}
    for path in (run/'glass.bin', initial_file, evolved[1], run/'eunha2.exe'):
        h = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(4*1024*1024), b''):
                h.update(block)
        digest[path.name] = h.hexdigest()
    result = dict(initial=str(initial_file), evolved=str(evolved[1]),
                  t0=meta0['t'], t1=meta1['t'], n=len(p0),
                  sha256=digest,
                  initial_position_error_max=position_initial_max,
                  initial_mass_relative_error_max=mass_initial_relative_error_max,
                  volume_storage_dtype=str(p0['volume'].dtype),
                  disk_cells=int(np.sum(disk)),
                  disk_orbit_error_p50_p99_max=np.quantile(orbit_error[disk],[.5,.99,1.]).tolist(),
                  disk_density_change_p50_p95_max=np.quantile(density_change,[.5,.95,1.]).tolist(),
                  initial_volume_min=float(np.min(p0['volume'])),
                  evolved_volume_min=float(np.min(p1['volume'])),
                  checks=checks,
                  scope='Native full RK4+CUDA preflight with opt-in polar generator positions and one-time exact initial cell-mass assignment; no long-time stability claim')
    with (run/'audit_t001.json').open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))
    if not all(checks.values()):
        raise SystemExit('Polar IC preflight gate failed')


if __name__ == '__main__':
    main()
