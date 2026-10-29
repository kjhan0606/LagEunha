"""Reconstruct instantaneous pressure work at the accepted early Q1 checkpoint."""
import hashlib
import json
from pathlib import Path
import socket
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from audit_kepler_bounded_pressure_force_t10_20261027 import finite_faces_touching_active
from audit_kepler_q1_t10_pressure_timescales_20261029 import describe

BASE = Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/ic_v1/khout.000001.dat')
END = Path('/gpfs/kjhan/LagForce/Kepler_Q1_OrbitalCap006_Ada_20261029/to0p4125/khout.005400.dat')
EXPECTED = 'df1eca2b0c341d0c4e61842c6e7ac8ec98ed88e24aea48ef63c1b19f6c57c924'
OUT = Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs/KEPLER_Q1_GRAVITY_CAP_2026-10-29/EARLY_WORK_RATE.json')


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as stream:
        for b in iter(lambda: stream.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def main():
    if socket.gethostname().split('.')[0].lower() in ('grammar', 'syntax'):
        raise RuntimeError('Run geometry analysis on a compute host')
    if sha(END) != EXPECTED:
        raise ValueError('Changed accepted endpoint')
    initial, im = ordered_global(BASE)
    p, meta = ordered_global(END, initial['u4if'] & ID_MASK)
    if im['t'] != 0 or meta['t'] != .4125 or not np.array_equal(initial['mass'], p['mass']) or np.any(p['w2'] != 0):
        raise ValueError('Unexpected input state')
    if np.any(p['ie'] <= 0) or np.any(~np.isfinite(p['ie'])):
        raise ValueError('Invalid thermal state')
    r0 = np.hypot(initial['x']-3., initial['y']-3.)
    xy = np.column_stack((p['x']-3., p['y']-3.))
    active = np.flatnonzero(r0 < 2.)
    flag = np.zeros(len(p), dtype=bool); flag[active] = True
    lookup = np.full(len(p), -1, dtype=int); lookup[active] = np.arange(len(active))
    area, geom, tri, centers = power_areas(xy, np.zeros(len(p)), active, lookup, return_geometry=True)
    areaerr = float(np.max(np.abs(area / p['volume'][active].astype(float)-1)))
    if geom['hidden_active'] or np.any(area <= 0) or areaerr >= 1e-5:
        raise ValueError('Independent volume mismatch')
    i, j, center, length = finite_faces_touching_active(tri, centers, flag)
    delta = xy[j]-xy[i]; dist = np.linalg.norm(delta, axis=1)
    if np.any(dist <= 0):
        raise ValueError('Coincident generators')
    normal = delta/dist[:, None]; offset = center-.5*(xy[i]+xy[j])
    vel = np.column_stack((p['vx'], p['vy']))
    speed = .5*np.einsum('ij,ij->i', vel[i]+vel[j], normal) - np.einsum('ij,ij->i', vel[j]-vel[i], offset)/dist
    flux = length*speed
    dv = np.zeros(len(p)); ii = flag[i]; jj = flag[j]
    np.add.at(dv, i[ii], flux[ii]); np.add.at(dv, j[jj], -flux[jj])
    udot = -p['pressure'].astype(float)*dv
    radius = np.hypot(xy[:, 0], xy[:, 1]); omega = (radius*radius+.01175**2)**(-.75)
    ratio = np.abs(udot)/(p['ie'].astype(float)*omega)
    mass = p['mass'].astype(float)
    regions = []
    for lo, hi in [(0., .1), (.1, .5), (.5, .75), (.75, 1.), (1., 2.), (.5, 2.)]:
        mask = (r0 >= lo) & (r0 < hi)
        regions.append({'rlo': lo, 'rhi': hi, 'cells': int(mask.sum()),
                        'thermal_fraction_per_dynamical_time': describe(ratio[mask], mass[mask])})
    result = {'time': meta['t'], 'cells': len(active), 'faces': len(i), 'max_volume_relative_error': areaerr,
              'regions': regions, 'initial_sha256': sha(BASE), 'endpoint_sha256': sha(END), 'script_sha256': sha(Path(__file__)),
              'scope': 'Instantaneous independent geometric reconstruction at an early screened checkpoint, not a native RHS comparison at this time, entropy production, a temporal convergence test, or a measured speedup. The same method was compared against eight native RHS samples at t10.'}
    with OUT.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
