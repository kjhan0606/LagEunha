"""Measure analytic Kepler face slip in one accepted fresh-weighted step."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_bounded_weighted_face_slip_t65_20261027 import diagram
from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    audit = json.loads(args.receipt.read_text())
    if len(audit['checks']) != 21 or not all(audit['checks'].values()):
        raise ValueError('Native weighted step is not accepted')
    parent, endpoint, initial = [Path(audit[key]) for key in
                                 ('parent', 'endpoint', 'initial')]
    if (digest(parent) != audit['parent_sha256'] or
            digest(endpoint) != audit['endpoint_sha256'] or
            digest(initial) != audit['initial_sha256']):
        raise ValueError('Native weighted snapshots changed')
    p0, _ = ordered(initial, N)
    a, ma = ordered(parent, N)
    b, mb = ordered(endpoint, N)
    dt = float(mb['t']) - float(ma['t'])
    if abs(dt - .001) > 1e-10 or not np.array_equal(a['mass'], b['mass']):
        raise ValueError('Not a single fixed-mass 0.001 step')
    mask = np.fromfile(args.receipt.parent / 'disk_active_mask.bin', dtype=np.uint8)
    if len(mask) != N or int(np.count_nonzero(mask)) != 49038:
        raise ValueError('Unexpected weighted active mask')
    active = np.flatnonzero(mask)
    lookup = np.full(N, -1, dtype=int)
    lookup[active] = np.arange(len(active))
    r0 = np.hypot(p0['x'].astype(float)-3., p0['y'].astype(float)-3.)
    xy = .5 * np.column_stack((a['x'].astype(float) + b['x'].astype(float)-6.,
                               a['y'].astype(float) + b['y'].astype(float)-6.))
    velocity = .5 * np.column_stack((a['vx'].astype(float)+b['vx'].astype(float),
                                     a['vy'].astype(float)+b['vy'].astype(float)))
    weight = .5*(a['w2'].astype(float)+b['w2'].astype(float))
    rate = (b['w2'].astype(float)-a['w2'].astype(float))/dt
    zero = np.zeros(N)
    modes = {
        'W0_same_generator_state': diagram(xy, velocity, r0, active, lookup,
                                            zero, zero),
        'power_frozen_weights': diagram(xy, velocity, r0, active, lookup,
                                         weight, zero),
        'power_moving_weights': diagram(xy, velocity, r0, active, lookup,
                                         weight, rate),
    }
    baseline = modes['W0_same_generator_state']['all_disk_disk']['moving']
    result = {
        'receipt': str(args.receipt),
        'receipt_sha256': digest(args.receipt),
        'start_sha256': audit['parent_sha256'],
        'end_sha256': audit['endpoint_sha256'],
        'start_time': float(ma['t']),
        'end_time': float(mb['t']),
        'dt': dt,
        'modes': modes,
        'moving_over_W0_all_disk_disk': {
            key: modes['power_moving_weights']['all_disk_disk']['moving'][key]/baseline[key]
            for key in ('integrated_absolute_slip', 'length_weighted_rms_slip')
        },
        'scope': 'Midpoint generator positions and velocities, secant weight rate, independently reconstructed power faces, analytic softened Kepler material velocity. This is a kinematic proxy, not native Riemann mass flux or stage-exact face transport.',
    }
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    for name, row in modes.items():
        values = row['all_disk_disk']['moving']
        print(name, values['faces'], values['integrated_absolute_slip'],
              values['length_weighted_rms_slip'])


if __name__ == '__main__':
    main()
