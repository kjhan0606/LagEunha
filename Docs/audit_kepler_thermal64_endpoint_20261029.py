"""Audit authoritative thermal64 state and wall budget for a saved native interval."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_q1_thermal64_promotion_20261029 import payload, TAG
from audit_kepler_wall_angular_budget_20261029 import budget


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def energy(p, u, eps):
    m = p['mass'].astype(np.float64)
    x, y = p['x'] - 3., p['y'] - 3.
    return float(np.sum(u + .5*m*(p['vx']**2 + p['vy']**2)
                        - m/np.sqrt(x*x + y*y + eps*eps), dtype=np.longdouble))


def run(initial, endpoint, logs, eps):
    p0, m0 = ordered_global(initial)
    p, meta = ordered_global(endpoint, p0['u4if'] & ID_MASK)
    q0, q = payload(initial, len(p0)), payload(endpoint, len(p))
    start_tagged = np.all(q0['tag'] == TAG)
    if not start_tagged and np.any(q0['tag'] == TAG):
        raise ValueError('Mixed initial thermal64 tags')
    u0 = q0['u'] if start_tagged else p0['ie'].astype(np.float64)
    if np.any(~np.isfinite(u0)) or np.any(u0 <= 0):
        raise ValueError('Invalid initial thermal energy')
    if start_tagged and (np.any(~np.isfinite(q0['v'])) or np.any(q0['v'] <= 0)):
        raise ValueError('Invalid initial thermal64 volume')
    if np.any(~np.isfinite(q['u'])) or np.any(~np.isfinite(q['v'])) or np.any(q['u'] <= 0) or np.any(q['v'] <= 0):
        raise ValueError('Invalid thermal64 U or V')
    gamma_minus_one = float(np.float32(1.666667)) - 1.
    rho = p['mass'].astype(np.float64)/q['v']
    pressure = gamma_minus_one*q['u']/q['v']
    rhoerr = float(np.max(np.abs(p['den']/rho - 1.)))
    perr = float(np.max(np.abs(p['pressure']/pressure - 1.)))
    angular = budget(initial, endpoint, logs)
    e0, e1 = energy(p0, u0, eps), energy(p, q['u'], eps)
    checks = {
        'time_advanced': bool(meta['t'] > m0['t']),
        'fixed_mass': bool(np.array_equal(p0['mass'], p['mass'])),
        'zero_weight': bool(np.all(p0['w2'] == 0) and np.all(p['w2'] == 0)),
        'thermal_tag': bool(np.all(q['tag'] == TAG)),
        'positive_finite_U_V': True,
        'float_U_V_mirrors': bool(np.array_equal(p['ie'], q['u'].astype(np.float32)) and np.array_equal(p['volume'], q['v'].astype(np.float32))),
        'density_eos_relative_1e-6': rhoerr < 1e-6,
        'pressure_eos_relative_1e-6': perr < 1e-6,
        'finite_total_energy': bool(np.isfinite(e0) and np.isfinite(e1) and e0 != 0),
        'wall_angular_budget': bool(angular['passed']),
    }
    return {'passed': all(checks.values()), 'checks': checks,
            'time_initial': m0['t'], 'time_final': meta['t'],
            'initial_thermal64_tagged': bool(start_tagged),
            'energy_initial': e0, 'energy_final': e1,
            'energy_relative_change': (e1-e0)/abs(e0),
            'density_max_relative_eos_error': rhoerr,
            'pressure_max_relative_eos_error': perr,
            'angular_budget': angular,
            'audit_sha256': sha(Path(__file__)),
            'scope': 'Endpoint EOS and full-interval angular accounting in bounded W0 fixed-mass run. Energy change is reported, not given an a priori pass threshold. No shock, MPI, core-accuracy, boundary-convergence or long-time certification.'}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--initial', required=True, type=Path)
    ap.add_argument('--endpoint', required=True, type=Path)
    ap.add_argument('--logs', required=True, nargs='+', type=Path)
    ap.add_argument('--output', required=True, type=Path)
    ap.add_argument('--eps', type=float, default=.01175)
    args = ap.parse_args()
    result = run(args.initial, args.endpoint, args.logs, args.eps)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))
    if not result['passed']:
        raise SystemExit(1)
