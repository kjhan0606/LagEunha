"""Audit zero-time thermal64 promotion; this does not validate time evolution."""
import argparse
import hashlib
import json
from pathlib import Path
import socket
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_q1_checkpoint_geometry_20261029 import independent

TAG = 2026102603.0

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''): h.update(b)
    return h.hexdigest()

def payload(path, n):
    dt = np.dtype({'names': ['id', 'u', 'v', 'rhs', 'tag', 'anchor'],
                   'formats': ['<u8'] + ['<f8'] * 5,
                   'offsets': [8, 472, 480, 488, 496, 504], 'itemsize': 576})
    p = np.memmap(path, mode='r', dtype=dt, offset=16, shape=(n,))
    return p[np.argsort(p['id'] & ID_MASK)]

def main():
    a = argparse.ArgumentParser()
    a.add_argument('--source', type=Path, required=True)
    a.add_argument('--promoted', type=Path, required=True)
    a.add_argument('--output', type=Path, required=True)
    args = a.parse_args()
    if socket.gethostname().split('.')[0].lower() in ('grammar', 'syntax'):
        raise RuntimeError('Run geometry analysis on a compute host')
    if args.output.exists(): raise FileExistsError(args.output)
    p0, m0 = ordered_global(args.source)
    p1, m1 = ordered_global(args.promoted, p0['u4if'] & ID_MASK)
    q0, q1 = payload(args.source, len(p0)), payload(args.promoted, len(p1))
    v, geometry = independent(np.column_stack((p1['x'], p1['y'])))
    gamma_minus_one = float(np.float32(1.666667)) - 1.
    rel_v = q1['v']/v - 1.
    rho = p1['mass'].astype(float)/q1['v']
    pressure = gamma_minus_one*q1['u']/q1['v']
    checks = {
        '90352_same_ids_at_t10': len(p0) == 90352 and m0['t'] == m1['t'] == 10.,
        'initial_payload_untagged': bool(np.all(q0['tag'] != TAG)),
        'coordinates_velocities_mass_weight_energy_unchanged': all(np.array_equal(p0[k], p1[k]) for k in ('x', 'y', 'vx', 'vy', 'mass', 'w2', 'ie')),
        'U_exactly_promoted': bool(np.array_equal(q1['u'], p0['ie'].astype(float))),
        'all_tagged_rhs_zero_anchor_U': bool(np.all(q1['tag'] == TAG) and np.all(q1['rhs'] == 0) and np.array_equal(q1['anchor'], q1['u'])),
        'U_V_positive_finite': bool(np.all(np.isfinite(q1['u'])) and np.all(np.isfinite(q1['v'])) and np.all(q1['u'] > 0) and np.all(q1['v'] > 0)),
        'float_mirrors_exact': bool(np.array_equal(p1['ie'], q1['u'].astype(np.float32)) and np.array_equal(p1['volume'], q1['v'].astype(np.float32))),
        'density_eos_relative_1e-6': bool(np.max(np.abs(p1['den']/rho-1)) < 1e-6),
        'pressure_eos_relative_1e-6': bool(np.max(np.abs(p1['pressure']/pressure-1)) < 1e-6),
        'independent_volume_relative_1e-6': bool(np.max(np.abs(rel_v)) < 1e-6),
    }
    out = {'passed': all(checks.values()), 'checks': checks, 'sites': len(p1),
           'independent_geometry': geometry, 'volume_max_relative_error': float(np.max(np.abs(rel_v))),
           'volume_median_relative_error': float(np.median(np.abs(rel_v))),
           'volume_promotion_max_log_change': float(np.max(np.abs(np.log(q1['v']/p0['volume'])))),
           'source_sha256': sha(args.source), 'promoted_sha256': sha(args.promoted),
           'audit_sha256': sha(Path(__file__)),
           'scope': 'Zero-time state promotion only. Volume promotion changes stored geometric precision; it is not physical heating or validation of wall-event thermal64 evolution.'}
    with args.output.open('x') as f: json.dump(out, f, indent=2); f.write('\n')
    print(json.dumps(out, indent=2))
    if not out['passed']: raise SystemExit(1)

if __name__ == '__main__': main()
