"""Check local total-PdV closure of two accepted native weighted RK4 steps."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
GAMMA_MINUS_ONE = float(np.float32(1.666667)) - 1.0


def weighted_mean(values, mass):
    return float(np.sum(values * mass, dtype=np.longdouble) /
                 np.sum(mass, dtype=np.longdouble))


def measure(receipt):
    audit = json.loads(receipt.read_text())
    if len(audit['checks']) != 21 or not all(audit['checks'].values()):
        raise ValueError(f'Native endpoint audit failed: {receipt}')
    start_path = Path(audit['parent'])
    end_path = Path(audit['endpoint'])
    if (digest(start_path) != audit['parent_sha256'] or
            digest(end_path) != audit['endpoint_sha256']):
        raise ValueError('Native snapshots changed')
    start, meta0 = ordered(start_path, N)
    end, meta1 = ordered(end_path, N)
    dt = float(meta1['t']) - float(meta0['t'])
    if (abs(dt - .001) > 1e-10 or
            not np.array_equal(start['mass'], end['mass'])):
        raise ValueError('Expected one fixed-mass 0.001 step')
    mask_path = receipt.parent / 'disk_active_mask.bin'
    active = np.fromfile(mask_path, dtype=np.uint8)
    if len(active) != N or int(np.count_nonzero(active)) != 49038:
        raise ValueError('Unexpected active mask')
    active = active.astype(bool)
    mass = start['mass'][active].astype(float)
    u0 = start['ie'][active].astype(float)
    u1 = end['ie'][active].astype(float)
    v0 = start['volume'][active].astype(float)
    v1 = end['volume'][active].astype(float)
    p0 = start['pressure'][active].astype(float)
    p1 = end['pressure'][active].astype(float)
    if not np.all((u0 > 0) & (u1 > 0) & (v0 > 0) & (v1 > 0)):
        raise ValueError('Nonpositive active thermodynamic state')
    pdv_residual = ((u1 - u0) + .5 * (p0 + p1) * (v1 - v0)) / u0
    entropy_log = np.log(u1 / u0) + GAMMA_MINUS_ONE * np.log(v1 / v0)
    return {
        'receipt': str(receipt),
        'receipt_sha256': digest(receipt),
        'start_sha256': audit['parent_sha256'],
        'end_sha256': audit['endpoint_sha256'],
        'start_time': float(meta0['t']),
        'end_time': float(meta1['t']),
        'dt': dt,
        'active_cells': int(active.sum()),
        'total_pdv_residual_over_initial_internal_energy': {
            'abs_p50': float(np.quantile(np.abs(pdv_residual), .5)),
            'abs_p99': float(np.quantile(np.abs(pdv_residual), .99)),
            'abs_max': float(np.max(np.abs(pdv_residual))),
            'abs_mass_mean': weighted_mean(np.abs(pdv_residual), mass),
            'signed_mass_mean': weighted_mean(pdv_residual, mass),
        },
        'entropy_log_step': {
            'abs_p50': float(np.quantile(np.abs(entropy_log), .5)),
            'abs_p99': float(np.quantile(np.abs(entropy_log), .99)),
            'abs_max': float(np.max(np.abs(entropy_log))),
            'abs_mass_mean': weighted_mean(np.abs(entropy_log), mass),
            'signed_mass_mean': weighted_mean(entropy_log, mass),
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--early', type=Path, required=True)
    parser.add_argument('--late', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = {
        'early': measure(args.early),
        'late': measure(args.late),
        'scope': 'Endpoint trapezoidal local internal-energy/PdV and entropy consistency. Does not resolve four-stage pairwise transport or physical material-face motion.',
    }
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
