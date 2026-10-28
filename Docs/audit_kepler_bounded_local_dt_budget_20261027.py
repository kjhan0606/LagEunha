"""Estimate bounded q2 active-cell arithmetic from native K4 CFL fields."""

import argparse
import hashlib
import json
from pathlib import Path
import struct

import numpy as np


N = 55351
STRIDE = 576
ETA = .06
EPS = .01175


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def analyze(path, span):
    with path.open('rb') as stream:
        head = struct.unpack('<4i', stream.read(16))
        stream.seek(16 + N * STRIDE)
        t = struct.unpack('<d', stream.read(8))[0]
    if head != (-1, N, STRIDE, 1) or path.stat().st_size != 16 + N * STRIDE + 52:
        raise ValueError(f'Incomplete or wrong native snapshot: {path}')
    dtype = np.dtype({'names': ('u4if', 'x', 'y', 'mass', 'hydro_dt'),
                      'formats': ('<u8', '<f8', '<f8', '<f8', '<f4'),
                      'offsets': (8, 24, 32, 48, 128), 'itemsize': STRIDE})
    p = np.memmap(path, dtype=dtype, mode='r', offset=16, shape=(N,))
    ids = p['u4if'].astype(np.uint64) & np.uint64((1 << 56) - 1)
    if not np.array_equal(np.sort(ids), np.arange(N, dtype=np.uint64)):
        raise ValueError('Native IDs are incomplete')
    r = np.hypot(p['x'].astype(float) - 3, p['y'].astype(float) - 3)
    hydro = p['hydro_dt'].astype(float)
    if not np.all(np.isfinite(hydro) & (hydro > 0)):
        raise ValueError('Nonpositive or nonfinite native hydro CFL')
    orbit = ETA * (r*r + EPS*EPS)**.75
    local = np.minimum(hydro, orbit)
    global_step = float(np.min(local))
    levels = np.floor(np.log2(local / global_step) + 1e-12).astype(int)
    block_steps = global_step * np.exp2(levels)
    global_updates = N * int(np.ceil(span / global_step))
    block_updates = int(np.sum(np.ceil(span / block_steps), dtype=np.int64))
    continuous_ratio = float(N / np.sum(np.exp2(-levels)))
    masks = {'core': r < .5, 'disk': (r >= .5) & (r < 2), 'outer': r >= 2}
    regions = {}
    for name, selected in masks.items():
        regions[name] = {
            'cells': int(np.sum(selected)),
            'mass': float(np.sum(p['mass'][selected], dtype=np.longdouble)),
            'minimum_radius': float(np.min(r[selected])),
            'minimum_orbital_dt': float(np.min(orbit[selected])),
            'minimum_hydro_dt': float(np.min(hydro[selected])),
            'minimum_combined_dt': float(np.min(local[selected])),
            'median_combined_dt': float(np.median(local[selected])),
        }
    return {
        'snapshot': str(path.resolve()), 'snapshot_sha256': digest(path),
        'time': t, 'cells': N, 'span_for_frozen_budget': span,
        'global_step_from_native_fields': global_step,
        'minimum_hydro_dt': float(np.min(hydro)),
        'minimum_orbital_dt': float(np.min(orbit)),
        'cells_within_2x_global_step': int(np.sum(local < 2*global_step)),
        'regions': regions,
        'power2_levels': [{'level': int(i), 'cells': int(v)}
                          for i, v in enumerate(np.bincount(levels)) if v],
        'all_cell_updates_at_global_step': global_updates,
        'power2_block_active_cell_updates': block_updates,
        'power2_active_cell_arithmetic_ratio': global_updates / block_updates,
        'continuous_time_arithmetic_ratio': continuous_ratio,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('snapshots', nargs='+', type=Path)
    ap.add_argument('--span', type=float, default=6.5)
    ap.add_argument('--output', required=True, type=Path)
    args = ap.parse_args()
    if args.span <= 0 or args.output.exists():
        raise ValueError('Positive span and fresh output path required')
    result = {
        'scope': 'Read-only native K4 per-cell CFL plus actual orbital limiter. '
                 'Frozen-timestep active-cell arithmetic only; no moving-face '
                 'flux synchronization, GPU geometry cost, or runtime claim.',
        'eta': ETA, 'epsilon': EPS,
        'cases': [analyze(path, args.span) for path in args.snapshots],
    }
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
