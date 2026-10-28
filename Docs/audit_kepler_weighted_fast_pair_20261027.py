"""Compare independently audited weighted native runs from one restart."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest

N = 55351
FIELDS = ('x', 'y', 'vx', 'vy', 'mass', 'volume', 'ie', 'den', 'pressure', 'w2')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('control', type=Path)
    ap.add_argument('candidate', type=Path)
    ap.add_argument('--control-binary-sha256', required=True)
    ap.add_argument('--candidate-binary-sha256', required=True)
    ap.add_argument('--parent-time', type=float, required=True)
    ap.add_argument('--time', type=float, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    a = json.loads(args.control.read_text())
    b = json.loads(args.candidate.read_text())
    if (len(a['checks']) != 21 or len(b['checks']) != 21 or
            not all(a['checks'].values()) or not all(b['checks'].values()) or
            a['binary_sha256'] != args.control_binary_sha256 or
            b['binary_sha256'] != args.candidate_binary_sha256 or
            a['initial_sha256'] != b['initial_sha256'] or
            a['parent_sha256'] != b['parent_sha256'] or
            any(abs(r['parent_time']-args.parent_time) > 1e-12 or
                abs(r['time']-args.time) > 1e-12 for r in (a, b))):
        raise ValueError('Native physics, lineage, binary or time check failed')
    arrays = {}
    for label, receipt in (('control', a), ('candidate', b)):
        for name in ('initial', 'parent', 'endpoint'):
            path = Path(receipt[name])
            if digest(path) != receipt[name+'_sha256']:
                raise ValueError(f'Changed {label} {name}')
            arrays[label, name], _ = ordered(path, N)
    comparisons = {}
    for phase in ('initial', 'parent', 'endpoint'):
        comparisons[phase] = {}
        for field in FIELDS:
            x, y = arrays['control', phase][field], arrays['candidate', phase][field]
            comparisons[phase][field] = {
                'bitwise_equal': bool(np.array_equal(x, y)),
                'max_abs_difference': float(np.max(np.abs(x.astype(float)-y.astype(float))))}
    checks = {'both_physics_gates': True,
              **{f'{phase}_fields_bitwise': all(v['bitwise_equal']
                 for v in comparisons[phase].values()) for phase in comparisons}}
    result = {'control_receipt': str(args.control.resolve()),
              'candidate_receipt': str(args.candidate.resolve()),
              'control_receipt_sha256': digest(args.control),
              'candidate_receipt_sha256': digest(args.candidate),
              'parent_time': args.parent_time, 'time': args.time,
              'control_binary_sha256': args.control_binary_sha256,
              'candidate_binary_sha256': args.candidate_binary_sha256,
              'fields': comparisons, 'checks': checks,
              'scope': 'Short same-input native weighted comparison. This does not by itself '
                       'establish long-time equivalence, physical closure or production speed.'}
    with args.output.open('x') as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write('\n')
    print(json.dumps({'checks': checks}, indent=2))
    if not all(checks.values()):
        raise SystemExit('Weighted physical fields differ')


if __name__ == '__main__':
    main()
