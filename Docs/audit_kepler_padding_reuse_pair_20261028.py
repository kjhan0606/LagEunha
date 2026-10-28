"""Compare the isolated native padding-reuse runs and retain failing evidence."""
import argparse
import json
from pathlib import Path
import re

import numpy as np
from audit_kepler_polar_long_gate_20261027 import ordered

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/tess_padding_reuse_pair_20261028')
FIELDS = ('x', 'y', 'vx', 'vy', 'mass', 'volume', 'ie', 'den', 'pressure', 'w2')


def load(name):
    path = ROOT / name
    gate = json.loads((path / 'audit_endpoint_gate.json').read_text())
    state, time = ordered(Path(gate['endpoint']))
    text = (path / 'log').read_text(errors='replace')
    trace = re.findall(r'^Time=\s*([0-9.eE+-]+)\s+step=\s*(\d+)\s+dt=\s*([0-9.eE+-]+)', text, re.M)
    ledgers = re.findall(r'^\[RK4E\].*$', text, re.M)
    injection_zero = bool(ledgers)
    for key in ('n_ie_le0', 'floor_inj', 'floor_cum', 'de_inj', 'de_cum', 'sfl_inj', 'sfl_cum'):
        values = re.findall(rf'\b{key}=([0-9.eE+-]+)', '\n'.join(ledgers))
        injection_zero &= len(values) == len(ledgers) and all(np.isfinite(float(v)) and float(v) == 0 for v in values)
    wall = re.findall(r'\[RUN_TIME\] wall=([0-9.]+)', (path / 'runner.err').read_text())
    if len(wall) != 1:
        raise ValueError(f'Missing unique native wall time in {name}')
    timings = re.findall(r'\[GPU step \d+\].*?SoA=([0-9.]+)ms.*?kernel=([0-9.]+)ms.*?total=([0-9.]+)ms', text)
    metrics = {'time': time, 'wall_seconds': float(wall[0]), 'steps': len(trace),
               'all_native_physics_checks': all(gate['checks'].values()),
               'every_ledger_positive_and_injection_free': bool(injection_zero),
               'reuse_ranks': sorted(set(map(int, re.findall(r'\[TESS_PADDING_REUSE\] rank=(\d+).*reuse=1', text)))),
               'GPU_assignments': re.findall(r'^\[GPU\] rank .*$', text, re.M),
               'force_wrapper_mean_ms': float(np.mean([float(t[2]) for t in timings])) if timings else None,
               'force_sync_SoA_mean_ms': float(np.mean([float(t[0]) for t in timings])) if timings else None}
    return state, trace, metrics


def compare(a_name, b_name, exact):
    a, ta, ma = load(a_name)
    b, tb, mb = load(b_name)
    checks = {'both_physics_gates': ma['all_native_physics_checks'] and mb['all_native_physics_checks'],
              'all_ledgers_positive_without_injection': ma['every_ledger_positive_and_injection_free'] and mb['every_ledger_positive_and_injection_free'],
              'same_endpoint': abs(ma['time'] - mb['time']) < 1e-12,
              'printed_step_trace_identical': bool(ta) and ta == tb}
    differences = {}
    for f in FIELDS:
        delta = b[f].astype(float) - a[f].astype(float)
        scale = max(float(np.sqrt(np.mean(a[f].astype(float) ** 2))), 1e-300)
        row = {'bitwise_equal': bool(np.array_equal(a[f], b[f])),
               'max_absolute_difference': float(np.max(abs(delta))),
               'normalized_rms_difference': float(np.sqrt(np.mean(delta ** 2))) / scale}
        differences[f] = row
        if exact or f in ('mass', 'w2'):
            checks[f + '_exact'] = row['bitwise_equal']
        else:
            checks[f + '_normalized_rms_below_1e_6'] = row['normalized_rms_difference'] < 1e-6
            if f in ('x', 'y', 'vx', 'vy'):
                checks[f + '_max_abs'] = row['max_absolute_difference'] < (1e-6 if f in ('x', 'y') else 1e-5)
            else:
                row['max_cell_relative_difference'] = float(np.max(abs(delta) / np.maximum(abs(a[f].astype(float)), 1e-300)))
                checks[f + '_max_cell_relative_below_1e_3'] = row['max_cell_relative_difference'] < 1e-3
    return {'control': a_name, 'candidate': b_name, 'checks': checks,
            'passed': bool(all(checks.values())), 'field_differences': differences,
            'metrics': [ma, mb], 'wall_speedup_control_over_candidate': ma['wall_seconds'] / mb['wall_seconds']}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--qualifier', action='store_true')
    args = ap.parse_args()
    if args.qualifier:
        result = compare('qualifier_old', 'qualifier_new', exact=True)
        output = ROOT / 'qualifier.json'
    else:
        one = compare('np1_control', 'np1_reuse', exact=True)
        two = compare('np2_control', 'np2_reuse', exact=True)
        rank = compare('np1_reuse', 'np2_reuse', exact=False)
        banners = one['metrics'][1]['reuse_ranks'] == [0] and two['metrics'][1]['reuse_ranks'] == [0, 1]
        result = {'scope': 'Short operator parity and throughput check at the preserved late checkpoint; no t>6 claim',
                  'one_rank_pair': one, 'two_rank_pair': two, 'candidate_rank_comparison': rank,
                  'reuse_executed_on_all_ranks': banners,
                  'candidate_MPI_speedup_np1_over_np2': rank['wall_speedup_control_over_candidate'],
                  'passed': bool(one['passed'] and two['passed'] and rank['passed'] and banners),
                  'production_promoted': False}
        output = ROOT / 'comparison.json'
    with output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))
    if not result['passed']:
        raise SystemExit('Padding reuse gate failed; result retained')


if __name__ == '__main__':
    main()
