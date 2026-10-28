"""Audit both preserved cap48/cap16 runs, retaining evidence on failure."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np

from audit_kepler_polar_long_gate_20261027 import ordered


Q1 = Path('/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1')
CONTROL = Q1 / 'dt_areaavg_audit_t150_to1515_a100_20261028'
CANDIDATE = Q1 / 'kernA_cap16_probe_t150_to1515_a100_20261028'
FIELDS = ('x', 'y', 'vx', 'vy', 'mass', 'volume', 'ie', 'den', 'pressure', 'w2')
TARGET = 1.515


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def inspect(run):
    log = (run / 'log').read_text(errors='replace') if (run / 'log').exists() else ''
    trace = re.findall(r'^Time=\s*([0-9.eE+-]+)\s+step=\s*(\d+)\s+dt=\s*([0-9.eE+-]+)', log, re.M)
    faces = re.findall(r'^\[GPU tess\].*?faces=(\d+).*?time=([0-9.eE+-]+)\s*ms', log, re.M)
    ledgers = re.findall(r'^\[RK4E\].*$', log, re.M)
    keys = ('n_ie_le0', 'floor_inj', 'floor_cum', 'de_inj', 'de_cum', 'sfl_inj', 'sfl_cum')
    maxima = {}
    for key in keys:
        matches = [re.search(rf'\b{key}=([0-9.eE+-]+)', line) for line in ledgers]
        values = [float(m[1]) for m in matches if m]
        maxima[key] = max(abs(v) for v in values) if values else None
    all_ledger_fields_present = bool(ledgers) and all(
        len(re.findall(rf'\b{key}=([0-9.eE+-]+)', '\n'.join(ledgers))) == len(ledgers)
        for key in keys)
    all_zero = all_ledger_fields_present and all(
        value is not None and np.isfinite(value) and value == 0.0
        for value in maxima.values())
    wall = re.findall(r'^\[PAIR_RUN_TIME\] wall_seconds=(\d+)$', log, re.M)
    matches = []
    errors = []
    for path in run.glob('khout.*.dat'):
        try:
            state, time = ordered(path)
            if abs(time - TARGET) < 1e-10:
                matches.append((path, state, time))
        except Exception as error:
            errors.append({'path': str(path), 'error': str(error)})
    state = matches[0][1] if len(matches) == 1 else None
    path = matches[0][0] if len(matches) == 1 else None
    result = {
        'run': str(run), 'exit_zero': bool(re.search(r'^EXIT:0$', log, re.M)),
        'accepted_steps': len(trace), 'last_time': float(trace[-1][0]) if trace else None,
        'ledger_count': len(ledgers), 'all_ledger_fields_present': all_ledger_fields_present,
        'ledger_max_absolute_values': maxima, 'all_ledger_injection_and_nonpositive_counts_zero': bool(all_zero),
        'capacity_abort_seen': '[GPU_CLIP_CAPACITY]' in log,
        'tessellation_invalid_seen': '[GPU_TESS_INVALID]' in log,
        'endpoint_count': len(matches), 'endpoint': str(path) if path else None,
        'endpoint_sha256': digest(path) if path else None,
        'wall_seconds': int(wall[0]) if len(wall) == 1 else None,
        'mean_tessellation_ms_excluding_first_four_calls': float(np.mean([float(ms) for _, ms in faces[4:]])) if len(faces) > 4 else None,
        'snapshot_read_errors': errors,
    }
    return result, state, trace, [int(count) for count, _ in faces]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    checks = {}
    inspected = []
    for label, root in (('cap48', CONTROL), ('cap16', CANDIDATE)):
        row, state, trace, faces = inspect(root)
        inspected.append((row, state, trace, faces))
        checks[label + '_exited_zero'] = row['exit_zero']
        checks[label + '_all_logged_injection_and_nonpositive_counts_zero'] = row['all_ledger_injection_and_nonpositive_counts_zero']
        checks[label + '_no_capacity_or_geometry_abort'] = not (row['capacity_abort_seen'] or row['tessellation_invalid_seen'])
        checks[label + '_unique_exact_endpoint'] = row['endpoint_count'] == 1
    checks['accepted_time_step_trace_identical'] = bool(inspected[0][2]) and inspected[0][2] == inspected[1][2]
    checks['all_tessellation_face_counts_identical'] = bool(inspected[0][3]) and inspected[0][3] == inspected[1][3]
    field_results = {}
    initial, time0 = ordered(Q1 / 'khout.000000.dat')
    checks['initial_time_zero'] = time0 == 0
    for label, (row, state, _, _) in zip(('cap48', 'cap16'), inspected):
        checks[label + '_fixed_mass'] = state is not None and bool(np.array_equal(state['mass'], initial['mass']))
        checks[label + '_W0'] = state is not None and bool(np.all(state['w2'] == 0))
        checks[label + '_finite'] = state is not None and all(bool(np.all(np.isfinite(state[f]))) for f in FIELDS)
        checks[label + '_positive_U_and_volume'] = state is not None and bool(np.all(state['ie'] > 0) and np.all(state['volume'] > 0))
    a, b = inspected[0][1], inspected[1][1]
    if a is not None and b is not None:
        for field in FIELDS:
            delta = b[field].astype(float) - a[field].astype(float)
            field_results[field] = {
                'bitwise_equal': bool(np.array_equal(a[field], b[field])),
                'different_cells': int(np.count_nonzero(a[field] != b[field])),
                'max_absolute_difference': float(np.max(abs(delta))),
                'rms_difference': float(np.sqrt(np.mean(delta * delta))),
            }
    checks['all_ten_endpoint_fields_exact'] = len(field_results) == len(FIELDS) and all(row['bitwise_equal'] for row in field_results.values())
    wall = [item[0]['wall_seconds'] for item in inspected]
    tess = [item[0]['mean_tessellation_ms_excluding_first_four_calls'] for item in inspected]
    result = {
        'scope': 'short exact-output capacity optimization audit; production requires full-capacity fallback and long evolution gates',
        'endpoint_time': TARGET, 'checks': checks, 'passed': bool(all(checks.values())),
        'runs': [row[0] for row in inspected], 'field_comparison': field_results,
        'wall_speedup_cap48_over_cap16': wall[0] / wall[1] if all(v is not None and v > 0 for v in wall) else None,
        'tessellation_speedup_cap48_over_cap16': tess[0] / tess[1] if all(v is not None and v > 0 for v in tess) else None,
    }
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))
    if not result['passed']:
        raise SystemExit('Extended cap16 comparison failed; evidence retained in output JSON')


if __name__ == '__main__':
    main()
