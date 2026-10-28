"""Compare short native bounded Kepler replays and partition RK4 wall time."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import statistics

from prepare_kepler_sparse_area_restart_20261027 import snapshot_time
from audit_kepler_sparse_area_t0825_20261027 import ordered
import numpy as np


COUNT = 55351
TARGET = 1.6785
BASE_SHA = '3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e'
PROFILE_SHA = 'e10cca50e3e158baeca28fea142b9f061ddafd17aade3756a8bc42242ae9761c'
FIELDS = ('wall', 'tracked', 'prelude', 'after_update', 'after_force',
          'after_post', 'gaps', 'closure')


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def endpoint(directory):
    candidates = [p for p in directory.glob('khout.*.dat')
                  if abs(snapshot_time(p, COUNT) - TARGET) < 1e-10]
    if len(candidates) != 1:
        raise ValueError(f'Expected one complete exact target in {directory}')
    return candidates[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, default=Path('audit_profile_corrected_v2.json'))
    args = ap.parse_args()
    root = Path.cwd()
    base, prof = root / 'base', root / 'profile'
    b_end, p_end = endpoint(base), endpoint(prof)
    b_state, b_meta = ordered(b_end, COUNT)
    p_state, p_meta = ordered(p_end, COUNT)
    mismatched_fields = [name for name in b_state.dtype.names
                         if not np.array_equal(b_state[name], p_state[name])]
    b_bytes = np.fromfile(b_end, dtype=np.uint8)
    p_bytes = np.fromfile(p_end, dtype=np.uint8)
    if b_bytes.size != 16 + COUNT * 576 + 52 or p_bytes.size != b_bytes.size:
        raise ValueError('Unexpected native snapshot byte length')
    b_records = b_bytes[16:-52].reshape(COUNT, 576)
    p_records = p_bytes[16:-52].reshape(COUNT, 576)
    pointer_rows_different = int(np.count_nonzero(
        np.any(b_records[:, :8] != p_records[:, :8], axis=1)))
    nonpointer_equal = bool(np.array_equal(b_records[:, 8:], p_records[:, 8:]))
    header_footer_equal = bool(np.array_equal(b_bytes[:16], p_bytes[:16])
                               and np.array_equal(b_bytes[-52:], p_bytes[-52:]))
    blog, plog = (base / 'log').read_text(errors='replace'), (prof / 'log').read_text(errors='replace')
    entries = []
    for line in plog.splitlines():
        if '[RK4 gap]' not in line:
            continue
        row = {'step': int(re.search(r'\bstep=(\d+)', line).group(1))}
        for name in FIELDS:
            match = re.search(rf'\b{name}=([0-9.eE+-]+)s?', line)
            if not match:
                raise ValueError(f'Missing {name} from {line}')
            row[name] = float(match.group(1))
        entries.append(row)
    base_steps = re.findall(r'\[CPU step (\d+)\]', blog)
    profile_steps = re.findall(r'\[CPU step (\d+)\]', plog)
    result = {
        'restart_time': 1.651374332206546,
        'target_time': TARGET,
        'base_endpoint': str(b_end), 'profile_endpoint': str(p_end),
        'base_endpoint_sha256': digest(b_end),
        'profile_endpoint_sha256': digest(p_end),
        'parsed_snapshot_field_count': len(b_state.dtype.names),
        'mismatched_parsed_fields': mismatched_fields,
        'metadata_equal': b_meta == p_meta,
        'serialized_bytes_identical': digest(b_end) == digest(p_end),
        'serialized_next_pointer_rows_different': pointer_rows_different,
        'all_nonpointer_record_bytes_identical': nonpointer_equal,
        'header_footer_bytes_identical': header_footer_equal,
        'base_binary_sha256': digest(base / 'eunha2.exe'),
        'profile_binary_sha256': digest(prof / 'eunha2.exe'),
        'base_steps': len(base_steps), 'profile_steps': len(profile_steps),
        'gap_rows': len(entries),
        'median_seconds': {name: statistics.median(row[name] for row in entries)
                           for name in FIELDS} if entries else {},
        'maximum_abs_closure_seconds': max((abs(row['closure']) for row in entries), default=None),
        'checks': {
            'accepted_base_binary': digest(base / 'eunha2.exe') == BASE_SHA,
            'isolated_profile_binary': digest(prof / 'eunha2.exe') == PROFILE_SHA,
            'same_accepted_restart': digest(base / 'khout.000001.dat') == digest(prof / 'khout.000001.dat'),
            'all_parsed_fields_bitwise_equal': not mismatched_fields,
            'snapshot_metadata_equal': b_meta == p_meta,
            'all_nonpointer_record_bytes_identical': nonpointer_equal,
            'header_footer_bytes_identical': header_footer_equal,
            'matched_steps': len(base_steps) == len(profile_steps) == len(entries) >= 10,
            'timing_closure': len(entries) >= 10 and all(abs(row['closure']) < 1e-3 for row in entries),
            'gap_partition': len(entries) >= 10 and all(abs(row['gaps'] - sum(row[n] for n in
                ('prelude', 'after_update', 'after_force', 'after_post'))) < 5e-6 for row in entries),
            'run_exit_zero': 'EXIT:0' in blog and 'EXIT:0' in plog,
        },
        'scope': '14-step one-rank A10 comparison; each 576-byte record matches except offset 0:8, the serialized linked-list next pointer; times partition RK4 wall only',
    }
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))
    if not all(result['checks'].values()):
        raise SystemExit('Bounded RK4 gap profile failed acceptance gate')


if __name__ == '__main__':
    main()
