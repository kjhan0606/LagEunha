"""Actual native GPU parity under an invertible global-ID relabeling."""
import argparse
import json
from pathlib import Path
import struct
import numpy as np
from audit_kepler_wall_events_20261029 import analyze
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK


def normalized_records(path, scale, offset):
    raw = path.read_bytes()
    header = raw[:16]
    magic, n, stride, stress = struct.unpack('<4i', header)
    if (magic, stride, stress) != (-1, 576, 1) or len(raw) != 16+n*576+52:
        raise ValueError('Not a complete native restart')
    ids = np.ndarray((n,), dtype='<u8', buffer=raw, offset=24, strides=(576,)) & ID_MASK
    if np.any(ids < offset) or np.any((ids-offset) % scale):
        raise ValueError('Unexpected remapped ID')
    canonical = (ids-offset)//scale
    order = np.argsort(canonical)
    if not np.array_equal(canonical[order], np.arange(n)):
        raise ValueError('Changed canonical ID set')
    records = np.frombuffer(raw, dtype=np.uint8, count=n*576, offset=16).reshape(n, 576)[order].copy()
    values = np.ndarray((n,), dtype='<u8', buffer=records, offset=8, strides=(576,))
    values[:] = (values & ~ID_MASK) | canonical[order]
    return header, records, raw[-52:]


def compare(a, b, remap=(1, 0)):
    ha, pa, fa = normalized_records(a, 1, 0)
    hb, pb, fb = normalized_records(b, *remap)
    changed = np.any(pa[:, 8:] != pb[:, 8:], axis=1)
    checks = {'same_header': ha == hb, 'same_full_records_after_ID_inverse_and_pointer_exclusion': not np.any(changed),
              'same_complete_restart_footer': fa == fb}
    return {'passed': all(checks.values()), 'checks': {k: bool(v) for k, v in checks.items()},
            'reference': str(a), 'candidate': str(b), 'ID_scale_offset': remap,
            'changed_records': int(changed.sum()),
            'differing_record_byte_columns': (np.flatnonzero(np.any(pa[:, 8:] != pb[:, 8:], axis=0))+8).tolist(),
            'scope': 'Sort by inverse-mapped global ID, preserve upper flag bits and all auxiliary/RK bytes; exclude only runtime next-pointer bytes 0:8.'}


def endpoint(root):
    files = [p for p in root.glob('khout.*.dat') if p.name != 'khout.011580.dat']
    if len(files) != 1:
        raise ValueError(f'Expected one endpoint: {root}')
    return files[0]


def audit(root):
    v1 = root.parent/'native_events_v1'
    pairs = [(endpoint(v1/'events_no_hit'), endpoint(root/'compact_no_hit'), (1, 0)),
             (endpoint(root/'compact_no_hit'), endpoint(root/'gapped_retry'), (2, 7)),
             (endpoint(v1/'crossing_to484'), endpoint(root/'compact_crossing'), (1, 0)),
             (endpoint(root/'compact_crossing'), endpoint(root/'gapped_crossing'), (2, 7))]
    parity = [compare(a, b, mapping) for a, b, mapping in pairs]
    expected = 2*np.arange(55351, dtype=np.uint64)+7
    results = {'compact': analyze(root/'compact_crossing'),
               'gapped': analyze(root/'gapped_crossing', reader=lambda p: ordered_global(p, expected))}
    retry = (root/'gapped_retry/log').read_text(errors='replace')
    checks = {'all_full_record_parities': all(p['passed'] for p in parity),
              'both_physical_crossings_pass_unchanged_gates': all(p['passed'] for p in results.values()),
              'actual_gapped_retry_restores_state_and_cache': '[WALL_ROLLBACK]' in retry and 'payload_counts_dtold_ledgers_exact=1 device_cache_invalidated=1' in retry,
              'all_native_exits_zero': all((root/name/'exit_status.txt').read_text().strip() == '0' for name in ['compact_no_hit','gapped_retry','compact_crossing','gapped_crossing'])}
    return {'passed': all(checks.values()), 'checks': checks, 'full_record_comparisons': parity,
            'physical_crossings': results, 'scope': 'Native v2 global-ID and endpoint-invariant regression, using the real late Q2 collision. Not a Q1 or long-time stability proof.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path);parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args();result = audit(args.root)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2);stream.write('\n')
    print(json.dumps(result, indent=2))
    if not result['passed']:
        raise SystemExit(1)
