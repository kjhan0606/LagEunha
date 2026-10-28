"""Exact full native restart comparison, excluding only the linked-list next pointer."""
import argparse
import hashlib
import json
from pathlib import Path
import struct

import numpy as np


def compare(a, b):
    a, b = Path(a), Path(b)
    hashes = [hashlib.sha256(), hashlib.sha256()]
    counts = np.zeros(576, dtype=np.int64)
    changed_rows = pointer_rows = 0
    first_rows = []
    with a.open('rb') as left, b.open('rb') as right:
        headers = [left.read(16), right.read(16)]
        if any(len(x) != 16 for x in headers):
            raise ValueError('Incomplete native header')
        layouts = [struct.unpack('<4i', x) for x in headers]
        for path, (magic, n, stride, stress) in zip([a, b], layouts):
            if magic != -1 or n <= 0 or stride != 576 or stress != 1:
                raise ValueError('Expected full native stress-RK4 ABI, stride 576')
            if path.stat().st_size != 16+n*stride+52:
                raise ValueError('Incomplete or trailing native restart data')
        if layouts[0] != layouts[1]:
            raise ValueError('Different native layouts or record counts')
        n = layouts[0][1]
        for digest, raw in zip(hashes, headers):
            digest.update(raw)
        for start in range(0, n, 8192):
            count = min(8192, n-start)
            raw = [left.read(count*576), right.read(count*576)]
            for digest, data in zip(hashes, raw):
                digest.update(data)
            x, y = [np.frombuffer(data, dtype=np.uint8).reshape(count, 576) for data in raw]
            delta = x != y
            counts += np.count_nonzero(delta, axis=0)
            rows = np.flatnonzero(np.any(delta[:, 8:], axis=1))
            changed_rows += len(rows)
            pointer_rows += int(np.count_nonzero(np.any(delta[:, :8], axis=1)))
            first_rows.extend((rows[:max(0, 16-len(first_rows))]+start).tolist())
        footers = [left.read(52), right.read(52)]
        for digest, data in zip(hashes, footers):
            digest.update(data)
    checks = {'native_headers_exact': headers[0] == headers[1],
              'all_record_bytes_after_next_pointer_exact': changed_rows == 0,
              'full_precision_restart_footer_exact': footers[0] == footers[1]}
    return {'passed': all(checks.values()), 'checks': checks,
            'paths': [str(a), str(b)], 'sha256': [d.hexdigest() for d in hashes],
            'records': n, 'stride': 576, 'changed_records_excluding_pointer': changed_rows,
            'first_changed_record_indices': first_rows,
            'differing_body_byte_columns': {str(i): int(counts[i]) for i in range(8, 576) if counts[i]},
            'records_with_different_next_pointer': pointer_rows,
            'footer_fields': [struct.unpack('<ddfiiddff', x) for x in footers],
            'scope': 'Same native ABI and record order; excludes only byte 0:8, the next pointer '
                     'declared first in treevorostressrk4particletype. All remaining record bytes, '
                     'including auxiliary/gradient/RK fields, and all 52 footer bytes are exact. '
                     'This compares stored data; it is not an executed restart bridge or long-time proof.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=Path)
    parser.add_argument('candidate', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.reference, args.candidate)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))
    if not result['passed']:
        raise SystemExit('Full restart bytes differ; evidence retained')
