"""Compare W0 and weighted bin-certificate on/off/fallback native replays."""

import argparse
import json
from pathlib import Path
import re

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest

N = 55351
FIELDS = ('x', 'y', 'vx', 'vy', 'mass', 'volume', 'ie', 'den', 'pressure', 'w2')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('root', type=Path)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    results = {}
    all_ok = True
    for case in ('w0', 'weighted'):
        rows = {}
        arrays = {}
        for mode in ('baseline', 'enabled', 'fallback'):
            run = args.root/case/mode
            receipt = json.loads((run/'audit.json').read_text())
            if not all(receipt['checks'].values()):
                raise ValueError(f'{case}/{mode} physics audit failed')
            endpoint = Path(receipt['endpoint'])
            parent = Path(receipt.get('restart', receipt.get('parent')))
            parent_sha = receipt.get('restart_sha256', receipt.get('parent_sha256'))
            if digest(endpoint) != receipt['endpoint_sha256'] or digest(parent) != parent_sha:
                raise ValueError(f'{case}/{mode} snapshot changed')
            arrays[mode], _ = ordered(endpoint, N)
            log = (run/'log').read_text(errors='replace')
            events = [(int(a), int(b), int(c)) for a,b,c in re.findall(
                r'\[GPU_SEED_BIN_CERT\] invalid=(\d+) skipped_bins=(\d+) forced=(\d+)', log)]
            if mode == 'enabled' and (not events or any(a or c for a,_,c in events)
                                      or not any(b>0 for _,b,_ in events)):
                raise ValueError(f'{case}: enabled certificate was not demonstrably active: {events}')
            if mode == 'fallback' and (not events or any(a!=1 or b!=0 or c!=1
                                                       for a,b,c in events)):
                raise ValueError(f'{case}: forced fallback did not disable all skips: {events}')
            rows[mode] = {'receipt': str((run/'audit.json').resolve()),
                          'receipt_sha256': digest(run/'audit.json'),
                          'restart_sha256': parent_sha, 'endpoint_sha256': receipt['endpoint_sha256'],
                          'telemetry': events}
        if len({row['restart_sha256'] for row in rows.values()}) != 1:
            raise ValueError(f'{case}: variant restart mismatch')
        comparison = {}
        for mode in ('enabled', 'fallback'):
            comparison[mode] = {}
            for name in FIELDS:
                x,y=arrays['baseline'][name],arrays[mode][name]
                equal=bool(np.array_equal(x,y))
                comparison[mode][name] = {'bitwise_equal': equal,
                    'max_abs_difference': float(np.max(np.abs(x.astype(float)-y.astype(float))))}
                all_ok &= equal
        results[case] = {'runs': rows, 'fields': comparison}
    result = {'cases': results, 'checks': {'all_native_gates': True,
              'enabled_and_forced_fallback_telemetry': True,
              'ten_fields_bitwise_equal': bool(all_ok)},
              'scope': 'Short same-checkpoint W0 and weighted GPU replays. '
                       'Forced fallback validates the flag route, not an actual corrupt CSR. '
                       'No long-time equivalence or performance conclusion.'}
    with args.output.open('x') as out:
        json.dump(result,out,indent=2,allow_nan=False)
        out.write('\n')
    print(json.dumps({'checks': result['checks']},indent=2))
    if not all_ok:
        raise SystemExit('Seed certificate changed physical fields')


if __name__ == '__main__':
    main()
