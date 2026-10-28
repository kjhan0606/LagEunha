"""Compare seed5 and seed3 native CUDA outputs and same-device timing."""

import argparse
import json
from pathlib import Path
import re
import statistics

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


FIELDS = ('x', 'y', 'vx', 'vy', 'mass', 'volume', 'ie', 'den', 'pressure', 'w2')
DEVICE = re.compile(r'\[GPU tess device profile\] call=(\d+) density_copy_A=([\d.]+).* ms')
TOTAL = re.compile(r'\[GPU tess profile\] call=(\d+).* total=([\d.]+) ms')


def summarize_profile(log):
    samples = {'kernel_A_ms': [], 'tess_total_ms': []}
    for line in log.splitlines():
        match = DEVICE.fullmatch(line)
        if match:
            samples['kernel_A_ms'].append(float(match[2]))
        match = TOTAL.fullmatch(line)
        if match:
            samples['tess_total_ms'].append(float(match[2]))
    return {name: {'calls': len(values),
                   'median_ms': statistics.median(values) if values else None,
                   'mean_ms': statistics.mean(values) if values else None,
                   'all_ms': values}
            for name, values in samples.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run', type=Path)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    manifest = json.loads((args.run / 'stage_manifest.json').read_text())
    out = {'run': str(args.run.resolve()), 'manifest': manifest,
           'modes': {}, 'checks': {},
           'scope': 'One native RK4 step on one H200 from the identical jitter t=6.5 state. Same-GPU timing is preliminary.'}
    for mode in ('w0', 'weighted'):
        rows = {}
        states = {}
        for version in ('v3', 'v4'):
            case = args.run / f'{mode}_{version}'
            audit = json.loads((case / 'audit_t6501.json').read_text())
            if not all(audit['checks'].values()):
                raise ValueError(f'{mode}_{version}: independent physics audit failed')
            endpoint = Path(audit['endpoint'])
            if (audit['restart_sha256'] != manifest['parent_sha256']
                    or digest(endpoint) != audit['endpoint_sha256']
                    or abs(snapshot_time(endpoint, 55351) - 6.501) > 1e-10):
                raise ValueError(f'{mode}_{version}: lineage or endpoint changed')
            states[version], _ = ordered(endpoint, 55351)
            rows[version] = {'audit': str(case / 'audit_t6501.json'),
                             'endpoint_sha256': audit['endpoint_sha256'],
                             'tess_profile': summarize_profile((case / 'log').read_text(errors='replace'))}
        comparison = {}
        for field in FIELDS:
            baseline = states['v3'][field]
            candidate = states['v4'][field]
            difference = np.abs(baseline.astype(float) - candidate.astype(float))
            comparison[field] = {'bitwise_equal': bool(np.array_equal(baseline, candidate)),
                                 'max_absolute_difference': float(np.max(difference)),
                                 'mean_absolute_difference': float(np.mean(difference)),
                                 'relative_to_max_field': float(np.max(difference) /
                                     max(float(np.max(np.abs(baseline.astype(float)))), 1e-300))}
        timing = {}
        for quantity in ('kernel_A_ms', 'tess_total_ms'):
            before = rows['v3']['tess_profile'][quantity]['median_ms']
            after = rows['v4']['tess_profile'][quantity]['median_ms']
            timing[quantity] = {'v3_median_ms': before, 'v4_median_ms': after,
                                'v3_over_v4': before / after if before and after else None}
        out['modes'][mode] = {'versions': rows, 'physical_comparison': comparison,
                              'timing': timing}
        out['checks'][f'{mode}_same_mass'] = comparison['mass']['bitwise_equal']
        out['checks'][f'{mode}_profiles_present'] = all(
            rows[version]['tess_profile'][quantity]['calls'] >= 4
            for version in ('v3', 'v4')
            for quantity in ('kernel_A_ms', 'tess_total_ms'))
    with args.output.open('x') as stream:
        json.dump(out, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({'checks': out['checks'],
                      'timing': {k: v['timing'] for k, v in out['modes'].items()},
                      'max_relative_field_difference': {
                          k: {field: values['relative_to_max_field']
                              for field, values in row['physical_comparison'].items()}
                          for k, row in out['modes'].items()}}, indent=2))
    if not all(out['checks'].values()):
        raise SystemExit('Matched seed3 comparison incomplete')


if __name__ == '__main__':
    main()
