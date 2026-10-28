"""Hash-guarded same-Ada comparison of native q2 tessellation branches."""

import json
from pathlib import Path

import numpy as np

from audit_kepler_q2_seed3_matched_t6501_20261027 import FIELDS, summarize_profile
from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


DOCS = Path(__file__).resolve().parent
ROOT = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
TRIPLE = ROOT / 'gpu_tess_triple_ada_t6501_v1'
LIMITER = ROOT / 'gpu_tess_limiter_v7_ada_t6501_v1'
OUTPUT = DOCS / 'KEPLER_Q2_TESS_ADA_MATCHED_T6501_2026-10-27.json'


def main():
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    staged = {name: json.loads((path / 'stage_manifest.json').read_text())
              for name, path in (('triple', TRIPLE), ('v7', LIMITER))}
    parent_sha = staged['triple']['parent_sha256']
    if staged['v7']['parent_sha256'] != parent_sha:
        raise ValueError('Different accepted checkpoint')
    result = {'source_runs': {'triple': str(TRIPLE), 'v7': str(LIMITER)},
              'parent_sha256': parent_sha,
              'time_start': 6.5, 'time_end': 6.501,
              'cases': {}, 'checks': {},
              'scope': 'One RK4 step from identical 5% site-jitter t=6.5 checkpoint on the same RTX 5000 Ada GPU; no long-time performance claim.'}
    for mode in ('w0', 'weighted'):
        states = {}
        rows = {}
        for version in ('v3', 'v4', 'v7'):
            base = TRIPLE if version != 'v7' else LIMITER
            case = base / f'{mode}_{version}'
            audit_name = ('audit_corrected_t6501.json'
                          if mode == 'weighted' and version in ('v3', 'v4')
                          else 'audit_t6501.json')
            audit_path = case / audit_name
            audit = json.loads(audit_path.read_text())
            if not all(audit['checks'].values()):
                raise ValueError(f'{mode}_{version}: one-step physics audit failed')
            endpoint = Path(audit['endpoint'])
            if (audit['restart_sha256'] != parent_sha
                    or digest(endpoint) != audit['endpoint_sha256']
                    or abs(snapshot_time(endpoint, 55351) - 6.501) > 1e-10):
                raise ValueError(f'{mode}_{version}: lineage or endpoint changed')
            states[version], _ = ordered(endpoint, 55351)
            log = case / 'log'
            profile = summarize_profile(log.read_text(errors='replace'))
            rows[version] = {'run': str(case), 'audit': str(audit_path),
                             'audit_sha256': digest(audit_path),
                             'endpoint_sha256': audit['endpoint_sha256'],
                             'log_sha256': digest(log),
                             'binary_sha256': staged['triple']['binary_sha256'][version]
                                if version != 'v7' else staged['v7']['binary_sha256'],
                             'checks_passed': len(audit['checks']),
                             'kernel_A_median_ms': profile['kernel_A_ms']['median_ms'],
                             'tess_median_ms': profile['tess_total_ms']['median_ms'],
                             'profile_calls': profile['kernel_A_ms']['calls']}
            if profile['kernel_A_ms']['calls'] != 4 or profile['tess_total_ms']['calls'] != 4:
                raise ValueError(f'{mode}_{version}: expected four tessellation timings')
        comparison = {}
        for version in ('v4', 'v7'):
            fields = {field: {'bitwise_equal': bool(np.array_equal(states['v3'][field],
                                                                    states[version][field])),
                              'max_absolute_difference': float(np.max(np.abs(
                                  states['v3'][field].astype(float) -
                                  states[version][field].astype(float))))}
                      for field in FIELDS}
            comparison[version] = {
                'physical_fields': fields,
                'all_ten_bitwise_equal': all(row['bitwise_equal'] for row in fields.values()),
                'kernel_A_v3_over_candidate': rows['v3']['kernel_A_median_ms'] /
                                              rows[version]['kernel_A_median_ms'],
                'tess_v3_over_candidate': rows['v3']['tess_median_ms'] /
                                          rows[version]['tess_median_ms'],
            }
            result['checks'][f'{mode}_{version}_bitwise'] = comparison[version]['all_ten_bitwise_equal']
        failed = TRIPLE / f'{mode}_v6'
        failed_log = (failed / 'log').read_text(errors='replace')
        v6_failure = {'run': str(failed),
                      'binary_sha256': staged['triple']['binary_sha256']['v6'],
                      'log_sha256': digest(failed / 'log'),
                      'geometry_abort': '[GPU_TESS_INVALID] volume=151 velocity_fit=0 geometry=151 halo=0' in failed_log,
                      'accepted_audit_exists': (failed / 'audit_t6501.json').exists()}
        result['checks'][f'{mode}_v6_rejected'] = bool(v6_failure['geometry_abort'] and
                                                      not v6_failure['accepted_audit_exists'])
        result['cases'][mode] = {'accepted_versions': rows,
                                 'versus_v3': comparison,
                                 'failed_v6': v6_failure}
    with OUTPUT.open('x') as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write('\n')
    print(json.dumps({'checks': result['checks'],
                      'timings': {mode: {name: {'kernel_A_ms': row['kernel_A_median_ms'],
                                               'tess_ms': row['tess_median_ms']}
                                        for name, row in case['accepted_versions'].items()}
                                  for mode, case in result['cases'].items()}}, indent=2))
    if not all(result['checks'].values()):
        raise SystemExit('Native tessellation comparison failed')


if __name__ == '__main__':
    main()
