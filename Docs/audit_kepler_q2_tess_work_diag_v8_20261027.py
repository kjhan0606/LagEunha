"""Compare diagnostic seed5/seed3 work counts against matched Ada controls."""

import argparse
import json
from pathlib import Path
import re

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


ROOT = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
DIAG = ROOT / 'gpu_tess_work_diag_v8_ada_t6501'
BASE = ROOT / 'gpu_tess_triple_ada_t6501_v1'
FIELDS = ('x', 'y', 'vx', 'vy', 'mass', 'volume', 'ie', 'den', 'pressure', 'w2')
WORK_KEYS = ('pre_visit', 'pre_clip', 'seed_visit', 'seed_clip',
             'adapt_visit', 'adapt_clip', 'wall_clip', 'clip_edges', 'box_retry')
WORK_RE = re.compile(r'\[GPU_TESS_WORK\] (.*)')
KERNEL_RE = re.compile(r'\[GPU tess device profile\] call=(\d+) density_copy_A=([\d.]+)')


def read_case(mode, ref):
    directory = DIAG / mode
    audit = json.loads((directory / 'audit_t6501.json').read_text())
    if len(audit['checks']) != 11 or not all(audit['checks'].values()):
        raise ValueError(f'{mode} short audit did not pass 11/11')
    log = (directory / 'log').read_text(errors='replace')
    work = []
    for match in WORK_RE.finditer(log):
        row = {k: int(v) for k, v in re.findall(r'(\w+)=(\d+)', match.group(1))}
        if set(row) != set(WORK_KEYS):
            raise ValueError(f'{mode} malformed work row: {row}')
        work.append(row)
    kernel_ms = [float(t) for _, t in KERNEL_RE.findall(log)]
    if len(work) != 4 or len(kernel_ms) != 4:
        raise ValueError(f'{mode}: expected four RK4 tessellation calls')
    endpoint = directory / 'khout.000002.dat'
    control = BASE / ref / 'khout.000002.dat'
    a, _ = ordered(control, 55351)
    b, _ = ordered(endpoint, 55351)
    equal = {key: bool(np.array_equal(a[key], b[key])) for key in FIELDS}
    if not all(equal.values()):
        raise ValueError(f'{mode} changed physical endpoint fields')
    return {
        'run': str(directory),
        'binary_sha256': digest(directory / 'eunha2.exe'),
        'log_sha256': digest(directory / 'log'),
        'audit_sha256': digest(directory / 'audit_t6501.json'),
        'endpoint_sha256': digest(endpoint),
        'control_endpoint_sha256': digest(control),
        'physical_fields_bitwise_equal': equal,
        'short_audit_checks_passed': len(audit['checks']),
        'kernel_A_ms': kernel_ms,
        'kernel_A_median_ms': float(np.median(kernel_ms)),
        'work_by_call': work,
        'work_mean': {key: float(np.mean([row[key] for row in work]))
                      for key in WORK_KEYS},
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    cases = {'seed5': read_case('seed5', 'w0_v3'),
             'seed3': read_case('seed3', 'w0_v4')}
    a, b = cases['seed5'], cases['seed3']
    ratio = {key: b['work_mean'][key] / a['work_mean'][key]
             for key in WORK_KEYS}
    result = {
        'time_start': 6.5,
        'time_end': 6.501,
        'parent_sha256': json.loads((DIAG / 'stage_manifest.json').read_text())['parent_sha256'],
        'cases': cases,
        'seed3_over_seed5_work_ratio': ratio,
        'seed3_over_seed5_kernel_A_median_ratio':
            b['kernel_A_median_ms'] / a['kernel_A_median_ms'],
        'scope': 'Same Ada GPU and W0 restart; diagnostic atomics affect absolute timings. Each instrumented binary matches its uninstrumented physical endpoint bitwise.'
    }
    with args.output.open('x') as out:
        json.dump(result, out, indent=2)
        out.write('\n')
    print(json.dumps({'seed5_kernel_ms': a['kernel_A_median_ms'],
                      'seed3_kernel_ms': b['kernel_A_median_ms'],
                      'kernel_ratio': result['seed3_over_seed5_kernel_A_median_ratio'],
                      'work_ratio': ratio}, indent=2))


if __name__ == '__main__':
    main()
