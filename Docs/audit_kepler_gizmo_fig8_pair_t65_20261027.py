"""Apply the published Fig. 8 annulus formulas to same-IC native t=6.5 states."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_gizmo_fig8_baseline_20261027 import annulus_metrics
from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


def read_case(path):
    audit = json.loads(path.read_text())
    if abs(audit['time'] - 6.5) > 1e-12 or not all(audit['checks'].values()):
        raise ValueError(f'Unaccepted t=6.5 native state: {path}')
    initial, endpoint = Path(audit['initial']), Path(audit['endpoint'])
    if (digest(initial) != audit['initial_sha256'] or
            digest(endpoint) != audit['endpoint_sha256']):
        raise ValueError(f'Changed native snapshot: {path}')
    p0, _ = ordered(initial, 55351)
    p, _ = ordered(endpoint, 55351)
    r = np.hypot(p['x'].astype(float)-3., p['y'].astype(float)-3.)
    r0 = np.hypot(p0['x'].astype(float)-3., p0['y'].astype(float)-3.)
    return audit, {
        'receipt': str(path),
        'receipt_sha256': digest(path),
        'endpoint_sha256': audit['endpoint_sha256'],
        'current_radius_annulus': annulus_metrics(p, p0, (r>.8)&(r<1.2)),
        'initial_radius_annulus': annulus_metrics(p, p0, (r0>.8)&(r0<1.2)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--w0', type=Path, required=True)
    ap.add_argument('--weighted', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    wa, w = read_case(args.w0)
    qa, q = read_case(args.weighted)
    if wa['initial_sha256'] != qa['initial_sha256']:
        raise ValueError('Different initial Kepler ICs')
    result = {
        'time': 6.5,
        'initial_sha256': wa['initial_sha256'],
        'paper_url': 'https://arxiv.org/pdf/1409.7395',
        'paper_annulus': 'current-position 0.8<R<1.2, arithmetic site means',
        'cases': {'W0': w, 'weighted': q},
        'scope': 'Published annulus formulas applied to a same-IC bounded polar q2 comparison, not the paper uniform-site or boundary setup and not a t=600 claim.',
    }
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    for name, row in result['cases'].items():
        print(name, row['current_radius_annulus'])


if __name__ == '__main__':
    main()
