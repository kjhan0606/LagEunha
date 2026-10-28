"""Measure signed and absolute disk entropy change in accepted q2 IC controls."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from audit_kepler_q2_three_ic_native_t65_20261027 import SOURCES
from prepare_kepler_sparse_area_restart_20261027 import digest


GAMMA_MINUS_ONE = float(np.float32(1.666667)) - 1.


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    cases = {}
    for label, source in SOURCES.items():
        d = json.loads(source.read_text())
        if d['time'] != 6.5 or not all(d['checks'].values()):
            raise ValueError(f'{label}: endpoint not accepted')
        p0_path, p1_path = Path(d['initial']), Path(d['endpoint'])
        if digest(p0_path) != d['initial_sha256'] or digest(p1_path) != d['endpoint_sha256']:
            raise ValueError(f'{label}: source snapshot changed')
        p0, _ = ordered(p0_path, 55351)
        p1, _ = ordered(p1_path, 55351)
        r0 = np.hypot(p0['x'].astype(float)-3., p0['y'].astype(float)-3.)
        disk = (r0 >= .5) & (r0 < 2.)
        m = p0['mass'][disk].astype(float)
        s = (np.log(p1['ie'][disk].astype(float)/p0['ie'][disk].astype(float))
             + GAMMA_MINUS_ONE * np.log(
                 p1['volume'][disk].astype(float)/p0['volume'][disk].astype(float)))
        mean = float(np.sum(s*m, dtype=np.longdouble)/np.sum(m, dtype=np.longdouble))
        mean_abs = float(np.sum(np.abs(s)*m, dtype=np.longdouble)/np.sum(m, dtype=np.longdouble))
        if abs(mean_abs-d['disk_entropy_log_abs_mass_weighted']) > 1e-12:
            raise ValueError(f'{label}: mismatch with endpoint audit')
        cases[label] = {'disk_sites': int(np.sum(disk)),
                        'signed_mass_weighted_mean': mean,
                        'absolute_mass_weighted_mean': mean_abs,
                        'rms_mass_weighted': float(np.sqrt(
                            np.sum(s*s*m, dtype=np.longdouble)/np.sum(m, dtype=np.longdouble))),
                        'fraction_positive_cells': float(np.mean(s > 0)),
                        'unweighted_p01': float(np.percentile(s, 1)),
                        'unweighted_p50': float(np.percentile(s, 50)),
                        'unweighted_p99': float(np.percentile(s, 99)),
                        'audit_sha256': digest(source)}
    result = {'time': 6.5, 'gamma_minus_one': GAMMA_MINUS_ONE,
              'cases': cases,
              'scope': 'Cell entropy proxy log(ie/ie0)+(gamma-1)log(V/V0) over initial disk cells; signed mean distinguishes net drift from absolute local changes, but does not by itself identify a unique dissipative operator.'}
    with args.output.open('x') as out:
        json.dump(result, out, indent=2)
        out.write('\n')
    print(json.dumps(cases, indent=2))


if __name__ == '__main__':
    main()
