"""Compare three accepted bounded q2 Kepler ICs at the same native time."""

import argparse
import json
from pathlib import Path

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
ROOT = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
SOURCES = {
    'aligned': DOCS / 'KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json',
    'ring_coherent_5pct': ROOT / 'ringcoherent005_to65_a10_v1/audit_t65.json',
    'independent_site_jitter_5pct': ROOT / 'jitter005_to65_a10_v1/audit_t65.json',
}
METRICS = (
    'disk_density_abs_mass_weighted',
    'disk_ballistic_density_abs_mass_weighted',
    'disk_native_minus_ballistic_abs_mass_weighted',
    'disk_entropy_log_abs_mass_weighted',
    'disk_position_error_abs_mass_weighted',
    'disk_radial_speed_rms_over_circular',
    'outer_density_abs_mass_weighted',
    'total_energy_relative',
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    cases = {}
    for label, source in SOURCES.items():
        d = json.loads(source.read_text())
        if (d['time'] != 6.5 or d['n'] != 55351 or d['disk_sites'] != 49038
                or len(d['checks']) != 17 or not all(d['checks'].values())):
            raise ValueError(f'{label}: incompatible or failed native t=6.5 audit')
        if digest(Path(d['endpoint'])) != d['endpoint_sha256']:
            raise ValueError(f'{label}: endpoint changed')
        cases[label] = {
            'source_audit': str(source.resolve()),
            'source_audit_sha256': digest(source),
            'endpoint': d['endpoint'],
            'endpoint_sha256': d['endpoint_sha256'],
            'checks_passed': len(d['checks']),
            'metrics': {key: d[key] for key in METRICS},
        }
    a = cases['aligned']['metrics']
    r = cases['ring_coherent_5pct']['metrics']
    j = cases['independent_site_jitter_5pct']['metrics']
    ratios = {key: {'ring_over_aligned': r[key]/a[key],
                    'jitter_over_ring': j[key]/r[key]}
              for key in METRICS}
    result = {
        'time': 6.5,
        'sites': 55351,
        'disk_sites': 49038,
        'cases': cases,
        'ratios': ratios,
        'scope': 'Three valid bounded-q2 fixed-mass W0 CUDA endpoints; ICs sample the same continuum profile but not identical discrete masses. This is a mechanism control, not a GIZMO-equivalent benchmark.'
    }
    with args.output.open('x') as out:
        json.dump(result, out, indent=2)
        out.write('\n')
    print(json.dumps({'density': {k: v['metrics']['disk_density_abs_mass_weighted']
                                  for k,v in cases.items()},
                      'native_minus_ballistic':
                          {k: v['metrics']['disk_native_minus_ballistic_abs_mass_weighted']
                           for k,v in cases.items()},
                      'ratios': {k: ratios[k] for k in (
                          'disk_density_abs_mass_weighted',
                          'disk_native_minus_ballistic_abs_mass_weighted',
                          'disk_position_error_abs_mass_weighted')}}, indent=2))


if __name__ == '__main__':
    main()
