"""Read-only native SSP diagnostics; run on lageunha, never evolve in Python."""
from pathlib import Path
import argparse
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from audit_snapshot import audit, read_snapshot

GAMMA = float(np.float32(1.666667))
FOCUS = (67935, 68001, 68526, 68592)


def ordered(path):
    p, meta = read_snapshot(path)
    return p[np.argsort(p['u4if'])], meta


def entropy_and_orbit(path):
    p, meta = ordered(path)
    initial = ROOT / 'core_refinement/fine_integrated/khout.000001.dat'
    q, mi = ordered(initial)
    assert np.array_equal(p['u4if'], q['u4if'])
    r0 = np.hypot(q['x'] - 12, q['y'] - 12)
    omega = (r0 * r0 + .047**2)**(-.75)
    theta = np.arctan2(q['y'] - 12, q['x'] - 12) + omega * (meta['t'] - mi['t'])
    vc = r0 * omega
    expected_x, expected_y = 12 + r0*np.cos(theta), 12 + r0*np.sin(theta)
    expected_vx, expected_vy = -vc*np.sin(theta), vc*np.cos(theta)
    uad = q['ie'].astype(float) * (q['volume'].astype(float)/p['volume'])**(GAMMA-1)
    ratio = p['ie']/uad
    ids = p['u4if'] & np.uint64((1 << 56) - 1)
    result = dict(snapshot=str(path), time=meta['t'], n=len(p),
                  mass_exact=np.array_equal(p['mass'], q['mass']),
                  W0_exact=bool(np.all(p['w2'] == 0)), min_U=float(p['ie'].min()),
                  entropy_reference='initial U times (initial V/current V)^(gamma-1); diagnostic adiabat, not a prescribed solution for Riemann exchanges',
                  orbital_reference='analytic circular particle trajectories in fixed Plummer potential; not a substitute for a continuum cell-average comparison',
                  focus={}, radial=[])
    for ident in FOCUS:
        hit = ids == ident
        result['focus'][str(ident)] = dict(U=float(p['ie'][hit][0]),
                Uad=float(uad[hit][0]), entropy_ratio=float(ratio[hit][0]),
                volume_ratio=float((p['volume'][hit]/q['volume'][hit])[0]))
    for lo, hi in ((0, .1), (.1, .2), (.2, .5), (.5, .65), (.65, 1), (1, 2), (2, 4), (4, 8)):
        s = (r0 >= lo) & (r0 < hi)
        if not s.any():
            continue
        mass = p['mass'][s]
        def mw(a):
            return float(np.sum(mass*a)/mass.sum())
        pos_err = np.hypot(p['x'][s]-expected_x[s], p['y'][s]-expected_y[s])/r0[s]
        vel_err = np.hypot(p['vx'][s]-expected_vx[s], p['vy'][s]-expected_vy[s])/vc[s]
        result['radial'].append(dict(initial_Rlo=lo, initial_Rhi=hi, n=int(s.sum()),
                entropy_error_mass_weighted=mw(abs(ratio[s]-1)),
                entropy_ratio_min=float(ratio[s].min()),
                entropy_ratio_p01=float(np.quantile(ratio[s], .01)),
                entropy_ratio_max=float(ratio[s].max()),
                cells_entropy_ratio_below_point1=int((ratio[s] < .1).sum()),
                orbital_position_error_mass_weighted=mw(pos_err),
                orbital_velocity_error_mass_weighted=mw(vel_err),
                internal_energy=float(p['ie'][s].sum())))
    return result


def endpoint_comparison(a_path, b_path):
    a, ma = ordered(a_path)
    b, mb = ordered(b_path)
    assert ma['t'] == mb['t'] and np.array_equal(a['u4if'], b['u4if'])
    result = dict(time=ma['t'], a=str(a_path), b=str(b_path), fields={})
    for key in ('x', 'y', 'vx', 'vy', 'volume', 'ie'):
        da = a[key].astype(float)
        db = b[key].astype(float)
        diff = abs(da-db)
        scale = max(float(abs(da).sum()), float(abs(db).sum()), 1e-300)
        rel = diff/np.maximum(np.maximum(abs(da), abs(db)), 1e-300)
        result['fields'][key] = dict(max_absolute=float(diff.max()),
                L1_relative_sum=float(diff.sum()/scale),
                max_pointwise_symmetric_relative=float(rel.max()),
                p99_pointwise_symmetric_relative=float(np.quantile(rel, .99)))
    return result


def checkpoint(manifest):
    return Path(manifest['final'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fresh', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit('Refuse to overwrite existing analysis')
    intervals = json.loads((ROOT/'rk4_diagnostics/ssprk3_intervals_v2/manifest.json').read_text())
    by_name = {c['name']: c for c in intervals['cases']}
    active = json.loads((ROOT/'rk4_diagnostics/ssprk3_fine_active_half_v2/manifest.json').read_text())['cases'][0]
    result = dict(scope='Short captured intervals and optional fresh run; no t6-10 or paper-duration claim',
                  fine_effective_half=endpoint_comparison(checkpoint(by_name['fine_base']), checkpoint(active)),
                  sedov_ceiling_comparison=endpoint_comparison(checkpoint(by_name['sedov_base']), checkpoint(by_name['sedov_half'])),
                  fine_base=entropy_and_orbit(checkpoint(by_name['fine_base'])),
                  fine_active_half=entropy_and_orbit(checkpoint(active)))
    if args.fresh:
        result['fine_fresh'] = entropy_and_orbit(args.fresh)
        result['fresh_geometry_audit'] = audit(args.fresh)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(output=str(args.output),
            fine_U_L1_relative=result['fine_effective_half']['fields']['ie']['L1_relative_sum'],
            sedov_U_L1_relative=result['sedov_ceiling_comparison']['fields']['ie']['L1_relative_sum'],
            fresh_time=result.get('fine_fresh', {}).get('time')), indent=2), flush=True)


if __name__ == '__main__':
    main()
