"""Read-only spin/energy ledger from the already captured t=.025 full RHS.

Uses no provider, GPU call, time update, or checkpoint rewrite. The raw 18
moment rates are the first 18 columns of the saved native stage output.
"""
import argparse
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

BASE = Path('/gpfs/kjhan/LagForce')
DEFAULT_CAPTURE = BASE / 'Codex_full_disk_saved_t0p025_flux_budget_20261019'
DEFAULT_CHECKPOINT = BASE / (
    'Codex_full_disk_continuous_two_dt_to_t1_20261019/'
    'coarse/checkpoints/coarse_0005/accepted.npz')
DEFAULT_INITIAL = BASE / 'Codex_full_curved_kepler_IC_20261019/initial.npz'


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def cross(a, b):
    return a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]


def summary(x):
    x = np.asarray(x, dtype=np.longdouble)
    return dict(sum=float(np.sum(x, dtype=np.longdouble)),
                abs_sum=float(np.sum(abs(x), dtype=np.longdouble)),
                max_abs=float(np.max(abs(x))),
                p95_abs=float(np.quantile(np.asarray(abs(x), dtype=float), .95)))


def radial_bounds(cell, radius):
    minimum, maximum, inside = [], [], True
    for line in cell.lines:
        a = line.A
        d = tuple(b - x for x, b in zip(a, line.B))
        norm2 = sum(x * x for x in d)
        normal = (d[1], -d[0])
        if sum(n * (-x) for n, x in zip(normal, a)) > 0:
            inside = False
        t = max(Fraction(0), min(Fraction(1), -sum(x * y for x, y in zip(a, d)) / norm2))
        point = tuple(x + t * y for x, y in zip(a, d))
        minimum.append(sum(x * x for x in point))
        maximum.extend(sum(x * x for x in q) for q in (a, line.B))
    if cell.arcs:
        minimum.append(Fraction(float(radius))**2)
        maximum.append(Fraction(float(radius))**2)
    if inside:
        minimum.append(Fraction(0))
    return float(np.sqrt(np.longdouble(min(minimum)))), float(np.sqrt(np.longdouble(max(maximum))))


def region_masks(checkpoint, initial_path, geometry_dir, ids, mass):
    sys.path.insert(0, str(BASE / 'Codex_circular_curved_modules_v3_20261019'))
    from curved_disk_cell import DiskCell
    with np.load(initial_path) as z:
        initial_ids, initial_mass = z['ids'], z['mass']
        initial_sites = z['state'][:-4].reshape(len(initial_ids), 20)[:, :2]
    assert np.array_equal(ids, initial_ids) and np.array_equal(mass, initial_mass)
    geometry = json.loads((geometry_dir / 'full_exact_supports.json').read_text())
    with np.load(geometry_dir / 'mesh.npz') as z:
        radius = float(z['radius'])
        arcs = z['arc_counts']
    cells = geometry['cells']
    assert len(cells) == len(ids) == len(arcs)
    bounds = np.array([radial_bounds(DiskCell.from_native(cell, radius), radius) for cell in cells])
    lo, hi = bounds.T
    pad = 32 * np.finfo(float).eps * np.maximum(1., hi)
    inner_contact = (lo - pad <= .5) & (hi + pad >= .5)
    outer_contact = (lo - pad <= 2.) & (hi + pad >= 2.)
    masks = dict(all=np.ones(len(ids), dtype=bool), inner=np.linalg.norm(initial_sites, axis=1) < .4,
        cavity=(hi < .5) & ~inner_contact,
        disk_core=(lo > .5) & (hi < 2.) & ~inner_contact & ~outer_contact,
        inner_contact=inner_contact, outer_contact=outer_contact,
        dynamic_buffer=(lo > 2.) & ~(arcs > 0) & ~outer_contact, wall=arcs > 0)
    return masks, radius


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--capture-root', type=Path, default=DEFAULT_CAPTURE)
    p.add_argument('--checkpoint', type=Path, default=DEFAULT_CHECKPOINT)
    p.add_argument('--initial', type=Path, default=DEFAULT_INITIAL)
    p.add_argument('--output', type=Path)
    args = p.parse_args()

    capture_report_path = args.capture_root / 'report.json'
    stage_path = args.capture_root / 'stage_compact.npz'
    capture_report = json.loads(capture_report_path.read_text())
    assert capture_report['status'] == 'completed_single_saved_state_rhs_budget_capture'
    assert capture_report['physical_providers'] == 1
    assert capture_report['accepted_time_steps'] == 0
    assert capture_report['original_sources_unchanged'] is True
    verified_sources = {}
    for name, expected in capture_report['source_sha256'].items():
        observed = sha256(name)
        assert observed == expected, (name, observed, expected)
        verified_sources[name] = observed

    with np.load(args.checkpoint) as z:
        state = z['state']
        mass = z['mass']
        ids = z['ids']
        time = float(z['time'])
    with np.load(stage_path) as z:
        stage = z['stage'].astype(np.longdouble)
        assert np.array_equal(ids, z['ids']) and np.array_equal(mass, z['mass'])
        input_state_sha = str(z['initial_state_sha256'].item())
        ledgers = z['ledgers'].reshape(3, 4).astype(np.longdouble)
    state_sha = hashlib.sha256(np.ascontiguousarray(state).view(np.uint8)).hexdigest()
    assert state_sha == input_state_sha
    n = len(mass)
    assert time == .025 and stage.shape == (n, 94) and n == 19945
    assert np.isfinite(stage[:, :18]).all()

    # Decode raw moments as [m,qx,qy,Qxx,Qxy,Qyx,Qyy,px,py,
    # Dxx,Dxy,Dyx,Dyy,Pixx,Pixy,Piyx,Piyy,E].
    raw = np.column_stack([mass, state[:-4].reshape(n, 20)[:, 3:]]).astype(np.longdouble)
    rate = stage[:, :18]
    m, md = raw[:, 0], rate[:, 0]
    q, qd = raw[:, 1:3], rate[:, 1:3]
    mom, momd = raw[:, 7:9], rate[:, 7:9]
    Q = raw[:, 3:7].reshape(n, 2, 2)
    Qd = rate[:, 3:7].reshape(n, 2, 2)
    D = raw[:, 9:13].reshape(n, 2, 2)
    Dd = rate[:, 9:13].reshape(n, 2, 2)
    Pi = raw[:, 13:17].reshape(n, 2, 2)
    Pid = rate[:, 13:17].reshape(n, 2, 2)
    outer = lambda a, b: np.einsum('ni,nj->nij', a, b)

    S = Q - outer(q, q) / m[:, None, None]
    Sdot = Qd - (outer(qd, q) + outer(q, qd)) / m[:, None, None]
    Sdot += outer(q, q) * md[:, None, None] / m[:, None, None]**2
    C = D - outer(q, mom) / m[:, None, None]
    Cdot = Dd - (outer(qd, mom) + outer(q, momd)) / m[:, None, None]
    Cdot += outer(q, mom) * md[:, None, None] / m[:, None, None]**2
    S = .5 * (S + S.transpose(0, 2, 1))
    det = S[:, 0, 0] * S[:, 1, 1] - S[:, 0, 1] * S[:, 1, 0]
    assert np.all(m > 0) and np.all(det > 0)
    invS = np.empty_like(S)
    invS[:, 0, 0] = S[:, 1, 1] / det
    invS[:, 1, 1] = S[:, 0, 0] / det
    invS[:, 0, 1] = -S[:, 0, 1] / det
    invS[:, 1, 0] = -S[:, 1, 0] / det
    B = np.einsum('nji,njk->nik', C, invS)  # C^T S^-1
    trace = lambda a: np.trace(a, axis1=1, axis2=2)

    I = trace(S)
    Idot = trace(Sdot)
    ell = C[:, 0, 1] - C[:, 1, 0]
    c = q / m[:, None]
    v = mom / m[:, None]
    orbital = cross(c, mom)
    ell_dot = (Dd[:, 0, 1] - Dd[:, 1, 0] - cross(qd, mom) / m
               - cross(c, momd) + orbital * md / m)
    cdot = (qd - c * md[:, None]) / m[:, None]
    orbital_dot = cross(cdot, mom) + cross(c, momd)
    raw_Ldot = Dd[:, 0, 1] - Dd[:, 1, 0]

    omega = ell / I
    Krot = .5 * ell * omega
    torque_power = omega * ell_dot
    inertia_power = -.5 * omega**2 * Idot
    Krot_dot = torque_power + inertia_power

    Kcm_dot = np.sum(v * momd, axis=1) - .5 * np.sum(v * v, axis=1) * md
    Kaffine_dot = (np.einsum('nji,nij->n', B, Cdot)
                   - .5 * np.einsum('nij,njk,nik->n', B, Sdot, B))
    T = Pi - outer(mom, mom) / m[:, None, None]
    Tdot = (Pid - (outer(momd, mom) + outer(mom, momd)) / m[:, None, None]
            + outer(mom, mom) * md[:, None, None] / m[:, None, None]**2)
    Kinternal_dot = .5 * trace(Tdot)
    Kresidual_dot = Kinternal_dot - Kaffine_dot
    thermal_dot = rate[:, 17] - .5 * trace(Pid)
    partition_dot = Kcm_dot + Kaffine_dot + Kresidual_dot + thermal_dot

    masks, radius = region_masks(args.checkpoint, args.initial, args.checkpoint.parent / 'geometry', ids, mass)
    regional = {}
    rates = dict(spin=ell_dot, centroid_orbital=orbital_dot, total_angular=raw_Ldot,
        rotational_projection=Krot_dot, torque_power=torque_power, inertia_change=inertia_power,
        centroid_kinetic=Kcm_dot, affine_kinetic=Kaffine_dot, residual_kinetic=Kresidual_dot,
        thermal=thermal_dot, hydro_energy=rate[:, 17])
    mass_ld = mass.astype(np.longdouble)
    for name, mask in masks.items():
        weight = mass_ld[mask]
        weighted_rms = ({key: float(np.sqrt(np.sum(weight * value[mask]**2, dtype=np.longdouble) / weight.sum()))
                         for key, value in rates.items()} if len(weight) else
                        {key: None for key in rates})
        regional[name] = dict(cells=int(mask.sum()), mass=float(weight.sum()),
            signed_rates={key: float(np.sum(value[mask], dtype=np.longdouble)) for key, value in rates.items()},
            abs_rate_sums={key: float(np.sum(abs(value[mask]), dtype=np.longdouble)) for key, value in rates.items()},
            mass_weighted_RMS=weighted_rms)

    angle_rate = np.sum(raw_Ldot, dtype=np.longdouble)
    angle_ledger = np.sum(ledgers[:, 3], dtype=np.longdouble)
    report = dict(
        status='completed_readonly_saved_full_moment_spin_budget',
        scope='One saved instantaneous RHS at t=.025; no provider replay or state update.',
        cells=n, time=time, physical_providers_replayed=0, accepted_steps=0,
        capture_report_status=capture_report['status'],
        angular_momentum=dict(
            internal_spin=sum_summary(ell),
            internal_spin_rate=sum_summary(ell_dot),
            centroid_orbital_rate=sum_summary(orbital_dot),
            raw_total_rate=sum_summary(raw_Ldot),
            spin_plus_centroid_closure_max_abs=float(np.max(abs(ell_dot + orbital_dot - raw_Ldot))),
            saved_global_ledger_rate=float(angle_ledger),
            global_rate_minus_ledger=float(angle_rate - angle_ledger)),
        rotational_projection=dict(
            energy=sum_summary(Krot), rate=sum_summary(Krot_dot),
            torque_power=sum_summary(torque_power), inertia_change=sum_summary(inertia_power),
            chain_rule_max_abs=float(np.max(abs(Krot_dot - torque_power - inertia_power))),
            thermal_rate=sum_summary(thermal_dot)),
        total_energy_partition_rate=dict(
            centroid_kinetic=sum_summary(Kcm_dot), affine_kinetic=sum_summary(Kaffine_dot),
            residual_kinetic=sum_summary(Kresidual_dot), thermal=sum_summary(thermal_dot),
            partition=sum_summary(partition_dot), stored_hydro=sum_summary(rate[:, 17]),
            max_cell_abs_difference=float(np.max(abs(partition_dot - rate[:, 17])))),
        radial_regions=dict(radius=radius, masks='Exact cell radial bounds from the saved curved mesh; fixed thresholds 0.4, 0.5, and 2.0.',
                            reports=regional),
        mass_rate=dict(max_abs=float(np.max(abs(md))), signed_sum=float(np.sum(md, dtype=np.longdouble))),
        input_sha256=dict(checkpoint=sha256(args.checkpoint), capture_report=sha256(capture_report_path),
                          stage_compact=sha256(stage_path), initial=sha256(args.initial),
                          disk_cell_geometry_module=sha256(BASE/'Codex_circular_curved_modules_v3_20261019/curved_disk_cell.py'),
                          saved_state_payload=state_sha,
                          verified_sources=verified_sources),
        limits=[
            'Instantaneous tangent only; no finite-step relaxation rate or t6-10 claim.',
            'Krot is a rigid-rotation projection inside affine kinetic energy, not additional stored energy.',
            'Its derivative is not heat: torque work and inertia-change work both contribute.',
            'Regional sums use overlapping named masks; they are diagnostics, not disjoint partitions.',
            'The saved RHS already existed; no new physical provider or accepted step was executed.'
        ])
    text = json.dumps(report, indent=2, allow_nan=False) + '\n'
    if args.output:
        assert not args.output.exists(), f'refuse overwrite: {args.output}'
        args.output.write_text(text)
    print(text, end='')


def sum_summary(x):
    return summary(x)


if __name__ == '__main__':
    main()
