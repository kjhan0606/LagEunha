"""Saved-state Eulerian Kepler reference; no fit, RHS, or evolution.

Run on lageunha. q32 refines interface-crossing cells and the union of
q16 maximum-error cells. Scientific judgments belong to the Astra review.
"""
import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

BASE = Path('/gpfs/kjhan/LagForce')
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--checkpoint', type=Path, required=True)
parser.add_argument('--initial', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
# Immutable geometry and field implementations used by the saved calculation.
sys.path[:0] = [str(BASE/'Codex_contact_prior_evolution_modules_20261018'),
               str(BASE/'Codex_steady_ale_reference_modules_20261018'), str(BASE)]
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
from steady_kepler_ale_reference import field_moments, polygon_integral
from audit_canonical_pair_history import gravity_coefficient
from run_canonical_av_branches_20261018 import read


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def polygon_area(poly):
    p = poly-poly[0]
    return .5*float(np.sum(p[:, 0]*np.roll(p[:, 1], -1)-p[:, 1]*np.roll(p[:, 0], -1)))


def crosses(poly, radius):
    # Convex polygon minimum distance includes edge interiors, not just vertices.
    d = np.roll(poly, -1, axis=0)-poly
    denom = np.sum(d*d, axis=1)
    t = np.zeros(len(d))
    good = denom > 0
    t[good] = np.clip(-np.sum(poly[good]*d[good], axis=1)/denom[good], 0, 1)
    minimum = np.linalg.norm(poly+t[:, None]*d, axis=1).min()
    side = d[:, 0]*(-poly[:, 1])-d[:, 1]*(-poly[:, 0])
    if np.all(side >= 0):
        minimum = 0.
    return bool(minimum <= radius <= np.linalg.norm(poly, axis=1).max())


def derived(raw, areas):
    assert np.isfinite(raw).all() and np.isfinite(areas).all()
    assert np.all(raw[:, 0] > 0) and np.all(areas > 0)
    m = raw[:, 0]
    q = raw[:, 1:3]
    Q = raw[:, 3:7].reshape(-1, 2, 2)
    p = raw[:, 7:9]
    D = raw[:, 9:13].reshape(-1, 2, 2)
    Pi = raw[:, 13:17].reshape(-1, 2, 2)
    S = Q-q[:, :, None]*q[:, None, :]/m[:, None, None]
    C = D-q[:, :, None]*p[:, None, :]/m[:, None, None]
    B = np.linalg.solve(S, C).transpose(0, 2, 1)
    R = Pi-p[:, :, None]*p[:, None, :]/m[:, None, None]-B@S@B.transpose(0, 2, 1)
    thermal = raw[:, 17]-.5*np.trace(Pi, axis1=1, axis2=2)
    result = dict(mass=m, q=q, Q=Q, momentum=p, D=D, Pi=Pi,
                energy=raw[:, 17], centroid=q/m[:, None], velocity=p/m[:, None],
                ell=C[:, 0, 1]-C[:, 1, 0], R=(R+R.transpose(0, 2, 1))/2,
                thermal=thermal, pressure=(gamma-1)*thermal/areas)
    assert all(np.isfinite(value).all() for value in result.values())
    return result


def statistics(delta, mask):
    indices = np.flatnonzero(mask)
    if not len(indices):
        return dict(cells=0)
    error = np.linalg.norm(delta.reshape(len(ids), -1), axis=1)
    k = indices[np.argmax(error[indices])]
    return dict(cells=len(indices), mass_weighted_RMS=float(np.sqrt(np.average(error[indices]**2, weights=mass[indices]))),
                maximum=float(error[k]), maximum_id=int(ids[k]))


saved = np.load(args.checkpoint)
initial = np.load(args.initial)
ids = saved['ids']
mass = saved['mass']
assert np.array_equal(ids, initial['ids']) and np.array_equal(mass, initial['mass'])
a = saved['state'][:-4].reshape(-1, 20)
a0 = initial['state'][:-4].reshape(-1, 20)
assert float(initial['time']) == 0.
frozen = np.linalg.norm(a0[:, :2], axis=1) >= .60
if 'frozen' in saved:
    assert np.array_equal(frozen, saved['frozen'])
ic, time0 = read(BASE/'Codex_canonical_initial_20261017/run/khout.000000.dat')
assert time0 == 0.
allx = np.column_stack([ic['x']-3, ic['y']-3])
assert np.array_equal(allx[ids], a0[:, :2])
env = json.loads((BASE/'Codex_canonical_initial_20261017/run/environment.json').read_text())
assert env['EUNHA_KEPLER_BG_ROTATE'] == '1' and env['EUNHA_KEPLER_PROFILE'] == 'pakmor'
rin, rout, rc, floor = [float(env['EUNHA_KEPLER_'+key]) for key in ('RIN', 'ROUT', 'RC', 'FLOOR')]
gamma = float(np.float32(1.666667))
pressure = float(ic['pressure'][0])
radii = tuple(sorted(set((rin, rout, rc/2, rc))))
clip_path = BASE/'Codex_weighted_compression_geometry_20261017/audit.py'
potential_path = BASE/'Codex_weight_remap_moments_centered_20261017/audit.py'
clip = load('saved_clip', clip_path)
potential = load('saved_potential', potential_path).potential
builder = EmbeddedLaguerreGeometry(allx, ids, clip.clip)
polys, edges = builder(a[:, :2], a[:, 2], float(saved['time']))
areas = np.array([polygon_area(p) for p in polys])
assert np.all(areas > 0)
raw = np.column_stack([mass, a[:, 3:]])
numerical = derived(raw, areas)
radius0 = np.linalg.norm(a0[:, :2], axis=1)
contact = np.array([crosses(p, rin) or crosses(p, rout) for p in polys])
split = np.array([any(crosses(p, r) for r in radii) for p in polys])
boundary = np.zeros(len(ids), dtype=bool)
for i, j, *_ in edges:
    if j >= 0 and frozen[i] != frozen[j]:
        boundary[i] = boundary[j] = True
regions = dict(all=np.ones(len(ids), dtype=bool), active=~frozen, frozen=frozen,
               inner=radius0 < .4, inner_transition=(radius0 >= .4)&(radius0 < .48),
               contact_band=(radius0 >= .48)&(radius0 < .52),
               outer_interior=(radius0 >= .52)&(radius0 < .55),
               active_edge=(radius0 >= .55)&(radius0 < .60),
               contact_cells=contact, active_frozen_boundary_cells=boundary)
for name in ('contact_cells', 'active_frozen_boundary_cells'):
    regions[name+'_active'] = regions[name] & ~frozen
    regions[name+'_frozen'] = regions[name] & frozen
args.output.mkdir(exist_ok=False)


def integrand(points):
    field, _ = field_moments(points, gravity_coefficient, potential, pressure, gamma, rin, rout, floor)
    return np.column_stack([field, np.ones(len(points))])


values = []
estimates = []
for i, poly in enumerate(polys):
    val, estimate = polygon_integral(poly, integrand, radii, 16)
    values.append(val)
    estimates.append(float(estimate))
    if (i+1) % 100 == 0:
        print(json.dumps(dict(order=16, cells=i+1)), flush=True)
values = np.array(values)
reference = derived(values[:, :18], values[:, 19])
factor_R = saved['factors']@saved['factors'].transpose(0, 2, 1)


def differences(ref):
    delta = {name: numerical[name]-val for name, val in ref.items()}
    delta['pressure_relative_P0'] = delta['pressure']/pressure
    delta['mass_relative_reference'] = delta['mass']/ref['mass']
    delta['factor_R'] = factor_R-ref['R']
    return delta


metrics = {name: {region: statistics(delta, mask) for region, mask in regions.items()}
           for name, delta in differences(reference).items()}
refine = split.copy()
for metric in metrics.values():
    for row in metric.values():
        if row['cells']:
            refine[np.flatnonzero(ids == row['maximum_id'])] = True
refined = values.copy()
q32_estimates = []
for k, i in enumerate(np.flatnonzero(refine)):
    refined[i], estimate = polygon_integral(polys[i], integrand, radii, 32)
    q32_estimates.append(float(estimate))
    if (k+1) % 100 == 0:
        print(json.dumps(dict(order=32, cells=k+1)), flush=True)
reference32 = derived(refined[:, :18], refined[:, 19])
metrics32 = {name: {region: statistics(delta, mask) for region, mask in regions.items()}
             for name, delta in differences(reference32).items()}
quad_delta = {name: reference32[name]-ref for name, ref in reference.items()}
quad_delta['pressure_relative_P0'] = quad_delta['pressure']/pressure
quad_delta['mass_relative_reference'] = quad_delta['mass']/reference['mass']
quad_delta['factor_R'] = quad_delta['R']
quadrature = {name: {region: statistics(delta, mask & refine)
                    for region, mask in regions.items()} for name, delta in quad_delta.items()}
offsets = np.r_[0, np.cumsum([len(p) for p in polys])]
np.savez(args.output/'reference.npz', ids=ids, mass=mass, raw=raw, reference_q16=values,
         reference_q32_selected=refined, refined_mask=refine, areas=areas, frozen=frozen,
         initial_radius=radius0, contact_mask=contact, boundary_mask=boundary,
         polygon_vertices=np.concatenate(polys), polygon_offsets=offsets,
         factors=saved['factors'], time=saved['time'])
report = dict(scope=__doc__, time=float(saved['time']), checkpoint=str(args.checkpoint),
              checkpoint_sha256=sha(args.checkpoint), initial=str(args.initial), initial_sha256=sha(args.initial),
              physical=dict(gamma=gamma, pressure=pressure, rin=rin, rout=rout, floor=floor, rc=rc, split_radii=radii),
              numerical_minus_q16=metrics, q32_minus_q16_selected=quadrature,
              numerical_minus_q32_selected_else_q16=metrics32,
              refined_cells=int(refine.sum()), interface_cells=int(split.sum()),
              q16_maximum_quad_vec_estimate=max(estimates), q32_maximum_quad_vec_estimate=max(q32_estimates, default=0.),
              geometric_area_difference_max=float(abs(values[:, 19]-areas).max()),
              raw_R_minus_factor_R_max=float(abs(numerical['R']-factor_R).max()),
              covariance_note='R metric derives raw S/C/Pi and includes cancellation; factor_R uses saved LL^T minus reference R. mass_relative_reference normalizes by q16 exact mass; pressure_relative_P0 normalizes by IC P0.',
              exact_field_potential_sum_q16={region: float(values[mask, 18].sum())
                                            for region, mask in regions.items()},
              exact_reference_pressure_minus_P0_max_q16=float(abs(reference['pressure']-pressure).max()),
              exact_reference_pressure_minus_P0_max_q32_selected=float(abs(reference32['pressure'][refine]-pressure).max()),
              numerical_potential_state='absent: exact-field Psi is a reference integral, not numerical evolved Psi',
              source_sha256={str(p): sha(p) for p in (Path(__file__), Path(sys.modules['steady_kepler_ale_reference'].__file__), Path(sys.modules['embedded_laguerre_geometry'].__file__), Path(sys.modules[gravity_coefficient.__module__].__file__), BASE/'Codex_canonical_initial_20261017/run/khout.000000.dat', BASE/'Codex_canonical_initial_20261017/run/environment.json', clip_path, potential_path)},
              limits='No material-centroid orbit scoring; fixed numerical mass is not substituted into reference. RMS uses numerical masses. q32 outside selected mask is a copy of q16, not independent integration. No order, long-time, closure-ranking or conservation verdict.')
(args.output/'report.json').write_text(json.dumps(report, indent=2))
print(json.dumps(dict(status='completed', report=str(args.output/'report.json'))), flush=True)
