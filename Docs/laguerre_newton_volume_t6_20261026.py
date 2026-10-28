"""Solve target power-cell volumes at fixed ballistic Kepler t=6 geometry.

Diagnostic only. The target is each cell's initial Voronoi area; outside the
active region (1 <= initial r < 8) weights remain zero. It does not implement
or validate weighted-face hydrodynamics or a fast evolving weight solver.
"""

from pathlib import Path
import json
import os
import sys
import time

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import LinearOperator, cg
from scipy.spatial import Voronoi

from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
sys.path.insert(0, str(ROOT / 'rk4_diagnostics/varpress_w0_geom_cfl_long_20261026'))
from audit_varpress_checkpoint_20261026 import load_complete


def matrix_from_faces(xy, lookup, nactive, pairs, lengths):
    """Exact local derivative d(power area)/d(site weight) for finite faces."""
    i, j = pairs[:, 0], pairs[:, 1]
    ii, jj = lookup[i], lookup[j]
    touching = (ii >= 0) | (jj >= 0)
    i, j, ii, jj, lengths = (v[touching] for v in (i, j, ii, jj, lengths))
    conductance = lengths/(2.*np.linalg.norm(xy[i]-xy[j], axis=1))
    diagonal = np.zeros(nactive)
    use_i = ii >= 0
    use_j = jj >= 0
    np.add.at(diagonal, ii[use_i], conductance[use_i])
    np.add.at(diagonal, jj[use_j], conductance[use_j])
    inside = use_i & use_j
    rows = np.concatenate((ii[inside], jj[inside]))
    cols = np.concatenate((jj[inside], ii[inside]))
    values = -np.tile(conductance[inside], 2)
    jac = coo_matrix((values, (rows, cols)), shape=(nactive, nactive)).tocsr()
    jac += diags(diagonal, format='csr')
    return jac, diagonal


def newton_jacobian(xy, lookup, nactive, triangles, centers):
    pairs = np.vstack((triangles[:, [0, 1]], triangles[:, [1, 2]],
                       triangles[:, [2, 0]]))
    ids = np.tile(np.arange(len(triangles)), 3)
    pairs.sort(axis=1)
    order = np.lexsort((pairs[:, 1], pairs[:, 0]))
    pairs = pairs[order]
    ids = ids[order]
    begins = np.r_[0, 1+np.flatnonzero(np.any(pairs[1:] != pairs[:-1], axis=1))]
    counts = np.diff(np.r_[begins, len(pairs)])
    touch = (lookup[pairs[begins, 0]] >= 0) | (lookup[pairs[begins, 1]] >= 0)
    if np.any(counts[touch] != 2):
        raise ValueError('An active power cell has a non-finite or nonmanifold face')
    begins = begins[counts == 2]
    lengths = np.linalg.norm(centers[ids[begins]]-centers[ids[begins+1]], axis=1)
    return matrix_from_faces(xy, lookup, nactive, pairs[begins], lengths)


def solve(jac, diagonal, rhs):
    iterations = [0]
    vector, info = cg(jac, rhs,
                      M=LinearOperator(jac.shape, matvec=lambda v:v/diagonal),
                      rtol=1e-10, atol=0., maxiter=2000,
                      callback=lambda _:iterations.__setitem__(0, iterations[0]+1))
    if info != 0:
        raise RuntimeError(f'CG failed: info={info}, iterations={iterations[0]}')
    return vector, iterations[0]


def metrics(areas, target, band):
    relative = abs(areas-target)/target
    return dict(relative_l2=float(np.linalg.norm(relative)),
                median_rel=float(np.median(relative)),
                p95_rel=float(np.quantile(relative, .95)),
                max_rel=float(relative.max()),
                inner_r2_2p25_median_rel=float(np.median(relative[band])),
                inner_r2_2p25_p95_rel=float(np.quantile(relative[band], .95)),
                inner_r2_2p25_max_rel=float(relative[band].max()))


def main():
    started = time.monotonic()
    target_time = float(os.environ.get('LAGUERRE_TIME', '6'))
    if target_time == 6.:
        time_label = 'T6'
    elif abs(target_time-6.250024729493809) < 1e-12:
        time_label = 'T625'
    else:
        raise ValueError('Diagnostic supports only exact t=6 or t=6.250024729493809')
    p0, meta = load_complete(ROOT / 'core_refinement/fine_integrated/khout.000001.dat')
    if abs(float(meta['t'])) > 1e-12:
        raise ValueError('Expected t=0 snapshot')
    xy0 = np.column_stack((p0['x'].astype(float)-12., p0['y'].astype(float)-12.))
    r0 = np.linalg.norm(xy0, axis=1)
    geometry_mode = os.environ.get('LAGUERRE_GEOMETRY', 'ballistic')
    if geometry_mode == 'ballistic':
        theta = target_time*(r0*r0+.047**2)**(-.75)
        xy = np.column_stack((xy0[:, 0]*np.cos(theta)-xy0[:, 1]*np.sin(theta),
                              xy0[:, 0]*np.sin(theta)+xy0[:, 1]*np.cos(theta)))
        source = 'exact softened circular-orbit generator positions'
    elif geometry_mode == 'actual':
        suffix_path = ('C010/to6/khout.032479.dat' if time_label == 'T6' else
                       'C010/to10/khout.036018.dat')
        checkpoint = ROOT / 'rk4_diagnostics/varpress_w0_geom_cfl_long_20261026' / suffix_path
        observed, checkpoint_meta = load_complete(checkpoint)
        if not np.array_equal(observed['u4if'], p0['u4if']):
            raise ValueError('Snapshot particle IDs differ from t=0')
        if abs(float(checkpoint_meta['t'])-target_time) > 1e-12:
            raise ValueError('Unexpected native checkpoint time')
        xy = np.column_stack((observed['x'].astype(float)-12.,
                              observed['y'].astype(float)-12.))
        source = str(checkpoint)
    else:
        raise ValueError('LAGUERRE_GEOMETRY must be ballistic or actual')
    inner_radius = float(os.environ.get('LAGUERRE_INNER_R', '1'))
    if inner_radius not in (1., 1.9, 2.):
        raise ValueError('Diagnostic supports inner radius 1, 1.9, or 2')
    active = np.flatnonzero((r0 >= inner_radius) & (r0 < 8.))
    lookup = np.full(len(xy), -1, dtype=int)
    lookup[active] = np.arange(len(active))
    band = (r0[active] >= 2.) & (r0[active] < 2.25)
    target = p0['volume'][active].astype(float)
    weights = np.zeros(len(xy))
    # Start with the exact Voronoi area, its edge conductances, and one CG step.
    vor = Voronoi(xy, qhull_options='Qbb Qc Qz')
    areas = np.empty(len(active))
    diagonal = np.zeros(len(active))
    rows, cols, values = [], [], []
    for k, site in enumerate(active):
        vertices = vor.vertices[vor.regions[vor.point_region[site]]]
        areas[k] = .5*abs(np.dot(vertices[:, 0], np.roll(vertices[:, 1], -1)) -
                            np.dot(vertices[:, 1], np.roll(vertices[:, 0], -1)))
    for (i, j), endpoints in zip(vor.ridge_points, vor.ridge_vertices):
        ii, jj = lookup[i], lookup[j]
        if ii < 0 and jj < 0:
            continue
        length = np.linalg.norm(vor.vertices[endpoints[0]]-vor.vertices[endpoints[1]])
        conductance = length/(2.*np.linalg.norm(xy[i]-xy[j]))
        if ii >= 0:
            diagonal[ii] += conductance
        if jj >= 0:
            diagonal[jj] += conductance
        if ii >= 0 and jj >= 0:
            rows.extend((ii, jj))
            cols.extend((jj, ii))
            values.extend((-conductance, -conductance))
    jac = coo_matrix((values, (rows, cols)), shape=(len(active), len(active))).tocsr()
    jac += diags(diagonal, format='csr')
    history = [dict(stage='unweighted', **metrics(areas, target, band))]
    warm_path = os.environ.get('LAGUERRE_WARM_START_PATH')
    if warm_path:
        with np.load(warm_path) as state:
            if not np.array_equal(state['active'], active):
                raise ValueError('Warm-start active set differs')
            weights[active] = state['weight']
        history.append(dict(stage='warm_start', source=warm_path))
    else:
        first, first_iters = solve(jac, diagonal, target-areas)
        weights[active] = first
        history.append(dict(stage='first_linear_step', cg_iterations=first_iters))
    for iteration in range(6):
        areas, meta, triangles, centers = power_areas(
            xy, weights, active, lookup, return_geometry=True)
        if meta['hidden_active'] or np.any(areas <= 0.):
            raise ValueError('Power cell disappeared or acquired invalid area')
        current = metrics(areas, target, band)
        history.append(dict(stage=f'nonlinear_{iteration}', **meta, **current))
        if current['max_rel'] < 1e-7:
            break
        jac, diagonal = newton_jacobian(xy, lookup, len(active), triangles, centers)
        correction, iterations = solve(jac, diagonal, target-areas)
        history[-1]['cg_iterations'] = iterations
        if iteration == 0:
            probe_scale = 1e-5
            probe = weights.copy()
            probe[active] += probe_scale*correction
            probe_areas, probe_meta = power_areas(xy, probe, active, lookup)
            if probe_meta['hidden_active'] or np.any(probe_areas <= 0.):
                raise ValueError('Finite-difference Jacobian probe invalid')
            linear_prediction = jac@correction
            measured = (probe_areas-areas)/probe_scale
            history[-1]['jacobian_fd_relative_l2'] = float(
                np.linalg.norm(measured-linear_prediction)/np.linalg.norm(linear_prediction))
            history[-1]['jacobian_fd_norm_ratio'] = float(
                np.linalg.norm(measured)/np.linalg.norm(linear_prediction))
            history[-1]['jacobian_fd_cosine'] = float(
                np.dot(measured, linear_prediction)/
                (np.linalg.norm(measured)*np.linalg.norm(linear_prediction)))
            history[-1]['jacobian_fd_inner_relative_l2'] = float(
                np.linalg.norm((measured-linear_prediction)[band])/
                np.linalg.norm(linear_prediction[band]))
            if os.environ.get('LAGUERRE_DEBUG_ONE'):
                break
        accepted = False
        for factor in (1., .5, .25, .125, .0625):
            candidate = weights.copy()
            candidate[active] += factor*correction
            candidate_areas, candidate_meta = power_areas(xy, candidate, active, lookup)
            if candidate_meta['hidden_active'] or np.any(candidate_areas <= 0.):
                continue
            candidate_metrics = metrics(candidate_areas, target, band)
            if candidate_metrics['relative_l2'] < current['relative_l2']:
                weights = candidate
                history[-1]['accepted_factor'] = factor
                accepted = True
                break
        if not accepted:
            history[-1]['accepted_factor'] = None
            break
    final_areas, _ = power_areas(xy, weights, active, lookup)
    final_relative = abs(final_areas-target)/target
    worst = np.argsort(final_relative)[-12:][::-1]
    result = dict(n_sites=len(xy), n_active=len(active), time=target_time,
                  geometry_mode=geometry_mode, source=source,
                  boundary_weights=f'zero outside {inner_radius:g} <= initial r < 8',
                  target='initial unweighted Voronoi area at each active site',
                  history=history, total_wall_seconds=time.monotonic()-started,
                  worst_cells=[dict(site=int(active[k]), initial_r=float(r0[active[k]]),
                                    relative_volume_error=float(final_relative[k]),
                                    target_area=float(target[k])) for k in worst],
                  limitation='Geometric t=6 existence test only; no time evolution, weighted-face conservation or near-AREPO runtime evidence')
    radius_label = 'R1P9' if inner_radius == 1.9 else f'R{int(inner_radius)}'
    suffix = f'{radius_label}_{time_label}'
    if geometry_mode == 'actual':
        suffix += '_ACTUAL'
    if os.environ.get('LAGUERRE_DEBUG_ONE'):
        suffix += '_DEBUG'
    out = Path(__file__).with_name(f'LAGUERRE_NEWTON_VOLUME_{suffix}_2026-10-26.json')
    with out.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    state_path = Path(__file__).with_name(f'LAGUERRE_NEWTON_VOLUME_{suffix}_2026-10-26.npz')
    with state_path.open('xb') as stream:
        np.savez_compressed(stream, active=active, weight=weights[active])
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
