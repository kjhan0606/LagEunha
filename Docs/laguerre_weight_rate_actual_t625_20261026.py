"""Reference dW/dt probe for fixed-volume Laguerre cells on native C010 t=6.25.

The generator velocity is taken from the native checkpoint. This geometry-only
probe does not update conserved variables or the hydrodynamic face flux.
"""

from pathlib import Path
import json
import sys
import time

import numpy as np
from scipy.sparse.linalg import LinearOperator, cg

from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from laguerre_newton_volume_t6_20261026 import newton_jacobian

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
RUN = ROOT / 'rk4_diagnostics/varpress_w0_geom_cfl_long_20261026'
sys.path.insert(0, str(RUN))
from audit_varpress_checkpoint_20261026 import load_complete


def relative_metrics(area, reference, band):
    relative = abs(area-reference)/reference
    return dict(median=float(np.median(relative)),
                p95=float(np.quantile(relative, .95)),
                max=float(relative.max()),
                inner_median=float(np.median(relative[band])),
                inner_p95=float(np.quantile(relative[band], .95)),
                inner_max=float(relative[band].max()))


def face_volume_rate(xy, velocity, weights, weight_rate, active, lookup,
                     triangles, centers):
    """Exact geometric face-normal speed from the differentiated power bisector."""
    pairs = np.vstack((triangles[:, [0, 1]], triangles[:, [1, 2]],
                       triangles[:, [2, 0]]))
    ids = np.tile(np.arange(len(triangles)), 3)
    pairs.sort(axis=1)
    order = np.lexsort((pairs[:, 1], pairs[:, 0]))
    pairs, ids = pairs[order], ids[order]
    starts = np.r_[0, 1+np.flatnonzero(np.any(pairs[1:] != pairs[:-1], axis=1))]
    counts = np.diff(np.r_[starts, len(pairs)])
    touching = ((lookup[pairs[starts, 0]] >= 0) |
                (lookup[pairs[starts, 1]] >= 0))
    if np.any(counts[touching] != 2):
        raise ValueError('Non-finite face touching active disk')
    starts = starts[touching]
    i, j = pairs[starts, 0], pairs[starts, 1]
    face_midpoint = .5*(centers[ids[starts]]+centers[ids[starts+1]])
    face_length = np.linalg.norm(centers[ids[starts]]-centers[ids[starts+1]],
                                 axis=1)
    distance = np.linalg.norm(xy[j]-xy[i], axis=1)
    fixed_normal_speed = (
        np.einsum('ij,ij->i', xy[j], velocity[j])-
        np.einsum('ij,ij->i', xy[i], velocity[i])-
        np.einsum('ij,ij->i', velocity[j]-velocity[i], face_midpoint))/distance
    if np.any(~np.isfinite(fixed_normal_speed)):
        raise ValueError('Nonfinite face-normal velocity')
    fixed = np.zeros(len(active))
    weighted = np.zeros(len(active))
    normal_speed = fixed_normal_speed + (weight_rate[i]-weight_rate[j])/(2*distance)
    ii, jj = lookup[i], lookup[j]
    ai, aj = ii >= 0, jj >= 0
    np.add.at(fixed, ii[ai], face_length[ai]*fixed_normal_speed[ai])
    np.add.at(fixed, jj[aj], -face_length[aj]*fixed_normal_speed[aj])
    np.add.at(weighted, ii[ai], face_length[ai]*normal_speed[ai])
    np.add.at(weighted, jj[aj], -face_length[aj]*normal_speed[aj])
    return fixed, weighted, len(starts)


def main():
    started = time.monotonic()
    p0, meta0 = load_complete(ROOT / 'core_refinement/fine_integrated/khout.000001.dat')
    current, meta = load_complete(RUN / 'C010/to10/khout.036018.dat')
    if abs(float(meta0['t'])) > 1e-12 or abs(float(meta['t'])-6.250024729493809) > 1e-12:
        raise ValueError('Wrong checkpoint time')
    if not np.array_equal(p0['u4if'], current['u4if']):
        raise ValueError('Checkpoint particle IDs differ')
    xy0 = np.column_stack((p0['x'].astype(float)-12., p0['y'].astype(float)-12.))
    r0 = np.linalg.norm(xy0, axis=1)
    xy = np.column_stack((current['x'].astype(float)-12.,
                          current['y'].astype(float)-12.))
    velocity = np.column_stack((current['vx'].astype(float),
                                current['vy'].astype(float)))
    with np.load(Path(__file__).with_name(
            'LAGUERRE_NEWTON_VOLUME_R2_T625_ACTUAL_2026-10-26.npz')) as state:
        active, local_weight = state['active'], state['weight']
    lookup = np.full(len(xy), -1, dtype=int)
    lookup[active] = np.arange(len(active))
    weight = np.zeros(len(xy))
    weight[active] = local_weight
    band = (r0[active] >= 2.) & (r0[active] < 2.25)
    area, geometry, triangles, centers = power_areas(
        xy, weight, active, lookup, return_geometry=True)
    if geometry['hidden_active'] or np.any(area <= 0.):
        raise ValueError('Invalid initial power diagram')
    jac, diagonal = newton_jacobian(xy, lookup, len(active), triangles, centers)
    h = 1e-4
    area_plus, _ = power_areas(xy+h*velocity, weight, active, lookup)
    area_minus, _ = power_areas(xy-h*velocity, weight, active, lookup)
    geometric_rate = (area_plus-area_minus)/(2*h)
    right_hand_side = -geometric_rate
    rhs_norm = np.linalg.norm(right_hand_side)
    snapshots = {}
    iterations = [0]

    def capture(vector):
        iterations[0] += 1
        if iterations[0] in (8, 16, 32, 64, 128, 256):
            snapshots[iterations[0]] = vector.copy()

    preconditioner = LinearOperator(jac.shape, matvec=lambda vector:vector/diagonal)
    rate, info = cg(jac, right_hand_side, M=preconditioner, rtol=1e-10,
                    atol=0., maxiter=2000, callback=capture)
    if info != 0:
        raise RuntimeError(f'Cold CG failed: info={info}')
    residuals = {str(k):float(np.linalg.norm(right_hand_side-jac@vector)/rhs_norm)
                 for k, vector in snapshots.items()}
    residuals[str(iterations[0])] = float(np.linalg.norm(right_hand_side-jac@rate)/rhs_norm)
    with np.load(Path(__file__).with_name(
            'LAGUERRE_NEWTON_VOLUME_R2_T6_ACTUAL_2026-10-26.npz')) as state:
        if not np.array_equal(state['active'], active):
            raise ValueError('t=6 active set differs')
        secant = (local_weight-state['weight'])/(6.250024729493809-6.)
    warm_iterations = [0]
    warm_rate, warm_info = cg(jac, right_hand_side, x0=secant,
                              M=preconditioner, rtol=1e-10, atol=0.,
                              maxiter=2000,
                              callback=lambda _:warm_iterations.__setitem__(0,
                                  warm_iterations[0]+1))
    if warm_info != 0:
        raise RuntimeError(f'Warm CG failed: info={warm_info}')
    if np.linalg.norm(warm_rate-rate) > 1e-6*np.linalg.norm(rate):
        raise ValueError('Cold and warm CG solutions disagree')
    all_rate = np.zeros(len(xy))
    all_rate[active] = rate
    fixed_face_rate, weighted_face_rate, nfaces = face_volume_rate(
        xy, velocity, weight, all_rate, active, lookup, triangles, centers)
    face_vs_finite_difference = float(
        np.linalg.norm(fixed_face_rate-geometric_rate)/np.linalg.norm(geometric_rate))
    face_rate_residual = float(
        np.linalg.norm(weighted_face_rate)/np.linalg.norm(geometric_rate))
    step = 1e-3
    next_xy = xy+step*velocity
    area_fixed, _ = power_areas(next_xy, weight, active, lookup)
    next_weight = weight.copy()
    next_weight[active] += step*rate
    area_rate, next_geometry = power_areas(next_xy, next_weight, active, lookup)
    if next_geometry['hidden_active'] or np.any(area_rate <= 0.):
        raise ValueError('Predicted diagram has invalid cells')
    truncated = {}
    for count in (16, 32, 64):
        candidate_weight = weight.copy()
        candidate_weight[active] += step*snapshots[count]
        candidate_area, candidate_geometry = power_areas(
            next_xy, candidate_weight, active, lookup)
        if candidate_geometry['hidden_active'] or np.any(candidate_area <= 0.):
            raise ValueError(f'Invalid {count}-iteration predicted diagram')
        truncated[str(count)] = relative_metrics(candidate_area, area, band)
    secant_weight = weight.copy()
    secant_weight[active] += step*secant
    secant_area, secant_geometry = power_areas(next_xy, secant_weight, active, lookup)
    if secant_geometry['hidden_active'] or np.any(secant_area <= 0.):
        raise ValueError('Invalid secant-predicted diagram')
    result = dict(t=float(meta['t']), n_active=len(active),
                  generator_displacement='native checkpoint velocity, forward Euler',
                  volume_derivative='centered finite difference at h=1e-4, fixed weights',
                  target='zero cell-volume derivative in the r>=2 disk',
                  graph_laplacian='current weighted power-face area/(2 generator distance)',
                  rhs_l2=float(rhs_norm), cg_iterations_cold=iterations[0],
                  cg_iterations_secant_warm=warm_iterations[0],
                  cg_relative_residual_by_iteration=residuals,
                  n_power_faces_touching_active=nfaces,
                  face_fixed_rate_vs_centered_fd_relative_l2=face_vs_finite_difference,
                  face_weighted_rate_over_fixed_rate_l2=face_rate_residual,
                  secant_rate_relative_l2=float(np.linalg.norm(secant-rate)/np.linalg.norm(rate)),
                  weight_rate_abs_quantiles=np.quantile(abs(rate), [.5, .95, .99, 1.]).tolist(),
                  step=step,
                  fixed_weight_area_change=relative_metrics(area_fixed, area, band),
                  truncated_cg_area_change=truncated,
                  secant_area_change=relative_metrics(secant_area, area, band),
                  predicted_weight_area_change=relative_metrics(area_rate, area, band),
                  wall_seconds=time.monotonic()-started,
                  limitation='Geometry-only, one Euler displacement; no hydro flux, conservation or production GPU performance result')
    output = Path(__file__).with_name('LAGUERRE_WEIGHT_RATE_ACTUAL_T625_V3_2026-10-26.json')
    with output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
