"""One-time power-diagram volume probe at exact ballistic Kepler t=6.

This is a geometry diagnostic. It does not evolve the hydro scheme or establish
that time-dependent weights can be used conservatively and efficiently.
"""

from pathlib import Path
import json
import sys
import time

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import LinearOperator, cg
from scipy.spatial import ConvexHull, Voronoi

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
sys.path.insert(0, str(ROOT / 'rk4_diagnostics/varpress_w0_geom_cfl_long_20261026'))
from audit_varpress_checkpoint_20261026 import load_complete


def power_areas(xy, weights, active, lookup, return_geometry=False):
    """Areas of bounded active power cells via the regular triangulation."""
    q = np.einsum('ij,ij->i', xy, xy) - weights
    lifted = np.column_stack((xy, q))
    hull = ConvexHull(lifted, qhull_options='Qx')
    triangles = hull.simplices[hull.equations[:, 2] < -1e-12]
    a, b, c = (xy[triangles[:, k]] for k in range(3))
    rhs_b = q[triangles[:, 1]] - q[triangles[:, 0]]
    rhs_c = q[triangles[:, 2]] - q[triangles[:, 0]]
    ab = b - a
    ac = c - a
    det = ab[:, 0]*ac[:, 1] - ab[:, 1]*ac[:, 0]
    good = abs(det) > 1e-13
    if not np.all(good):
        triangles = triangles[good]
        ab, ac, rhs_b, rhs_c, det = (v[good] for v in
                                     (ab, ac, rhs_b, rhs_c, det))
    centers = np.column_stack(((rhs_b*ac[:, 1] - rhs_c*ab[:, 1])/(2*det),
                               (ab[:, 0]*rhs_c - ac[:, 0]*rhs_b)/(2*det)))
    incident = [[] for _ in active]
    for triangle_id, triple in enumerate(triangles):
        for site in triple:
            local_id = lookup[site]
            if local_id >= 0:
                incident[local_id].append(triangle_id)
    areas = np.full(len(active), np.nan)
    for local_id, ids in enumerate(incident):
        if len(ids) < 3:
            continue
        vertices = centers[ids]
        angles = np.arctan2(vertices[:, 1]-xy[active[local_id], 1],
                            vertices[:, 0]-xy[active[local_id], 0])
        vertices = vertices[np.argsort(angles)]
        areas[local_id] = .5*abs(np.dot(vertices[:, 0], np.roll(vertices[:, 1], -1)) -
                                np.dot(vertices[:, 1], np.roll(vertices[:, 0], -1)))
    metadata = dict(lower_facets=len(triangles),
                    hidden_active=int(np.count_nonzero(~np.isfinite(areas))))
    if return_geometry:
        return areas, metadata, triangles, centers
    return areas, metadata


def main():
    clock = time.monotonic()
    p0, meta = load_complete(ROOT / 'core_refinement/fine_integrated/khout.000001.dat')
    if abs(float(meta['t'])) > 1e-12:
        raise ValueError('Expected t=0 snapshot')
    xy0 = np.column_stack((p0['x'].astype(float)-12., p0['y'].astype(float)-12.))
    r0 = np.linalg.norm(xy0, axis=1)
    theta = 6.*(r0*r0 + .047**2)**(-.75)
    xy = np.column_stack((xy0[:, 0]*np.cos(theta)-xy0[:, 1]*np.sin(theta),
                          xy0[:, 0]*np.sin(theta)+xy0[:, 1]*np.cos(theta)))
    active = np.flatnonzero((r0 >= 1.) & (r0 < 8.))
    lookup = np.full(len(xy), -1, dtype=int)
    lookup[active] = np.arange(len(active))
    vor = Voronoi(xy, qhull_options='Qbb Qc Qz')
    volume = np.empty(len(active))
    for k, site in enumerate(active):
        ids = vor.regions[vor.point_region[site]]
        if not ids or -1 in ids:
            raise ValueError(f'Unbounded active Voronoi cell {site}')
        v = vor.vertices[ids]
        volume[k] = .5*abs(np.dot(v[:, 0], np.roll(v[:, 1], -1)) -
                             np.dot(v[:, 1], np.roll(v[:, 0], -1)))
    target = p0['volume'][active].astype(float)
    rhs = target-volume
    diagonal = np.zeros(len(active))
    rows, cols, values = [], [], []
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
    jacobian = coo_matrix((values, (rows, cols)), shape=(len(active), len(active))).tocsr()
    jacobian += diags(diagonal, format='csr')
    weight, info = cg(jacobian, rhs, M=LinearOperator(jacobian.shape,
                      matvec=lambda v:v/diagonal), rtol=1e-10, atol=0., maxiter=2000)
    if info != 0:
        raise RuntimeError(f'PCG info={info}')
    all_weights = np.zeros(len(xy))
    checks = []
    for scale in (0., .25, .5, 1.):
        all_weights[active] = scale*weight
        begin = time.monotonic()
        areas, mesh = power_areas(xy, all_weights, active, lookup)
        valid = np.isfinite(areas) & (areas > 0.)
        rel = (areas[valid]-target[valid])/target[valid]
        band = valid & (r0[active] >= 2.) & (r0[active] < 2.25)
        check = dict(scale=scale, wall_seconds=time.monotonic()-begin, **mesh,
                     finite_positive_fraction=float(np.mean(valid)),
                     median_abs_target_relative=float(np.median(abs(rel))),
                     p95_abs_target_relative=float(np.quantile(abs(rel), .95)),
                     inner_r2_2p25_median_abs_target_relative=float(
                         np.median(abs(areas[band]-target[band])/target[band])))
        if scale == 0.:
            check['baseline_vs_scipy_voronoi_median_rel'] = float(
                np.median(abs(areas[valid]-volume[valid])/volume[valid]))
            check['baseline_vs_scipy_voronoi_p99_rel'] = float(
                np.quantile(abs(areas[valid]-volume[valid])/volume[valid], .99))
        checks.append(check)
    result = dict(n_sites=len(xy), n_active=len(active),
                  boundary_weights='zero for r<1 and r>=8',
                  first_order_step='PCG of exact unweighted Voronoi volume Jacobian',
                  checks=checks, total_wall_seconds=time.monotonic()-clock,
                  limitation='Frozen t=6 geometry; no nonlinear Newton iteration, hydro evolution, conservative weighted-face proof, or GPU speed test')
    output = Path(__file__).with_name('LAGUERRE_NONLINEAR_VOLUME_PROBE_T6_2026-10-26.json')
    with output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
