"""Audit signed power weights, common gauge shift, and native containment limit."""

from pathlib import Path
import json
import sys

import numpy as np

from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
sys.path.insert(0, str(ROOT / 'rk4_diagnostics/varpress_w0_geom_cfl_long_20261026'))
from audit_varpress_checkpoint_20261026 import load_complete


def main():
    p0, _ = load_complete(ROOT / 'core_refinement/fine_integrated/khout.000001.dat')
    current, meta = load_complete(ROOT / 'rk4_diagnostics/varpress_w0_geom_cfl_long_20261026/C010/to10/khout.036018.dat')
    if abs(float(meta['t'])-6.250024729493809) > 1e-12:
        raise ValueError('Wrong checkpoint time')
    if not np.array_equal(p0['u4if'], current['u4if']):
        raise ValueError('Snapshot IDs changed')
    xy = np.column_stack((current['x'].astype(float)-12.,
                          current['y'].astype(float)-12.))
    r0 = np.hypot(p0['x'].astype(float)-12., p0['y'].astype(float)-12.)
    with np.load(Path(__file__).with_name('LAGUERRE_NEWTON_VOLUME_R2_T625_ACTUAL_2026-10-26.npz')) as state:
        active, local_weight = state['active'], state['weight']
    lookup = np.full(len(xy), -1, dtype=int)
    lookup[active] = np.arange(len(active))
    signed_weight = np.zeros(len(xy))
    signed_weight[active] = local_weight
    offset = -float(signed_weight.min())
    shifted_weight = signed_weight+offset
    signed_area, signed_geometry, triangles, _ = power_areas(
        xy, signed_weight, active, lookup, return_geometry=True)
    shifted_area, shifted_geometry = power_areas(xy, shifted_weight, active, lookup)
    _, _, voronoi_triangles, _ = power_areas(
        xy, np.zeros(len(xy)), active, lookup, return_geometry=True)
    if signed_geometry['hidden_active'] or shifted_geometry['hidden_active']:
        raise ValueError('Active cell disappeared under gauge shift')
    face_pairs = np.vstack((triangles[:, [0, 1]], triangles[:, [1, 2]],
                            triangles[:, [2, 0]]))
    face_pairs.sort(axis=1)
    face_pairs = np.unique(face_pairs, axis=0)
    voronoi_pairs = np.vstack((voronoi_triangles[:, [0, 1]],
                               voronoi_triangles[:, [1, 2]],
                               voronoi_triangles[:, [2, 0]]))
    voronoi_pairs.sort(axis=1)
    voronoi_pairs = np.unique(voronoi_pairs, axis=0)
    cell_index = np.floor((xy+12.)/(24./256.)).astype(int)
    def candidate_span_summary(pairs):
        touching = ((lookup[pairs[:, 0]] >= 0) |
                    (lookup[pairs[:, 1]] >= 0))
        pairs = pairs[touching]
        span = np.max(abs(cell_index[pairs[:, 0]]-cell_index[pairs[:, 1]]), axis=1)
        return dict(n_faces=len(pairs), max_cell_index_span=int(span.max()),
                    n_span_over_2=int(np.count_nonzero(span > 2)),
                    n_span_over_4=int(np.count_nonzero(span > 4)))
    i, j = face_pairs[:, 0], face_pairs[:, 1]
    distance2 = np.sum((xy[i]-xy[j])**2, axis=1)
    native_style_ceiling = np.full(len(xy), np.inf)
    np.minimum.at(native_style_ceiling, i, distance2+shifted_weight[j])
    np.minimum.at(native_style_ceiling, j, distance2+shifted_weight[i])
    violation = shifted_weight > native_style_ceiling+1e-12
    offenders = np.flatnonzero(violation)
    result = dict(moment=float(meta['t']),
                  n_sites=len(xy), n_active=len(active),
                  active_signed_weight_min=float(local_weight.min()),
                  active_signed_weight_max=float(local_weight.max()),
                  negative_active_fraction=float(np.mean(local_weight < 0.)),
                  common_offset=offset,
                  shifted_weight_min=float(shifted_weight.min()),
                  shifted_weight_max=float(shifted_weight.max()),
                  gauge_area_max_relative_difference=float(
                      np.max(abs(signed_area-shifted_area)/signed_area)),
                  candidate_span_voronoi=candidate_span_summary(voronoi_pairs),
                  candidate_span_laguerre=candidate_span_summary(face_pairs),
                  n_power_faces=len(face_pairs),
                  n_sites_above_native_style_power_ceiling=int(len(offenders)),
                  offending_sites=[dict(site=int(k), initial_r=float(r0[k]),
                                        shifted_weight=float(shifted_weight[k]),
                                        face_ceiling=float(native_style_ceiling[k]),
                                        excess=float(shifted_weight[k]-native_style_ceiling[k]))
                                   for k in offenders],
                  limitation='Power-face neighbor ceiling from independent regular triangulation; not a native CUDA execution or face-flux test')
    output = Path(__file__).with_name('LAGUERRE_GAUGE_CONTRACT_ACTUAL_T625_V2_2026-10-26.json')
    with output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
