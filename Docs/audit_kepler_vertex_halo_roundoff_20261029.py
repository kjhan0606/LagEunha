"""Localize the stored halo20 failure; exact rational references are diagnostic only."""
import hashlib
import json
from pathlib import Path
from fractions import Fraction as F

import numpy as np

from audit_kepler_geometry_force_dump_20261029 import read_dump

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/vertex_halo_pair_t2654918_ada_20261029')


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def analyze(call):
    paths = [ROOT/case/f'geom_force_rank000_call{call:03d}.bin' for case in ['control40', 'vertex20']]
    a, b = map(read_dump, paths)
    if not (np.array_equal(a['offset'], b['offset']) and
            np.array_equal(a['id'][:a['nbp']], b['id'][:b['nbp']])):
        raise ValueError('Unexpected topology/identity mismatch')
    fields, affected = {}, set()
    for name in ['c1x', 'c1y', 'c2x', 'c2y']:
        x, y = np.asarray(a[name]), np.asarray(b[name])
        indices = np.flatnonzero(x.view(np.uint64) != y.view(np.uint64))
        owners = np.searchsorted(a['offset'], indices, side='right')-1
        affected.update(map(int, a['id'][owners]))
        fields[name] = {'bit_differences': len(indices),
                        'numeric_differences': int(np.count_nonzero(x != y)),
                        'signed_zero_only': int(np.count_nonzero((x[indices] == 0) & (y[indices] == 0))),
                        'max_abs_difference': float(np.max(abs(x-y)))}
    if affected != {262}:
        raise ValueError('This bounded diagnostic expects the observed one-cell failure')
    i = int(np.flatnonzero(a['id'][:a['nbp']] == 262)[0])
    lo, hi = map(int, a['offset'][i:i+2])
    px, py = F(float(a['x'][i])), F(float(a['y'][i]))
    planes = []
    rounded_direction_errors = []
    for face in range(lo, hi):
        j = int(a['neighbor'][face])
        if j < 0:
            raise ValueError('Unexpected missing support plane')
        rx, ry = F(float(a['x'][j]))-px, F(float(a['y'][j]))-py
        planes.append((rx, ry, (rx*rx+ry*ry)/2))
        rounded_direction_errors.append(max(
            float(abs(F(float(a['x'][j]-a['x'][i]))-rx)),
            float(abs(F(float(a['y'][j]-a['y'][i]))-ry))))
    vertices = []
    for k, face in enumerate(range(lo, hi)):
        ax, ay, c = planes[k-1]
        bx, by, d = planes[k]
        det = ax*by-ay*bx
        if not det:
            vertices.append({'face': face, 'singular_adjacent_planes': True})
            continue
        exact = ((c*by-ay*d)/det, (ax*d-c*bx)/det)
        norm_a = max(abs(ax)+abs(ay), abs(bx)+abs(by))
        norm_inv = max(abs(by)+abs(ay), abs(bx)+abs(ax))/abs(det)
        other_slack = [z-x*exact[0]-y*exact[1] for n, (x, y, z) in enumerate(planes)
                       if n not in ((k-1) % len(planes), k)]
        row = {'face': face, 'singular_adjacent_planes': False,
               'reference_xy': [float(x) for x in exact],
               'abs_determinant': float(abs(det)), 'condition_inf': float(norm_a*norm_inv),
               'minimum_other_recorded_support_slack': float(min(other_slack))}
        for label, data in [('halo40', a), ('halo20', b)]:
            v = [F(float(data['c1x'][face])), F(float(data['c1y'][face]))]
            row[label] = {'max_abs_exact_intersection_error': max(float(abs(v[n]-exact[n])) for n in range(2)),
                          'max_abs_exact_adjacent_support_residual': max(
                              float(abs(x*v[0]+y*v[1]-z)) for x, y, z in [planes[k-1], planes[k]])}
        vertices.append(row)
    return {'force_call': call, 'affected_owned_IDs': sorted(affected),
            'owned_cells': a['nbp'], 'owner_xy': [float(px), float(py)], 'owner_faces': hi-lo,
            'coordinate_differences': fields, 'ideal_binary_position_plane_reference': vertices,
            'max_native_direction_subtraction_error': max(rounded_direction_errors),
            'dump_sha256': {str(path): sha(path) for path in paths}}


if __name__ == '__main__':
    prior = json.loads((ROOT/'comparison.json').read_text())
    result = {'original_strict_seven_case_passed': prior['passed'],
              'force_calls': [analyze(i) for i in [2, 4]],
              'reference_method': 'Exact Fraction arithmetic from the represented generator positions and '
                                  'recorded periodic supporting images, using adjacent W0 bisectors. '
                                  'Only the reported summaries are converted to binary64.',
              'native_rounded_plane_coefficients_serialized': False,
              'candidate_clipping_order_serialized': False,
              'independent_native_operation_path_roundoff_bound_established': False,
              'all_unrecorded_generator_constraints_checked': False,
              'new_geometry_equivalence_gate_accepted': False,
              'scope': 'Failure localization and reference errors, not proof of clipping-order causality '
                       'or native rounded-plane equivalence. Original FAIL remains; halo20 is not promoted.',
              'auditor_sha256': sha(Path(__file__))}
    output = Path(__file__).parent/'KEPLER_VERTEX_HALO_ADA_ROUNDOFF_2026-10-29.json'
    with output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'output': str(output), 'original_strict_seven_case_passed': result['original_strict_seven_case_passed'],
                      'affected_IDs': [x['affected_owned_IDs'] for x in result['force_calls']],
                      'new_geometry_equivalence_gate_accepted': False}, indent=2))
