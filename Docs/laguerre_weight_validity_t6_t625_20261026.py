"""Audit power-face weight differences against the native |dW|/d^2 guard."""

from pathlib import Path
import json
import os
import sys

import numpy as np

from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
sys.path.insert(0, str(ROOT / 'rk4_diagnostics/varpress_w0_geom_cfl_long_20261026'))
from audit_varpress_checkpoint_20261026 import load_complete


def main():
    geometry_mode = os.environ.get('LAGUERRE_GEOMETRY', 'ballistic')
    if geometry_mode not in ('ballistic', 'actual'):
        raise ValueError('LAGUERRE_GEOMETRY must be ballistic or actual')
    p0, _ = load_complete(ROOT / 'core_refinement/fine_integrated/khout.000001.dat')
    xy0 = np.column_stack((p0['x'].astype(float)-12., p0['y'].astype(float)-12.))
    r0 = np.linalg.norm(xy0, axis=1)
    out = {}
    for label, moment in (('T6', 6.), ('T625', 6.250024729493809)):
        geometry_suffix = '_ACTUAL' if geometry_mode == 'actual' else ''
        state_path = Path(__file__).with_name(
            f'LAGUERRE_NEWTON_VOLUME_R2_{label}{geometry_suffix}_2026-10-26.npz')
        with np.load(state_path) as state:
            active, local_weight = state['active'], state['weight']
        lookup = np.full(len(xy0), -1, dtype=int)
        lookup[active] = np.arange(len(active))
        weights = np.zeros(len(xy0))
        weights[active] = local_weight
        if geometry_mode == 'ballistic':
            theta = moment*(r0*r0+.047**2)**(-.75)
            xy = np.column_stack((xy0[:, 0]*np.cos(theta)-xy0[:, 1]*np.sin(theta),
                                  xy0[:, 0]*np.sin(theta)+xy0[:, 1]*np.cos(theta)))
        else:
            suffix_path = ('C010/to6/khout.032479.dat' if label == 'T6' else
                           'C010/to10/khout.036018.dat')
            snapshot, meta = load_complete(
                ROOT / 'rk4_diagnostics/varpress_w0_geom_cfl_long_20261026' / suffix_path)
            if not np.array_equal(snapshot['u4if'], p0['u4if']):
                raise ValueError('Snapshot particle IDs differ from t=0')
            if abs(float(meta['t'])-moment) > 1e-12:
                raise ValueError('Unexpected snapshot time')
            xy = np.column_stack((snapshot['x'].astype(float)-12.,
                                  snapshot['y'].astype(float)-12.))
        _, geometry, triangles, _ = power_areas(xy, weights, active, lookup,
                                                return_geometry=True)
        faces = np.vstack((triangles[:, [0, 1]], triangles[:, [1, 2]],
                           triangles[:, [2, 0]]))
        faces.sort(axis=1)
        faces = np.unique(faces, axis=0)
        touches_active = ((lookup[faces[:, 0]] >= 0) |
                          (lookup[faces[:, 1]] >= 0))
        faces = faces[touches_active]
        distances_squared = np.sum((xy[faces[:, 0]]-xy[faces[:, 1]])**2, axis=1)
        ratios = abs(weights[faces[:, 0]]-weights[faces[:, 1]])/distances_squared
        inner = (((r0[faces[:, 0]] >= 2.) & (r0[faces[:, 0]] < 2.25)) |
                 ((r0[faces[:, 1]] >= 2.) & (r0[faces[:, 1]] < 2.25)))
        flagged = np.flatnonzero(ratios > .8)
        out[label] = dict(n_power_faces_touching_active=len(faces),
                          n_active=int(len(active)),
                          hidden_active=geometry['hidden_active'],
                          quantiles=np.quantile(ratios, [.5, .95, .99, .999, 1]).tolist(),
                          n_over_0p8=int(len(flagged)),
                          n_over_1=int(np.count_nonzero(ratios > 1.)),
                          inner_faces=int(np.count_nonzero(inner)),
                          inner_n_over_0p8=int(np.count_nonzero(inner & (ratios > .8))),
                          flagged_faces=[dict(site_i=int(faces[k, 0]),
                                              site_j=int(faces[k, 1]),
                                              initial_r_i=float(r0[faces[k, 0]]),
                                              initial_r_j=float(r0[faces[k, 1]]),
                                              weight_difference_over_d2=float(ratios[k]))
                                         for k in flagged])
    suffix = '_ACTUAL' if geometry_mode == 'actual' else ''
    output = Path(__file__).with_name(f'LAGUERRE_WEIGHT_VALIDITY_FACES{suffix}_2026-10-26.json')
    with output.open('x') as stream:
        json.dump(out, stream, indent=2)
        stream.write('\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
