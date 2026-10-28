"""Try simple radial taper of exact-volume disk weights at t=6.25.

This is a static screening study, not a conservative evolution method.
"""

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
    xy0 = np.column_stack((p0['x'].astype(float)-12., p0['y'].astype(float)-12.))
    r0 = np.linalg.norm(xy0, axis=1)
    t = 6.250024729493809
    theta = t*(r0*r0+.047**2)**(-.75)
    xy = np.column_stack((xy0[:, 0]*np.cos(theta)-xy0[:, 1]*np.sin(theta),
                          xy0[:, 0]*np.sin(theta)+xy0[:, 1]*np.cos(theta)))
    with np.load(Path(__file__).with_name('LAGUERRE_NEWTON_VOLUME_R2_T625_2026-10-26.npz')) as state:
        active, solved_weights = state['active'], state['weight']
    target = p0['volume'][active].astype(float)
    lookup = np.full(len(xy), -1, dtype=int)
    lookup[active] = np.arange(len(active))
    band = (r0[active] >= 2.) & (r0[active] < 2.25)
    results = []
    for width in (0., .025, .05, .1, .2, .3, .5):
        weights = np.zeros(len(xy))
        factor = (np.ones(len(active)) if width == 0. else
                  np.clip((r0[active]-2.)/width, 0., 1.))
        weights[active] = solved_weights*factor
        area, meta, triangles, _ = power_areas(xy, weights, active, lookup,
                                               return_geometry=True)
        faces = np.vstack((triangles[:, [0, 1]], triangles[:, [1, 2]],
                           triangles[:, [2, 0]]))
        faces.sort(axis=1)
        faces = np.unique(faces, axis=0)
        touches = (lookup[faces[:, 0]] >= 0) | (lookup[faces[:, 1]] >= 0)
        faces = faces[touches]
        d2 = np.sum((xy[faces[:, 0]]-xy[faces[:, 1]])**2, axis=1)
        ratio = abs(weights[faces[:, 0]]-weights[faces[:, 1]])/d2
        inner_volume = abs(area[band]-target[band])/target[band]
        inner_density = abs(target[band]/area[band]-1.)
        results.append(dict(taper_width=width, hidden_active=meta['hidden_active'],
                            nonpositive_active=int(np.count_nonzero(area <= 0.)),
                            n_over_0p8=int(np.count_nonzero(ratio > .8)),
                            n_over_1=int(np.count_nonzero(ratio > 1.)),
                            max_abs_dw_over_d2=float(ratio.max()),
                            inner_volume_median=float(np.median(inner_volume)),
                            inner_volume_p95=float(np.quantile(inner_volume, .95)),
                            inner_volume_max=float(inner_volume.max()),
                            inner_density_median=float(np.median(inner_density)),
                            inner_density_p95=float(np.quantile(inner_density, .95))))
    result = dict(moment=t, active_region='2 <= initial r < 8',
                  taper='w_i = exact-volume w_i * clip((r_i-2)/width,0,1)',
                  results=results,
                  limitation='Static ad hoc weight taper; no native face speed, mass or energy update')
    output = Path(__file__).with_name('LAGUERRE_BOUNDARY_TAPER_T625_2026-10-26.json')
    with output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
