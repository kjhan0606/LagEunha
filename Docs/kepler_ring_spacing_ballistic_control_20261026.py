"""Static exact-orbit IC control for radial versus azimuthal seed spacing.

No hydro forces, remeshing adaptation, mass transfer or native GPU execution.
Compare the geometric fixed-mass density change of designed polar point sets.
"""

from pathlib import Path
import json
import time

import numpy as np
from scipy.spatial import Voronoi

H = 24./256.
EPS = .047
REFERENCE_RADIUS = 2.125
TIMES = (6.500079194680843, 10., 17.778891990936625)


def polar_points(aspect, variable_count):
    radial_spacing = H*np.sqrt(aspect)
    azimuthal_target = H/np.sqrt(aspect)
    rings = np.arange(1.25+.5*radial_spacing, 4.5, radial_spacing)
    xy = []
    radius = []
    ring_id = []
    for j, r in enumerate(rings):
        nphi = int(round(2*np.pi*(r if variable_count else REFERENCE_RADIUS)/azimuthal_target))
        angle = (np.arange(nphi)+.5*(j%2))*2*np.pi/nphi
        xy.append(np.column_stack((r*np.cos(angle), r*np.sin(angle))))
        radius.append(np.full(nphi, r))
        ring_id.append(np.full(nphi, j))
    return (np.vstack(xy), np.concatenate(radius), np.concatenate(ring_id),
            radial_spacing, azimuthal_target, len(rings))


def areas_at(points, selected):
    vor = Voronoi(points, qhull_options='Qbb Qc Qz')
    area = np.empty(len(selected))
    for k, site in enumerate(selected):
        region = vor.regions[vor.point_region[site]]
        if not region or -1 in region:
            raise ValueError(f'Unbounded target-ring cell {site}')
        polygon = vor.vertices[region]
        area[k] = .5*abs(np.dot(polygon[:, 0], np.roll(polygon[:, 1], -1))-
                      np.dot(polygon[:, 1], np.roll(polygon[:, 0], -1)))
    return area


def main():
    started = time.monotonic()
    cases = []
    for variable_count in (False, True):
        for aspect in (.5, 1., 2.):
            xy0, radius, ring, spacing_r, spacing_phi, nrings = polar_points(
                aspect, variable_count)
            selected = np.flatnonzero((radius >= 2.) & (radius < 2.25))
            area0 = areas_at(xy0, selected)
            mass = 2./radius[selected]*area0
            values = []
            for t in TIMES:
                angle = t*(radius*radius+EPS*EPS)**(-.75)
                xy = np.column_stack((xy0[:, 0]*np.cos(angle)-xy0[:, 1]*np.sin(angle),
                                      xy0[:, 0]*np.sin(angle)+xy0[:, 1]*np.cos(angle)))
                area = areas_at(xy, selected)
                relative_density = abs(area0/area-1.)
                values.append(dict(time=t,
                    density_change_mass_weighted=float(np.sum(mass*relative_density)/np.sum(mass)),
                    density_change_median=float(np.median(relative_density)),
                    density_change_p95=float(np.quantile(relative_density, .95)),
                    density_change_max=float(relative_density.max()),
                    min_area_ratio=float(np.min(area/area0))))
            cases.append(dict(variable_azimuthal_count=variable_count,
                              radial_over_azimuthal_spacing=aspect,
                              radial_spacing=spacing_r,
                              azimuthal_target_spacing=spacing_phi,
                              total_points=len(radius), target_ring_cells=len(selected),
                              radial_rings=nrings, target_radii=np.unique(radius[selected]).tolist(),
                              results=values))
    baseline = json.loads(Path(__file__).with_name(
        'BALLISTIC_VORONOI_CONTROL_T65_2026-10-26.json').read_text())
    native_bin = next(b for b in baseline['bins'] if b['rlo']==2 and b['rhi']==2.25)
    output = dict(description='Independent polar point sets at similar nominal area H^2; all generators follow exact softened Kepler orbits and keep initial cell masses',
                  box_reference_H=H, potential_softening=EPS,
                  target_ring='2 <= initial r < 2.25',
                  original_C010_exact_orbit_t65_mass_weighted_density_change=native_bin['ballistic_density_MW'],
                  cases=cases, wall_seconds=time.monotonic()-started,
                  limitation='Geometry-only IC comparison; ring symmetry, unequal target cell counts, boundaries, and later hydrodynamic instability are not controlled by this experiment')
    path = Path(__file__).with_name('KEPLER_RING_SPACING_BALLISTIC_CONTROL_2026-10-26.json')
    with path.open('x') as stream:
        json.dump(output, stream, indent=2)
        stream.write('\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
