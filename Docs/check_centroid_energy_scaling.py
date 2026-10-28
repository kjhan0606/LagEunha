"""Screen a conservative centroid penalty; not a hydrodynamic evolution.

Periodic 8x8 disordered motif, repeated under refinement. Unit rho/P at
zero shear, fixed cell masses and entropy during each shear perturbation.
E_D = sum U0 * |centroid-generator|^2 / V (beta=1, two dimensions).
"""
import json
import numpy as np
from scipy.spatial import Voronoi

N = 8
GAMMA = 5 / 3
rng = np.random.default_rng(20261018)
points = (np.array([(i, j) for i in range(N) for j in range(N)])
          + .5 + rng.uniform(-.2, .2, (N*N, 2))) / N
shifts = np.array([(i, j) for i in range(-2, 3) for j in range(-2, 3)])
center = int(np.where(np.all(shifts == 0, axis=1))[0][0])


def geometry(scale, shear):
    matrix = np.array([[1., shear], [0., 1.]])
    cloud = ((points[None, :, :] + shifts[:, None, :]) @ matrix.T * scale).reshape(-1, 2)
    mesh = Voronoi(cloud)
    volumes, shape = [], []
    for a in range(center*N*N, (center+1)*N*N):
        region = mesh.regions[mesh.point_region[a]]
        assert region and min(region) >= 0
        q = mesh.vertices[region] - cloud[a]
        b = np.roll(q, -1, axis=0)
        cross = q[:, 0]*b[:, 1] - q[:, 1]*b[:, 0]
        area = cross.sum()/2
        centroid = ((q+b)*cross[:, None]).sum(axis=0)/(6*area)
        volumes.append(abs(area))
        shape.append(centroid @ centroid / abs(area))
    volumes = np.array(volumes)
    assert abs(volumes.sum()/scale**2-1) < 1e-11
    return volumes, np.array(shape)


rows = []
for scale in [1., .5, .25]:
    v0, _ = geometry(scale, 0.)
    def energy(shear):
        v, shape = geometry(scale, shear)
        thermal = v0/(GAMMA-1)*(v0/v)**(GAMMA-1)
        return float(np.dot(thermal, shape)/scale**2)
    derivatives = []
    for eps in [1e-3, 5e-4, 2.5e-4]:
        derivatives.append({'epsilon': eps,
                            'shear_stress': (energy(eps)-energy(-eps))/(2*eps),
                            'shear_curvature': (energy(eps)+energy(-eps)-2*energy(0.))/eps**2})
    rows.append({'scale': scale, 'added_energy_density': energy(0.),
                 'derivatives': derivatives})
print(json.dumps({'scope': 'Static thermodynamic fixed-entropy shear derivative; not Kepler, force validation, or a universal impossibility proof.',
                  'family': 'Same relative disorder under uniform refinement, tiled to a fixed domain.',
                  'candidate': 'beta=1 centroid-distance-squared / cell-area multiplier of thermal energy',
                  'rows': rows}, indent=2))
