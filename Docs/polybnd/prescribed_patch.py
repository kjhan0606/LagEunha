"""Open, prescribed four-node patch: geometry forcing, NOT a fluid solution.

Affine incompressible shear is an exact constant-pressure invariant. Circular
Kepler trajectories are non-affine: chord-triangle areas need not be constant.
No pressure feedback, gravity integration, taper, or frozen wall is involved.
"""
import json
import numpy as np
from scipy.spatial import Delaunay
from stag_checks import signed_area, check_orientation, overlap_remap


def mesh(x):
    return sorted(tuple(sorted(tri)) for tri in Delaunay(x).simplices)


def trajectory(t, kind, size):
    x = np.array([[.8-size/2, -size/2], [.8+size/2, -size/2],
                  [.8+size/2, size/2], [.8-size/2, size/2]])
    if kind == 'affine':
        x[:, 0] += t*x[:, 1]
    else:
        r = np.linalg.norm(x, axis=1)
        th = np.arctan2(x[:, 1], x[:, 0])+t/r**1.5
        x = np.column_stack((r*np.cos(th), r*np.sin(th)))
    return x


def run(kind, dt, density_ratio, size=.2, remap=True):
    gamma = 5/3; p0 = 2.78e-5
    x = trajectory(-.2, kind, size); keys = mesh(x)
    area = np.array([abs(signed_area(x[list(k)])) for k in keys])
    mass = area*np.array([1., density_ratio]); internal = p0*area/(gamma-1)
    entropy = (gamma-1)*internal/area/(mass/area)**gamma
    m0 = mass.sum(); flips = 0; noise = 0.; remapU = 0.; pressure_jump = 0.
    for t in np.linspace(-.2, .2, round(.4/dt)+1)[1:]:
        xn = trajectory(t, kind, size)
        for k in keys: check_orientation(xn[list(k)], signed_area(x[list(k)]))
        area = np.array([abs(signed_area(xn[list(k)])) for k in keys])
        internal = entropy*mass**gamma*area**(1-gamma)/(gamma-1)
        noise = max(noise, float(np.max(abs((gamma-1)*internal/area/p0-1))))
        newkeys = mesh(xn)
        if remap and newkeys != keys:
            old_u = internal.sum()
            mass, internal, _ = overlap_remap([xn[list(k)] for k in keys],
                                             [xn[list(k)] for k in newkeys], mass, internal)
            remapU += internal.sum()-old_u
            keys = newkeys; area = np.array([abs(signed_area(xn[list(k)])) for k in keys])
            pressure_jump = max(pressure_jump, float(np.max(abs((gamma-1)*internal/area/p0-1))))
            entropy = (gamma-1)*internal/area/(mass/area)**gamma
            flips += 1
        x = xn
    result = dict(kind=kind, dt=dt, density_ratio=density_ratio, size=size, remap=remap, flips=flips,
                  max_pressure_deviation=noise, post_flip_deviation=pressure_jump,
                  mass_relative=mass.sum()/m0-1, remap_internal_change=remapU)
    if kind == 'affine':
        if remap: assert flips > 0, 'control must exercise a topology change'
        assert noise < 1e-10, result
    assert abs(result['mass_relative']) < 1e-9, result
    return result


if __name__ == '__main__':
    for kind in ('affine', 'kepler'):
        for dt in (.01, .005):
            for ratio in (1., .01):
                print(json.dumps(run(kind, dt, ratio)), flush=True)
    for size in (.2, .1, .05):
        for dt in (.005, .0025):
            for remap in (False, True):
                print(json.dumps(run('kepler', dt, 1., size, remap)), flush=True)
