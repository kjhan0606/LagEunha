"""Necessary energy feasibility for a conservative local nodal-mass remap.

At fixed positions, prescribed positive nodal masses, P and L impose K>=Kmin.
If H-Phi_new-Pext*A-Kmin<=0, no velocity assignment can also retain U>0.
This checks feasibility only; it does not define an accepted remap algorithm.
"""
import json
import numpy as np
from feedback_patch import initial, mesh, geom, circular, cross, GAMMA
from stag_checks import overlap_remap


def lump(keys, mass):
    result = np.zeros(4)
    np.add.at(result, np.asarray(keys).ravel(), np.repeat(mass/3, 3))
    return result


def minimum_kinetic(x, mass, momentum, angular):
    total = mass.sum(); center = np.sum(mass[:, None]*x, axis=0)/total
    rel = x-center; inertia = np.sum(mass*np.sum(rel*rel, axis=1))
    spin = angular-cross(center, momentum)
    omega = spin/inertia
    velocity = momentum/total+omega*np.column_stack((-rel[:, 1], rel[:, 0]))
    kmin = np.dot(momentum, momentum)/(2*total)+spin*spin/(2*inertia)
    np.testing.assert_allclose(np.sum(mass[:, None]*velocity, axis=0), momentum, atol=1e-14)
    np.testing.assert_allclose(np.sum(mass*cross(x, velocity)), angular, atol=1e-14)
    np.testing.assert_allclose(.5*np.sum(mass*np.sum(velocity*velocity, axis=1)), kmin, rtol=1e-12)
    return kmin


def audit(h, density_ratio, p0, time=-.2, reverse=False):
    x = initial(h); old = mesh(x)
    radius = np.linalg.norm(x, axis=1)
    theta = np.arctan2(x[:, 1], x[:, 0])+(time+.2)/radius**1.5
    x = np.column_stack((radius*np.cos(theta), radius*np.sin(theta)))
    shared = set(old[0]) & set(old[1]); outer = list(set(range(4))-shared)
    new = [tuple(outer+[i]) for i in sorted(shared)]
    # Preserve a positive orientation for the alternate diagonal.
    for i, k in enumerate(new):
        a, b, c = x[list(k)]
        if cross(b-a, c-a) < 0: new[i] = k[::-1]
    area, _ = geom(x, old); rho = np.array([1., density_ratio])
    if reverse: rho = rho[::-1]
    mass = area*rho; u = p0*area/(GAMMA-1)
    mn, un, _ = overlap_remap(x[np.asarray(old)], x[np.asarray(new)], mass, u)
    before = lump(old, mass); after = lump(new, mn); v = circular(x)
    momentum = np.sum(before[:, None]*v, axis=0); angular = np.sum(before*cross(x, v))
    k0 = .5*np.sum(before*np.sum(v*v, axis=1)); radius = np.linalg.norm(x, axis=1)
    phi0 = -np.sum(before/radius); phi1 = -np.sum(after/radius)
    kmin = minimum_kinetic(x, after, momentum, angular)
    available_u = k0+u.sum()+phi0-phi1-kmin
    return dict(h=h, density_ratio=density_ratio, p0=p0, time=time, reverse=reverse,
                nodal_mass_l1=float(np.abs(after-before).sum()/before.sum()),
                original_internal=float(u.sum()), maximum_possible_internal=float(available_u),
                feasible_positive_internal=bool(available_u > 0),
                maximum_internal_excluding_gravity=float(k0+u.sum()-kmin),
                naive_momentum_change=float(np.linalg.norm(np.sum((after-before)[:, None]*v, axis=0))),
                naive_angular_change=float(np.sum((after-before)*cross(x, v))),
                naive_energy_change=float(np.sum((after-before)*(np.sum(v*v, axis=1)/2-1/radius))))


if __name__ == '__main__':
    for h in (.2, .1, .05):
        for rho in (1., .01):
            for p0 in (2.78e-5, .00278, .278):
                print(json.dumps(audit(h, rho, p0)), flush=True)
    # At t=0 the four nodes are a square, i.e. the exact cocircular flip event.
    for h in (.2, .1, .05):
        for rho in (1., .01):
            for reverse in (False, True):
                print(json.dumps(audit(h, rho, 2.78e-5, 0., reverse)), flush=True)
