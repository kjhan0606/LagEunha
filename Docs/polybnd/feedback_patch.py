"""Four-node pressure response experiment, not a global Kepler stability test.

Material mass and nodal inertia are held fixed across each +/- seed pair.
Uniform external pressure P0 acts on the material polygon: H=KE+U+Phi+P0*A.
feedback=0 disables the net pressure force and books the compensating work.
No taper, frozen wall, mass source, dissipation, or pressure floor is used.
"""
import json
import numpy as np
from scipy.spatial import Delaunay
from stag_checks import overlap_remap

GAMMA = 5/3
OMEGA = .8**-1.5


def cross(a, b):
    return a[..., 0]*b[..., 1]-a[..., 1]*b[..., 0]


def mesh(x):
    keys = []
    for tri in Delaunay(x).simplices:
        a, b, c = x[tri]
        if cross(b-a, c-a) < 0: tri = tri[::-1]
        tri = np.roll(tri, -np.argmin(tri))
        keys.append(tuple(int(i) for i in tri))
    return sorted(keys)


def geom(x, keys):
    ids = np.asarray(keys); p = x[ids]; a, b, c = p[:, 0], p[:, 1], p[:, 2]
    area = cross(b-a, c-a)/2
    if not np.isfinite(area).all() or (area <= 0).any():
        raise RuntimeError('inverted or nonfinite patch')
    differences = np.stack((c-b, a-c, b-a), axis=1)
    gradient = np.stack((-differences[..., 1], differences[..., 0]), axis=-1)/2
    return area, gradient


def pressure_force(x, keys, mass, entropy, p0):
    area, gradient = geom(x, keys)
    pressure = entropy*(mass/area)**GAMMA
    force = np.zeros_like(x)
    np.add.at(force, np.asarray(keys).ravel(), ((pressure-p0)[:, None, None]*gradient).reshape(-1, 2))
    return force, pressure, area


def circular(x):
    r = np.linalg.norm(x, axis=1)
    return np.column_stack((-x[:, 1], x[:, 0]))/r[:, None]**1.5


def initial(h):
    x = np.array([[.8-h/2, -h/2], [.8+h/2, -h/2], [.8+h/2, h/2], [.8-h/2, h/2]])
    r = np.linalg.norm(x, axis=1); th = np.arctan2(x[:, 1], x[:, 0])-.2/r**1.5
    return np.column_stack((r*np.cos(th), r*np.sin(th)))


def run(h, p0, dt, feedback, remap, seed=0., duration=.4):
    xb = initial(h); keys = mesh(xb); mass, _ = geom(xb, keys)  # rho=1
    inertia = np.zeros(4); np.add.at(inertia, np.asarray(keys).ravel(), np.repeat(mass/3, 3))
    tangent = np.column_stack((-xb[:, 1], xb[:, 0]))/np.linalg.norm(xb, axis=1)[:, None]
    x = xb + seed*h*np.array([1., -1., 1., -1.])[:, None]*tangent
    v = circular(x); xinitial = x.copy(); vinitial = v.copy()
    area, _ = geom(x, keys); entropy = p0/(mass/area)**GAMMA
    y = np.r_[x.ravel(), v.ravel(), 0.]

    def energy(y):
        x = y[:8].reshape(4, 2); v = y[8:16].reshape(4, 2)
        area, _ = geom(x, keys)
        u = entropy*mass**GAMMA*area**(1-GAMMA)/(GAMMA-1)
        return np.sum(inertia*((v*v).sum(1)/2-1/np.linalg.norm(x, axis=1)))+u.sum()+p0*area.sum()

    def rhs(y):
        x = y[:8].reshape(4, 2); v = y[8:16].reshape(4, 2)
        force, _, _ = pressure_force(x, keys, mass, entropy, p0)
        acc = -x/np.linalg.norm(x, axis=1)[:, None]**3 + feedback*force/inertia[:, None]
        work = -(1-feedback)*np.sum(force*v)
        return np.r_[v.ravel(), acc.ravel(), work]

    e0 = energy(y); thermal0 = p0*area.sum()/(GAMMA-1)
    l0 = np.sum(inertia*cross(x, v)); mass0 = mass.sum()
    maxp = 0.; flips = 0; max_remap_force = 0.; remap_u = 0.
    samples = [y.copy()]; initial_state = y.copy()
    nstep = round(duration/dt)
    assert abs(nstep*dt-duration) < 1e-12
    for _ in range(nstep):
        k1 = rhs(y); k2 = rhs(y+dt*k1/2); k3 = rhs(y+dt*k2/2); k4 = rhs(y+dt*k3)
        y = y+dt*(k1+2*k2+2*k3+k4)/6
        x = y[:8].reshape(4, 2)
        before, pressure, area = pressure_force(x, keys, mass, entropy, p0)
        maxp = max(maxp, float(np.max(abs(pressure/p0-1))))
        if remap:
            newkeys = mesh(x)
            if newkeys != keys:
                u = pressure*area/(GAMMA-1)
                mass, un, _ = overlap_remap(x[np.asarray(keys)], x[np.asarray(newkeys)], mass, u)
                remap_u += un.sum()-u.sum()
                keys = newkeys; area, _ = geom(x, keys)
                entropy = (GAMMA-1)*un/area/(mass/area)**GAMMA
                after, pressure, _ = pressure_force(x, keys, mass, entropy, p0)
                max_remap_force = max(max_remap_force, float(np.linalg.norm(after-before)/(p0*h)))
                maxp = max(maxp, float(np.max(abs(pressure/p0-1))))
                flips += 1
        samples.append(y.copy())
    x = y[:8].reshape(4, 2); v = y[8:16].reshape(4, 2)
    radius = np.linalg.norm(xinitial, axis=1)
    theta = np.arctan2(xinitial[:, 1], xinitial[:, 0])+duration/radius**1.5
    exact = np.column_stack((radius*np.cos(theta), radius*np.sin(theta)))
    ef = energy(y); lz = np.sum(inertia*cross(x, v))
    projected = np.zeros(4); np.add.at(projected, np.asarray(keys).ravel(), np.repeat(mass/3, 3))
    result = dict(h=h, p0=p0, dt=dt, feedback=feedback, remap=remap, seed=seed, duration=duration,
                  max_pressure_deviation=maxp, position_error_h=float(np.linalg.norm(x-exact)/h),
                  energy_residual=float((ef-e0-y[-1])/abs(e0)),
                  thermal_scaled_residual=float((ef-e0-y[-1])/thermal0),
                  control_work=float(y[-1]), lz_relative=float((lz-l0)/abs(l0)),
                  mass_relative=float(mass.sum()/mass0-1), flips=flips,
                  remap_internal_change=float(remap_u), remap_force_jump=max_remap_force,
                  inertia_l1=float(np.abs(projected-inertia).sum()/inertia.sum()))
    return result, np.asarray(samples), initial_state


def phase_norm(delta):
    return np.sqrt(np.sum(delta[..., :8]**2, axis=-1)+np.sum(delta[..., 8:16]**2, axis=-1)/OMEGA**2)


def suite():
    for p0 in (2.78e-5, .00278, .278):
        for dt in (.0025, .00125):
            for feedback in (0, 1):
                for remap in (False, True):
                    base, _, _ = run(.1, p0, dt, feedback, remap)
                    for eps in (1e-5, 1e-6):
                        plus, yp, _ = run(.1, p0, dt, feedback, remap, eps)
                        minus, ym, _ = run(.1, p0, dt, feedback, remap, -eps)
                        gain = phase_norm(yp-ym)/phase_norm(yp[0]-ym[0])
                        record = dict(base, seed_amplitude=eps, gain_final=float(gain[-1]),
                                      gain_max=float(gain.max()),
                                      pair_flips=[plus['flips'], minus['flips']],
                                      pair_energy_max=max(abs(plus['energy_residual']), abs(minus['energy_residual'])))
                        print(json.dumps(record), flush=True)


if __name__ == '__main__':
    suite()
