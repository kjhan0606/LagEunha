"""Instantaneous ALE balance on fixed-mass power cells in an exact steady flow.

rho=1, v=J*x, Phi=|x|^2/2, constant P: circular flow in harmonic gravity.
The finite polygon domain is open: boundary fluxes are included explicitly.
No time-stepping simulation or generator velocity correction is performed.
"""
import numpy as np
from shapely.geometry.polygon import orient
from weighted_mass_patch import CapacityPatch
from subcell_exact_geometry import Polygon
from subcell_moments import polygon_moments


def rotation(x):
    return np.stack((-x[..., 1], x[..., 0]), axis=-1)


def integrated(cells, pressure):
    state = []
    for c in cells:
        a, q, Q = polygon_moments(Polygon(list(c.exterior.coords)[:-1]))
        state.append(np.r_[a, rotation(q), np.trace(Q)+pressure*a/(5/3-1)])
    return np.array(state)


def audit():
    domain = np.array([[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]])
    sites = .8*domain+np.array([[.03, .01], [-.06, .04], [.02, -.02], [.01, .05]])
    velocity = rotation(sites); pressure = .1
    patch = CapacityPatch(domain); w, _, _ = patch.solve(sites)
    mass, area, cells = patch.evaluate(sites, w)
    eps = 2e-5
    # Derivative of weights follows the constant-mass constraint; generators
    # already move at the exact fluid velocity and receive no correction.
    wp, _, _ = patch.solve(sites+eps*velocity)
    wm, _, _ = patch.solve(sites-eps*velocity)
    wdot = (wp-wm)/(2*eps)
    cp = patch.evaluate(sites+eps*velocity, wp)[2]
    cm = patch.evaluate(sites-eps*velocity, wm)[2]
    derivative = (integrated(cp, pressure)-integrated(cm, pressure))/(2*eps)
    nodes, quadrature_weights = np.polynomial.legendre.leggauss(2)
    rhs = np.zeros((4, 4)); advective = np.zeros((4, 2)); absolute_mass_exchange = np.zeros(4)
    for i, c in enumerate(cells):
        p = np.array(orient(c, sign=1).exterior.coords[:-1])
        for a, b in zip(p, np.roll(p, -1, axis=0)):
            tangent = b-a; length = np.linalg.norm(tangent)
            if length == 0: continue
            normal = np.array([tangent[1], -tangent[0]])/length; mid = (a+b)/2
            neighbor = None
            for j in range(4):
                if j == i: continue
                aa = 2*(sites[j]-sites[i])
                bb = np.dot(sites[j], sites[j])-np.dot(sites[i], sites[i])+w[i]-w[j]
                if abs(np.dot(aa, mid)-bb) < 1e-10*patch.scale:
                    neighbor = j; break
            for q, weight in zip(nodes, quadrature_weights):
                x = mid+q*tangent/2; v = rotation(x); wn = 0.
                if neighbor is not None:
                    j = neighbor
                    ad = 2*(velocity[j]-velocity[i])
                    bd = 2*np.dot(sites[j], velocity[j])-2*np.dot(sites[i], velocity[i])+wdot[i]-wdot[j]
                    wn = (bd-np.dot(ad, x))/(2*np.linalg.norm(sites[j]-sites[i]))
                factor = weight*length/2; relative = np.dot(v, normal)-wn
                rhs[i, 0] -= relative*factor
                transport = v*relative*factor
                advective[i] += transport
                rhs[i, 1:3] -= transport+pressure*normal*factor
                energy = np.dot(x, x)+pressure/(5/3-1)  # KE + Phi + internal
                rhs[i, 3] -= (energy*relative+pressure*np.dot(v, normal))*factor
                absolute_mass_exchange[i] += abs(relative)*factor
        _, first, _ = polygon_moments(Polygon(p))
        rhs[i, 1:3] -= first  # integrated harmonic gravity
    return dict(max_balance_error=float(np.max(abs(derivative-rhs))),
                mass_flux_max=float(np.max(abs(rhs[:, 0]))),
                advective_momentum_max=float(np.max(np.linalg.norm(advective, axis=1))),
                mass_exchange_sum=float(absolute_mass_exchange.sum()),
                weight_rate_norm=float(np.linalg.norm(wdot)),
                omit_advective_error=float(np.max(abs(derivative[:, 1:3]-(rhs[:, 1:3]+advective)))),
                derivative=derivative.tolist(), rhs=rhs.tolist(), weight_rate=wdot.tolist())
