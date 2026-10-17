"""Static fixed-mass support and Laguerre capacity reference (four sites).

Weights solve prescribed cell mass in a GIVEN spatial density field.
This supplies geometry, not a compressible-fluid evolution law.
"""
import numpy as np
from shapely.geometry import Polygon
from subcell_exact_geometry import Polygon as ExactPolygon, corner_polygons


OLD = [(0, 1, 2), (0, 2, 3)]
NEW = [(0, 1, 3), (1, 2, 3)]


def extended_support(x, rho):
    """Retain material labels across a flip, including nonincident support."""
    corners, labels = corner_polygons(np.asarray(x), OLD)
    targets = [ExactPolygon(np.asarray(x)[list(k)]) for k in NEW]
    allocation = np.zeros((4, 2)); fixed = np.zeros(4)
    for corner, (t, i) in zip(corners, labels):
        fixed[i] += rho[t]*corner.area
        for j, target in enumerate(targets):
            allocation[i, j] += rho[t]*corner.intersection(target).area
    nonincident = sum(allocation[i, t] for t, ids in enumerate(NEW) for i in range(4) if i not in ids)
    return dict(fixed=fixed, allocation=allocation, nonincident=nonincident)


def clip_halfplane(p, normal, bound):
    if len(p) == 0: return p
    result = []
    for a, b in zip(np.roll(p, 1, axis=0), p):
        da = np.dot(normal, a)-bound; db = np.dot(normal, b)-bound
        if (da <= 0) != (db <= 0):
            result.append(a+(b-a)*da/(da-db))
        if db <= 0: result.append(b)
    return np.asarray(result).reshape(-1, 2)


def power_cells(sites, weights, domain):
    cells = []
    for i, xi in enumerate(sites):
        p = np.asarray(domain, float).copy()
        for j, xj in enumerate(sites):
            if i == j: continue
            p = clip_halfplane(p, 2*(xj-xi), np.dot(xj, xj)-np.dot(xi, xi)+weights[i]-weights[j])
        cells.append(Polygon(p) if len(p) >= 3 else Polygon())
    return cells


class CapacityPatch:
    def __init__(self, domain, rho=(1., 1.)):
        self.domain = np.asarray(domain, float)
        self.donors = [Polygon(self.domain[list(k)]) for k in OLD]
        self.rho = np.asarray(rho)
        self.total = sum(p.area*r for p, r in zip(self.donors, rho))
        self.fixed = np.zeros(4)
        for ids, p, r in zip(OLD, self.donors, rho): self.fixed[list(ids)] += p.area*r/3
        self.scale = np.max(np.ptp(self.domain, axis=0))**2

    def evaluate(self, sites, weights):
        cells = power_cells(sites, weights, self.domain)
        mass = np.array([sum(c.intersection(p).area*r for p, r in zip(self.donors, self.rho)) for c in cells])
        volume = np.array([c.area for c in cells])
        return mass, volume, cells

    def solve(self, sites, target=None, initial=None, rtol=1e-11):
        target = self.fixed if target is None else np.asarray(target)
        if (target <= 0).any() or not np.isclose(target.sum(), self.total, rtol=1e-12, atol=0):
            raise ValueError('positive targets with correct total required')
        w = np.zeros(4) if initial is None else np.asarray(initial, float).copy()
        w -= w[-1]  # additive weight gauge
        eps = self.scale*1e-6
        for iteration in range(60):
            mass, _, _ = self.evaluate(sites, w)
            error = np.max(abs(mass-target)/target)
            if error < rtol: return w, iteration, error
            if (mass <= 0).any(): raise RuntimeError('Newton initial cell has zero mass')
            jac = np.zeros((3, 3))
            for j in range(3):
                dw = np.zeros(4); dw[j] = eps
                jac[:, j] = (self.evaluate(sites, w+dw)[0][:3]-self.evaluate(sites, w-dw)[0][:3])/(2*eps)
            step = np.r_[np.linalg.solve(jac, target[:3]-mass[:3]), 0.]
            for k in range(40):
                trial = w+step*2.**(-k); mt = self.evaluate(sites, trial)[0]
                if (mt > 0).all() and np.max(abs(mt-target)/target) < error:
                    w = trial; break
            else: raise RuntimeError('capacity Newton line search failed')
        raise RuntimeError('capacity Newton did not converge')

    def volume_jacobian(self, sites, weights, eps=None):
        eps = self.scale*1e-6 if eps is None else eps
        jac = np.zeros((4, 4))
        for j in range(4):
            dw = np.zeros(4); dw[j] = eps
            jac[:, j] = (self.evaluate(sites, weights+dw)[1]-self.evaluate(sites, weights-dw)[1])/(2*eps)
        return jac
