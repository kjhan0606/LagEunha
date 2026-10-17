"""Research reference: exact fragment transfer of piecewise-affine subcell state.

Authoritative state is the fragment list, not a uniform receiver average.
No nodal scatter, time evolution, or gravitational closure is claimed here.
"""
from dataclasses import dataclass, replace
import numpy as np
from subcell_exact_geometry import Polygon, corner_polygons

GAMMA = 5/3


@dataclass
class Fragment:
    polygon: object
    rho: float
    internal_density: float
    b: np.ndarray                 # v(x)=b+B*x
    B: np.ndarray


def polygon_moments(poly, origin=None):
    """Integrals of 1, x-origin, and (x-origin)(x-origin)^T; degree 2 exact."""
    if poly.geom_type != 'Polygon' or not poly.is_valid or poly.interiors or poly.area <= 0:
        raise ValueError('expected positive simple polygon')
    origin = np.zeros(2) if origin is None else np.asarray(origin)
    p = np.array(poly.vertices, dtype=float)-origin
    area = 0.; first = np.zeros(2); second = np.zeros((2, 2))
    for i in range(1, len(p)-1):
        v = p[[0, i, i+1]]; a, b, c = v
        d = b-a; e = c-a; at = (d[0]*e[1]-d[1]*e[0])/2
        total = v.sum(0)
        area += at; first += at*total/3
        second += at*(v.T@v+np.outer(total, total))/12
    if area <= 0: raise ValueError('nonpositive integrated area')
    return area, first, second


def validate_partition(polygons, rtol):
    for i, p in enumerate(polygons):
        if not p.is_valid or p.area <= 0: raise RuntimeError('invalid partition polygon')
        for j in range(i):
            if p.intersection(polygons[j]).exact_area > 0:
                raise RuntimeError('overlapping partition')


def transfer(fragments, receivers, rtol=1e-9):
    """Intersect, retaining donor rho,U,b,B on every fragment: no averaging."""
    donors = [f.polygon for f in fragments]
    validate_partition(donors, rtol); validate_partition(receivers, rtol)
    covered = [0 for _ in donors]; result = []
    receiver_error = 0.
    for receiver in receivers:
        pieces = []; area = 0
        for i in range(len(donors)):
            overlap = receiver.intersection(donors[i])
            if overlap.area == 0: continue
            if overlap.geom_type != 'Polygon': raise RuntimeError('unexpected nonconvex intersection')
            pieces.append(replace(fragments[i], polygon=overlap))
            covered[i] += overlap.exact_area; area += overlap.exact_area
        if area != receiver.exact_area: raise RuntimeError('incomplete subcell coverage')
        receiver_error = max(receiver_error, float(abs(area/receiver.exact_area-1)))
        result.append(pieces)
    if any(covered[i] != p.exact_area for i, p in enumerate(donors)):
        raise RuntimeError('incomplete subcell coverage')
    donor_error = max(float(abs(covered[i]/p.exact_area-1)) for i, p in enumerate(donors))
    if max(donor_error, receiver_error) > rtol:
        raise RuntimeError('incomplete subcell coverage')
    return result, max(donor_error, receiver_error)


def moments(fragments):
    """Mass-weighted affine projection + residual velocity covariance.

    C=integral rho (x-c)(v-vbar)^T. B=C^T S^-1 preserves P and L.
    Residual tensor R is retained as kinetic information, never added to U.
    """
    mass = 0.; first = np.zeros(2); momentum = np.zeros(2); internal = 0.; area = 0.
    for f in fragments:
        if f.rho <= 0 or f.internal_density <= 0: raise ValueError('nonpositive donor state')
        a, q, _ = polygon_moments(f.polygon)
        m = f.rho*a; mass += m; first += f.rho*q
        momentum += m*f.b+f.rho*f.B@q
        internal += f.internal_density*a; area += a
    center = first/mass; mean = momentum/mass
    S = np.zeros((2, 2)); C = np.zeros((2, 2)); T = np.zeros((2, 2))
    for f in fragments:
        a, q, Q = polygon_moments(f.polygon, center)
        w = f.b+f.B@center-mean
        S += f.rho*Q
        C += f.rho*(np.outer(q, w)+Q@f.B.T)
        T += f.rho*(a*np.outer(w, w)+np.outer(f.B@q, w)+np.outer(w, f.B@q)+f.B@Q@f.B.T)
    B = np.linalg.solve(S, C).T
    residual = T-B@S@B.T
    residual = (residual+residual.T)/2
    eig = np.linalg.eigvalsh(residual)
    scale = max(np.trace(T), mass*np.dot(mean, mean), np.finfo(float).tiny)
    if eig.min() < -512*np.finfo(float).eps*scale:
        raise RuntimeError('non-realizable residual covariance')
    spin = C[0, 1]-C[1, 0]
    angular = center[0]*momentum[1]-center[1]*momentum[0]+spin
    kinetic = .5*(mass*np.dot(mean, mean)+np.trace(T))
    return dict(mass=mass, first=first, momentum=momentum, internal=internal, area=area,
                center=center, mean=mean, S=S, C=C, B=B, residual=residual,
                residual_energy=.5*np.trace(residual), kinetic=kinetic,
                angular=angular, spin=spin, pressure=(GAMMA-1)*internal/area,
                residual_min_eigenvalue=eig.min())


def invariants(fragments):
    m = moments(fragments)
    return np.r_[m['mass'], m['first'], m['momentum'], m['angular'], m['kinetic'], m['internal']]
