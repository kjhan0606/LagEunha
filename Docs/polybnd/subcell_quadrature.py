"""Positive finite realization of 16 subcell moments, at a fixed geometry.

Internal quadrature points are not moving particles or Voronoi generators.
This is a representation/compression test, not a transport or dynamics closure.
"""
import numpy as np
from dataclasses import replace
from subcell_exact_geometry import Polygon


def subdivide(fragments):
    """Exact geometric refinement without changing donor fields (test fixture)."""
    result = []
    for f in fragments:
        p = f.polygon.vertices
        c = tuple(sum(v[k] for v in p)/len(p) for k in (0, 1))
        for a, b in zip(p, p[1:]+p[:1]):
            poly = Polygon([c, a, b])
            if poly.exact_area > 0: result.append(replace(f, polygon=poly))
    return result


def quadrature(fragments):
    points = []; weights = []
    bary = np.array([[2/3, 1/6, 1/6], [1/6, 2/3, 1/6], [1/6, 1/6, 2/3]])
    for f in fragments:
        p = f.polygon.vertices
        for i in range(1, len(p)-1):
            a, b, c = p[0], p[i], p[i+1]
            area = ((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))/2
            if area == 0: continue
            if area < 0: raise RuntimeError('quadrature expects convex polygon')
            x = bary@np.array([a, b, c], dtype=float)
            v = f.b+x@f.B.T
            for xi, vi in zip(x, v):
                points.append(np.r_[xi, vi, f.internal_density/f.rho])
                weights.append(float(area)*f.rho/3)
    return np.asarray(points), np.asarray(weights)


def features(points, potential=None):
    x, y, vx, vy, e = points.T
    result = np.array([np.ones(len(x)), x, y, x*x, x*y, y*y, vx, vy,
                       x*vx, x*vy, y*vx, y*vy, vx*vx, vx*vy, vy*vy, e])
    if potential is not None: result = np.vstack((result, potential(points[:, :2])))
    return result


def compress(points, weights, potential=None):
    """Eliminate affine dependence of positive weighted moment vectors.

    At most 16 (17 with a fixed potential) retained points; no thermal correction.
    Numerical SVD means moment preservation is tested to floating precision.
    """
    points = points.copy(); weights = weights.copy()
    if (weights <= 0).any() or (points[:, -1] <= 0).any():
        raise ValueError('expected positive quadrature')
    limit = 16 if potential is None else 17
    block = limit+1
    while len(weights) > limit:
        f = features(points[:block], potential); scales = np.max(abs(f), axis=1)
        scales[scales == 0] = 1.  # normalization of identically zero moment rows
        _, _, vh = np.linalg.svd(f/scales[:, None], full_matrices=True)
        z = vh[-1].astype(np.longdouble)
        w = weights[:block].astype(np.longdouble)
        positive = np.flatnonzero(z > 0)
        if not len(positive): z = -z; positive = np.flatnonzero(z > 0)
        pivot = positive[np.argmin(w[positive]/z[positive])]
        new = (w*z[pivot]-w[pivot]*z)/z[pivot]
        new[pivot] = 0.  # the chosen exact elimination, not a physical floor
        if (new < 0).any(): raise RuntimeError('negative compression weight')
        weights[:block] = new.astype(float)
        keep = weights > 0
        points = points[keep]; weights = weights[keep]
    return points, weights
