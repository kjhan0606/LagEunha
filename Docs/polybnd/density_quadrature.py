"""Explicit positive triangle quadrature for exponential density fitting."""
from functools import lru_cache
import numpy as np


def validate_order(value):
    order = int(value)
    if str(order) != str(value) or not 6 <= order <= 64:
        raise ValueError('density quadrature order must be an integer from 6 to 64')
    return order


@lru_cache(maxsize=59)
def nodes(order):
    z, w = np.polynomial.legendre.leggauss(order)
    return (z+1)/2, w/2


def rule(poly, order):
    order = validate_order(order)
    z, w = nodes(order)
    center = poly.mean(0); points = []; weights = []
    for a, b in zip(poly, np.roll(poly, -1, axis=0)):
        u = a-center; v = b-center
        determinant = abs(u[0]*v[1]-u[1]*v[0])
        for r, wr in zip(z, w):
            points.extend(center+r*u+((1-r)*z)[:, None]*v)
            weights.extend(determinant*(1-r)*wr*w)
    return np.asarray(points), np.asarray(weights)
