"""Candidate copower-circle times on polynomial site/weight paths.

Four copower-circular sites need not form a feasible mesh event. Callers must
verify adjacency and actual geometry; this module never changes a trajectory.
"""
import itertools
import numpy as np
from numpy.polynomial import Polynomial


def determinant(paths):
    paths = np.asarray(paths,dtype=float)
    if paths.ndim!=3 or paths.shape[:2]!=(4,3) or not np.all(np.isfinite(paths)):
        raise ValueError('four finite x,y,weight polynomial paths required')
    rows = []
    for site in paths[1:]:
        dx,dy,dw = (Polynomial(v) for v in site-paths[0])
        rows.append((dx,dy,dx*dx+dy*dy-dw))
    result = Polynomial([0.])
    for permutation in itertools.permutations(range(3)):
        inversions = sum(permutation[i]>permutation[j] for i in range(3) for j in range(i+1,3))
        term = Polynomial([(-1.)**inversions])
        for row,col in enumerate(permutation):
            term *= rows[row][col]
        result += term
    return result


def candidate_times(paths):
    polynomial = determinant(paths)
    scale = np.max(abs(polynomial.coef))
    if scale==0:
        raise ValueError('persistent copower degeneracy requires separate handling')
    roots = (polynomial/scale).roots()
    return sorted(float(z.real) for z in roots
                  if abs(z.imag)<1e-9*(1+abs(z.real)) and 0<z.real<=1)
