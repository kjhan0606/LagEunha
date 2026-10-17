"""Necessary local feasibility gate for fixed generator masses.

Corner masses must be nonnegative, sum to the fixed mass of each generator,
and to the remapped mass of each triangle. This is NOT a dynamics closure.
"""
import numpy as np
from scipy.optimize import linprog
from subcell_exact_geometry import Polygon


def audit(x, rho=(1., .01)):
    x = np.asarray(x, float)
    old = [(0, 1, 2), (0, 2, 3)]; new = [(0, 1, 3), (1, 2, 3)]
    oldpolys = [Polygon(x[list(k)]) for k in old]
    newpolys = [Polygon(x[list(k)]) for k in new]
    oldmass = np.array([p.area*r for p, r in zip(oldpolys, rho)])
    fixed = np.zeros(4)
    for ids, m in zip(old, oldmass): fixed[list(ids)] += m/3
    target = np.array([sum(p.intersection(q).area*r for q, r in zip(oldpolys, rho)) for p in newpolys])
    edges = [(i, t) for t, ids in enumerate(new) for i in ids]
    A = np.zeros((6, 6))
    for k, (i, t) in enumerate(edges): A[i, k] = 1; A[4+t, k] = 1
    rhs = np.r_[fixed, target]
    result = linprog(np.zeros(6), A_eq=A, b_eq=rhs, bounds=(0, None), method='highs')
    out = dict(feasible=bool(result.success), status=int(result.status), fixed_mass=fixed.tolist(),
               triangle_mass=target.tolist(), degree_one_deficit=float(max(fixed[0]-target[0], fixed[2]-target[1])))
    if result.success:
        out['maximum_constraint_residual'] = float(np.max(abs(A@result.x-rhs)))
        out['corner_mass'] = result.x.tolist()
    return out
