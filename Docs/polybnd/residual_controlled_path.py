"""Recursive interval acceptance using the dense-path mass residual.

The checker is independent of the ODE's embedded error estimate. Rejected
intervals are bisected; the accepted left endpoint feeds the right interval.
This is a sampled residual controller, not a rigorous continuous bound.
"""
import numpy as np
from scipy.integrate import solve_ivp


def integrate(fun, y0, a, b, check, target, rtol=1e-12, max_depth=14,
              progress=None):
    attempts = []

    def visit(left, right, initial, depth):
        sol = solve_ivp(fun, (left, right), initial, method='DOP853',
                        rtol=rtol, atol=rtol*.01, max_step=right-left,
                        dense_output=True)
        if not sol.success:
            raise RuntimeError(sol.message)
        segments = sol.sol.interpolants
        defects = np.asarray([check(s) for s in segments])
        if not np.all(np.isfinite(defects)):
            raise RuntimeError('nonfinite dense residual')
        defect = float(np.max(defects))
        accepted = target is None or defect <= target
        row = dict(a=left, b=right, depth=depth, residual=defect,
                   accepted=accepted, segments=len(segments), nfev=sol.nfev)
        attempts.append(row)
        if progress is not None:
            progress(row)
        if accepted:
            return segments, sol.y[:, -1]
        if depth >= max_depth:
            raise RuntimeError('residual refinement depth exhausted')
        middle = (left+right)/2
        if middle in (left, right):
            raise RuntimeError('interval cannot be bisected')
        ls, mid = visit(left, middle, initial, depth+1)
        rs, end = visit(middle, right, mid, depth+1)
        return ls+rs, end

    segments, end = visit(a, b, np.array(y0, copy=True), 0)
    return segments, end, attempts
