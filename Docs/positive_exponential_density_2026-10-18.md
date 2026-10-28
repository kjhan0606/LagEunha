# Positive exponential density realizes the five rejected spatial states

Test only the five saved cell states rejected by the quadratic density in the
repeated path: N4 IDs665,691,808 and N8 IDs691,808. Same transported m,q,Q.
Use local coordinates z=(x-X)/sqrt(V) and normalized density

    rho=(m/V) exp(theta.f(z)-log Z),
    f=(zx,zy,zx^2,zx*zy,zy^2), Z=(1/V) integral exp(theta.f) dV.

Solve five normalized moment equations with analytic covariance Jacobian and
stable log normalization. Positivity follows from the finite exponential form.
The minimum is evaluated by minimizing its quadratic exponent over the polygon.

| State | Minimum rho / mean rho | Order16 fit checked with order32: maximum normalized mass/moment error |
|---|---:|---:|
| N4 ID665 | .314635 | 4.45e-16 |
| N4 ID691 | .212869 | 5.83e-15 |
| N4 ID808 | .229532 | 1.20e-15 |
| N8 ID691 | .290271 | 5.97e-16 |
| N8 ID808 | .296075 | 8.89e-16 |

Order8 fits checked using order16 have errors 3.90e-9 to 2.79e-8. Therefore
small nonlinear residual on one quadrature alone is not sufficient evidence
of matching the continuous moments. Order32 fits checked with order48 have
maximum error4.92e-15. Independent checks retain the fitted normalization,
so they test mass as well as normalized moments, rather than hiding mass error
by renormalizing on the check grid.

All five solves complete. Reported per-cell time for all three fits plus
independent checks is .014--.159 seconds on this run, including first-use
overhead. This is not a production cost estimate, GPU benchmark or a reason
to apply high-order quadrature to every cell at every RK stage.

## Next gate

A positive fallback exists for these rejected states without changing mass
or spatial moments. Test a hybrid reconstruction on the same saved path:
use the quadratic only when continuously positive; otherwise use the positive
exponential with independent quadrature checks. Integrate exponential donors
accurately on overlap fragments, and verify global moments without silently
renormalizing transfers. Any solve/positivity/conservation failure must stop
the case, not alter energy or weights.

This does not yet complete a weight path, establish spatial flux accuracy,
resolve residual covariance transport as gas physics, or prove equal-weight
fixed-mass reduction. No physical trajectory or production code was changed.
Long Kepler stability and full CUDA/OpenMP/MPI consistency remain open.

Evidence: `canonical_positive_exponential_density_2026-10-18.json`.
Script: `audit_positive_exponential_density.py`.
Remote: `/gpfs/kjhan/LagForce/Codex_positive_exponential_density_20261018`.
