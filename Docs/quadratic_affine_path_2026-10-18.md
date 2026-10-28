# Repeated quadratic-density affine remap stops on density positivity

Connected the preceding density fit and true-spatial-moment affine recovery in
a diagnostic on the same weight path 0 -> .25 -> 0. Use 4 or 8 subdivisions
per leg. All 17 moment slots are carried (potential slot is zero/unused).
Mean velocity is affine; residual kinetic covariance is transported separately
as a constant specific tensor per donor. Thermal energy density is constant
per donor, independently remapped. Positive order4 quadrature integrates the
polynomial moment products. No physical time advancement or production edit.

| Subdivisions per leg | Completed transitions | Weight at rejection | Failing density cells | Minimum rho / mean |
|---|---:|---:|---:|---:|
| 4 | 2 | .125 | 3 | -.289126 |
| 8 | 3 | .09375 | 2 | -.106401 |

Neither case reaches the outward endpoint or completes the round trip.
Reconstruction fails before using a negative-density state for the next remap.
No flooring, clipping, extra energy, or fallback projection is applied.

The failed cells are 665,691,808 for N4 and 691,808 for N8. They all have
positive thermal energies and positive residual eigenvalues. Across the patch,
minimum R eigenvalues are -5.94e-18 and -6.35e-18, within the specified floating
tolerance. Maximum aggregate moment changes before rejection are 1.37e-14 and
2.16e-14. Thus the encountered rejection is the density polynomial, not the
earlier inconsistent spatial-moment reset or an aggregate conservation failure.

## Decision

Small first steps are insufficient to guarantee positivity along the path.
Do not promote the quadratic density candidate to production or keep reducing
step size without addressing its finite-shape limitation. A quadratic polynomial
can represent the correct moments yet become negative inside a valid cell.
The moments have a positive incoming overlap representation; the polynomial
restriction is what fails in these cases.

Next evaluate a strictly positive density representation constrained by the
same spatial moments (for example an exponential quadratic with normalization),
first on the saved rejected states. Check independent quadrature convergence
and moment residuals; successful nonlinear fitting alone is insufficient.
Retain the overlap representation as reference and preserve failed outputs.

This does not imply that exponential reconstruction is affordable or that a
moment Euler closure is ready. Coupled face fluxes, original fixed-mass
equal-weight reduction, physical pressure/gravity, backend/MPI consistency,
and long Kepler stability remain open. Existing original runs are unchanged.

Evidence: `canonical_quadratic_affine_path_2026-10-18.json`.
Script: `audit_quadratic_affine_path.py`.
Remote: `/gpfs/kjhan/LagForce/Codex_quadratic_affine_path_20261018`.
