# Quadratic positive-density reconstruction: first-remap gate

Fit rho=(m/V)*b.f(z), z=(x-X_generator)/sqrt(V), with six basis functions
f=(1,zx,zy,zx^2,zx*zy,zy^2). Solve the normalized six-moment Gram system
using the actual transported m,q,Q from the preceding uniform donor remaps.
Degree-four products are integrated with positive order4 polygon quadrature.
No mass centroid or covariance is reset to uniform geometric values.

Check the continuous minimum of the quadratic over each convex polygon:
vertices, convex edge stationary points, and an interior positive-definite
stationary point. Interior minima with singular Hessian extend to a boundary,
so they are covered by the boundary candidates. This is not merely a sampled
quadrature-point positivity check. Numerical comparisons use stated tolerances.

| Weight fraction | Negative cells /1590 | Minimum rho / cell mean rho |
|---|---:|---:|
| .25 | 8 | -.745480 |
| .125 | 3 | -.288715 |
| .0625 | 0 | .139644 |
| .03125 | 0 | .484596 |

Maximum normalized moment-equation residual is 4.45e-16; maximum Gram condition
number is 361.41. The two larger-step failures are not roundoff or failed linear
solves. No coefficients are clipped and no mass/spatial moments are relaxed.

## Decision

The six-moment quadratic density is a usable candidate for the next diagnostic
at increments .0625 and .03125, with positivity checked anew after every
transition. It is not universally positive, and the successful first remaps do
not prove the full fixed path will pass. Use the actual q,Q and the previously
tested G=C^T S^-1; keep residual kinetic covariance and thermal energy separate.
Check m,p,D,Pi,U plus represented angular momentum through the full path.

This is a compact density reconstruction gate only. It does not define a
physical face closure, prove original point-state reduction, or validate a
production solver. Existing long Kepler trajectories are not changed.

Evidence: `canonical_quadratic_density_reconstruction_2026-10-18.json`.
Script: `audit_quadratic_density_reconstruction.py`.
Remote: `/gpfs/kjhan/LagForce/Codex_quadratic_density_reconstruction_20261018`.
