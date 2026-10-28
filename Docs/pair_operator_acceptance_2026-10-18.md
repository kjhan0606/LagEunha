# Static pair operator: what a converged fit would establish

The current refinement continues the saved interior-constrained solution;
it does not restart a physical simulation. Column scaling is applied only to
the least-squares correction variables, with zero columns left unchanged.
The required maximum interior constraint residual remains 1e-8.

## Negative coefficients are not by themselves an acoustic instability proof

For a fixed geometry and compatible semidiscrete equations

    Vdot = B v,   M vdot = B^T P,   Udot_i = -P_i (Bv)_i,

linearize about a stationary state with constant pressure, positive volumes,
and positive masses. Ignore gravity, AV, moving geometry, and boundary exchange
for this limited argument. With D=diag(gamma P/V), pressure perturbations obey
delta P=-D delta V and

    delta Vddot = -B M^-1 B^T D delta V.

The matrix is similar to the positive-semidefinite matrix
D^(1/2) B M^-1 B^T D^(1/2). Thus negative scalar pair coefficients do not alone
imply growing frozen-geometry acoustic modes. Null modes may still exist, and
positivity/monotonicity, shocks, geometry dependence, rotation and boundary work
are separate requirements. This argument requires an equilibrium and cannot
excuse a residual force under uniform pressure.

## Remaining gates even if the algebraic fit converges

1. Constant pressure closure, affine moments, uniform translation and rigid
   rotation must pass in the intended region; outer reaction/boundary work must
   be explicit. An interior fit is not a full-domain solver.
2. Establish the physical volume associated with B. Substituting Bv for the
   derivative of actual Voronoi volume is not automatically consistent with
   rho=m/V. The original pressure-volume adjoint already passes its own tests.
3. Any Laguerre realization must satisfy the approved fixed-mass reduction for
   persistently equal weights and account for weight work/transport once.
4. Verify physical-time nonlinear evolution, positivity, shocks, and the
   historical t6–10 and longer intervals with the original implementation.

This static investigation addresses an unresolved feasibility question. It is
not authorization to replace geometric thermodynamics merely because the fit
reduces circular divergence, nor evidence of long-term Kepler stability.
