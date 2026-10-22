# Fixed P3 material-map bundle

2026-10-21. Records the existing Astra planner's bundled decision after native
hybrid v1-v5. The hybrid ledger is retained. Its initial thermal curvature
error must be addressed before curved Kepler remap or evolution promotion.

## Selected spaces only

Use complete C0 cubic material-reference polynomials plus two original input
velocity residuals v0_k-I3(v0_k). Residual coefficients are ordinary dynamic
kinematic unknowns. No gravity, gravity derivatives, or future exact orbits
are added to the basis. Initial x, rho, v, physical e and conserved integrals
remain the original continuous fields. No thermal correction is applied.

Evaluate exactly two representations:

1. Original96 triangles, P3.
2. Each triangle split into4, same P3 (h/2).

Original24 physical cells, generator positions, weights, initial cell masses
and shared physical boundaries remain unchanged. This is internal kinematic
resolution, not a new graded-generator IC. Preserve the old inner-refinement
evidence separately.

Keep M a=F+B^T lambda, B a=0, and physical Udot=-int P div(v). Internal mesh
interfaces must be continuous. Derive each physical-face trace Gram/rank from
actual P3 and velocity-residual functions; do not extend axis-specific
incidence exceptions. Choose rank/gauge without loads. Explicitly reject
uncertain rank. Physical-cell bubbles must not be tied unnecessarily.

## One bounded execution bundle

Keep gravity-off, linear gravity, actual softened Kepler, and Kepler+uniform
velocity-boost controls. Record at each representation: native MPI1/4,
momentum/work/angular and bubble/face ledgers, mass positivity/conditioning,
a_h-g and its full gradient, curvature C_h, point maxima and locations, cell
volume L1, quadrature sensitivity, DOFs/memory/assembly/solve costs. The dense
reference Schur is not the planned full71080-cell production solver.

Apply the already validated finite-map derivative check once to the final
space. Do not run a timestep search or increase to P4/P5 if this bundle fails.

## Fixed engineering promotion gate (not a long-time guarantee)

For actual Kepler define chi=abs(C_h)/Omega^2. Final h/2 requires:

- maximum chi <=2e-3;
- every cell's volume mean chi <=1e-3;
- quadrature/linear-solver uncertainties separately small against the gate;
- both curvature maximum and volume L1 improve at least4x under h/2, except
  quantities already at numerical uncertainty;
- existing face, conservative work/angular, positive-mass and MPI gates hold.

Two levels do not establish an asymptotic convergence order. The chi bound
controls the initial pressure second-Taylor contribution on Omega^-1 only;
it does not bound actual evolved pressure over that time or forecast t10.

If failed, preserve locations, derivative error and costs and close this
bundle; do not bypass by changing gravity, heat, timestep or higher degrees.
If passed, next is actual nonuniform Kepler curved material-map to Laguerre
remap and cold re-expression, then composed density-wave/sliding and shocks,
then native scalable CUDA evolution through t6-10 and longest-paper duration.
