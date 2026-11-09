# Multi-cell compatibility of independently solved material faces

Diagnostic continuation of the user's heated-pair example. Reuses existing
exact Riemann C header and native1024-cell Gresho geometry. Actual cell faces,
neighbor inference checks and least-squares velocity reconstruction are as in
WEIGHT_VOLUME_PROJECTION_SCREEN_2026-11-06.md. Density, pressure and velocity
are linearly reconstructed at each face, with positive traces asserted; no
limiter or new physics added. This is a smooth-state algebraic screen.

For each face solve the independent 2D normal Riemann problem, yielding u*.
With prescribed generator velocities, geometric face speed is Nv+K Wdot.
Choose the least-squares weight rate in the face norm diag(A*d):

    L Wdot = D (u* - Nv), L=DK.

All quantities use the same native faces. This matches net target volume
rate to <=4.70e-18 in the three states, but cannot match all target face speeds.

| state | post-fit speed mismatch RMS | maximum |
|---|---:|---:|
| initial | .00355966 | .0162570 |
| deformed frozen .1 | .00443698 | .0298616 |
| deformed dynamic .1 | .00507543 | .0326122 |

RMS is weighted by A*d. At the exact contact speed, the existing ALE fan has
mass flux <=2.23e-16. Sampling the same fan at the fitted geometric boundary
instead gives max face mass flux .016209,.029523,.032540 respectively. Sum
absolute face mass flux times area is .138619,.167806,.169605. Net cell mass
flux is much smaller (~1e-5), illustrating why net-cell checks can conceal
opposing face transfers. These are hypothetical flux evaluations, NOT actual
mass changes in the fixed-mass native trajectory; no evolution was performed.

Reason: a scalar weight rate defines an edge difference field. Independent
pair requirements need not satisfy the cycle sums of a single scalar field.
A volume projection removes only the divergence mismatch; a circulation-like
face residual remains. The scalar least-squares solution is exact within its
stated linear fixed-geometry problem, so tightening that solver cannot make
all these residuals vanish. This is not a theorem excluding every coupled
multidimensional Lagrangian discretization or all dynamic-weight laws.

Design implication: do not independently compute contact speeds, project them,
and then claim every face remains an exact material contact with the same
Riemann solution. Do not add mass transfer to hide the mismatch under the
user's fixed-mass scope. Next formulation must couple the geometrically
representable boundary motion and pressure/work conditions, while quantifying
its approximation to the local Riemann solution and entropy generation.
A simple pressure averaging or an extra weight inertia does not provide this
closure. Production algorithms and previous successful1D/Kepler results are
unchanged. No3D test, full2D qualification or performance claim follows.

Artifact: Research/MaterialVolume_20261105/material_face_compatibility_20261106_v1.
riemann.c uses the existing exact solver; audit.py/results.json retain each
state and both contact-speed/geometric-speed flux batches. All tools exited0.
