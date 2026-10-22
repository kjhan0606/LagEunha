# Astra decision: common enriched Euler representation

2026-10-21. One bundled review of `DG_PROJECTION_RESULTS.md`, independent
central-cell integration, angular identities and legacy compatibility. Driver
implements/runs; Astra plans and evaluates bundles. This supersedes the pure
polynomial choice in `ASTRA_EULER_SPATIAL_PLAN.md`. No enriched evolution has
been implemented or certified by this decision.

## Selected space and physical scope

Use the same scalar space for all four Euler components:

    span{P2, rho_eq, (rho vx)_eq, (rho vy)_eq, E_eq}.

All coefficients, including the enrichment coefficients, evolve through the
full Euler weak form. Initial equilibrium data are disclosed basis functions,
not immutable state coefficients. No stationary-residual subtraction, state
reset, independent kinetic reservoir, thermal refill, viscosity or spin
damping. Original initial fields, epsilon and layout remain unchanged.
This known-field enrichment is not evidence of generic Kepler accuracy when
the initial field is unavailable. A common scalar space is needed for linear
Galilean transformations of the conservative variables.

Stop raising global pure polynomial degree: P2/P4/P6 improve most cells but
all fail the cold-pressure gate on this original initial field. The narrow
central density structure is one measured projection obstacle, not an
alternative explanation that erases the historical native work error.

## Conditioning and integration

Normalize physical units and orthogonalize reference functions against P2.
Near-polynomial functions require explicit conditioning diagnostics. No ridge
or diagonal regularization. Rank exclusion must quantify omitted conserved
field and cold-pressure representation error, not merely a small coefficient.

Use origin-aware segmented positive quadrature consistently in Gram, initial
moments, full volume flux/source and relevant faces. Selected-cell pressure
convergence near 1e-9 does not certify arbitrary moving-stage quadrature.
The central radial derivative is undefined at one measure-zero point; nearby
gradients are real. Do not smooth the initial field to conceal them.

## Evolving moments and moving cells

State is M_ia = integral_cell U phi_ia. Full Euler face flux and actual
gravitational source determine its RHS. A time-dependent basis requires
integral U partial_t(phi). Prefer physical-coordinate basis fixed within a
timestep, with a consistent conservative transformation/recentering at step
boundaries. Never convex-mix coefficients expressed in different bases.

Use one canonical shared-face quadrature and actual face speed. Couple the
weight-mass solve to the actual Euler mass flux. Check the periodic outer
boundary: a radial reference field is not automatically a global stationary
solution on a periodic square.

## Positivity and evidence requirements

Conservative scaling U_theta = U_mean + theta (U-U_mean), 0<=theta<=1, applies
to every nonconstant mode, including enrichment. Preserve mean mass,p,E.
An inadmissible mean rejects the stage. Record limiter cells, scale and
kinetic/thermal redistribution. Initial exact-field reconstruction must not
need this limiter. Positivity at quadrature points is not a whole-cell proof,
and a positivity limiter alone does not establish entropy stability.

## Continuous next bundle

Driver status2026-10-21: `ENRICHED_EULER_RESULTS.md` implements native CUDA
initialization and fixed-cell Euler helper checks. Full moving-cell RHS,
weight/mass coupling, runtime state and evolution gates below are still open.
The new initial ledger discloses continuous-field versus legacy point-state
momentum/E differences. This partial implementation does not change the
binding sequence or certify a long run.

1. Same71080-layout enriched initialization AND complete Euler RHS: conservative
   means, independent cold-pressure/core4 quadrature, inactive initial limiter,
   complete weak form/gravity/weight-mass relation, MPI1/4, checkpoint/restart,
   boosts/origin transformations. Explicitly expose conditioning failures.
   Projection alone is insufficient to start evolution.
2. A disclosed nonstationary perturbation: polynomial and enrichment
   coefficients must change through Euler RHS, without state reset; compare
   response and ledgers with independent full Euler flux/source.
3. Same Euler time/flux/positivity path on Sod/Sedov with Kepler enrichment OFF:
   shock speed, conserved ledgers, heating and limiter diagnostics.
4. Only after these pass, one fresh original refined Kepler continuation to
   t=.2, then t6-10 and the actual longest GIZMO orbital duration.

Stop this candidate if it needs energy repair/IC changes, uncontrolled
conditioning, broken moving-basis or weight-mass relations, equilibrium
forcing, incorrect shock heating, thermal depletion, timestep collapse,
secular weights or lost valid geometry. Preserve outputs and restarts.
