# Astra revised spatial decision: conservative moving-cell P2 Euler

2026-10-20. This records one bundled initialization review and a correction
after driver checked the proposed tensor energy identity. Driver codes and
runs the native implementation; Astra plans/evaluates coherent bundles.

## Binding decision and rejected alternatives

The initial affine state m,q,Q,p,D,E is useful algebra/diagnostics but is not
an adequate evolving thermal closure for the current cold curved Kepler IC.
Full E-K_affine incorrectly assigns up to10.02% extra initial thermal energy.
Constant-specific-u introduces P proportional to varying rho inside a cell.

Astra first proposed independent Pi=int rho vv^T with zero third central
velocity moment and no isotropization. That decision is **withdrawn**: driver
showed that E-half-trace(Pi) then has only thermal advection and reversible
pressure-divergence work, including under a common augmented HLL flux.
Euler shock heating/entropy weak-limit consistency is not established. Do not
implement that model or disguise the omission with Pi→U relaxation.

Do not delete initial kinetic energy through an undisclosed projection or add
energy to depleted checkpoints. Do not evolve Pi with a different numerical
flux from E as an ad hoc correction. Re=infinity permits shock entropy
production while excluding an invented physical viscosity.

## Selected implementation

Implement a conservative **moving-cell P2 DG Euler representation** in the
existing C/OpenMP+MPI+CUDA body. Density, both momentum components and total
energy all use the same spatial order and one canonical Euler Riemann flux.
Existing actual face-speed, MPI ownership and weight–mass relation are reused;
this is not Python evolution, a new standalone production solver or removal
of the old native code. Pi, R, spin/strain energy are reconstructed diagnostics,
not an independently advected collisionless reservoir.

The polynomial representation must include momentum AND energy spatial
information. Merely adding scalar spin or quadratic velocity while retaining
an incompatible energy representation is insufficient. Internal angular
momentum is contained in the spatial momentum moments.

P2 is the first candidate, not a guarantee of adequate cold-flow accuracy.
Finite degree, exact conservative cell averages and exact pointwise pressure
cannot generally all match a nonpolynomial continuous IC. Keep cell averages
conservative and report pressure approximation error directly; do not repair
it by deleting kinetic energy or adding thermal energy.

## First projection gate

Use the existing71080 initial layout, same continuous rho,v,P, epsilon and
parameters. Do not introduce graded layout, limiter or physical changes in
the same comparison. Independently control quadrature error.

1. Project Euler conserved fields to P2 and preserve cell m,p,E integrals.
2. Reconstruct rho and pressure at independent volume/face points; require
   positivity and measure location-dependent cold-pressure error.
3. Separate projection-order error from integration error. Compare against
   the quantified affine10.02% initial bias without claiming that any positive
   result is automatically accurate.
4. If this representation fails the cold-pressure gate, change spatial order
   or representation before evolution. P2 alone is not a physics pass.

## Full flux/source and runtime gates

- Same canonical quadrature points/fluxes for Euler and spatial moments.
- Moving-basis/geometric conservation and uniform/affine ALE checks.
- Mass, both momentum components, total energy and angular ledgers; MPI
  shared-face agreement and complete migration/halo/checkpoint state.
- Positive rho and P at every stage; conservative positivity processing
  preserves cell means and records its activation. Full rollback includes
  all polynomial/weight/geometry state.
- Establish shock treatment/entropy behavior, then existing Sod/Sedov shock
  speeds, conservation and thermal generation. Euler flux alone does not
  prove a high-order method entropy-stable.
- Couple the actual Euler mass flux to the approved weight–mass relation.
  Complete Jacobian/monotonicity tests including HLL fallback; initial W0
  identity and the1206 HLLC samples are not a dynamic-weight solver proof.
- Only after those gates, one existing-layout fresh Kepler path to the first
  continuation gate. Preserve restarts and extend passing physics through
  t6–10 and the actual longest GIZMO orbital duration. Short positivity or
  initialization results never satisfy the final goal.

Status update 2026-10-21: P2/P4/P6 Euler projection is implemented and measured
in `DG_PROJECTION_RESULTS.md`; every degree fails the original initial-field
cold-pressure gate. DG time transport remains unimplemented. The next spatial
decision is `ASTRA_ENRICHED_EULER_PLAN.md`. No long run is launched from either
the rejected affine thermal closure or the failed polynomial projection.
