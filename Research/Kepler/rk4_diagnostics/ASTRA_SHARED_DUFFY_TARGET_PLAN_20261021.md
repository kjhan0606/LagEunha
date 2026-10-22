# Astra end review: shared canonical function -> complete target68526

2026-10-21. ONE end-of-function-connection bundle evaluation. Driver codes
and runs. Astra found no algebraic defect in material_shared_function.h,
material_shared_native_ball.h or material_shared_real_bound.h that reverses
the current results. Same-function Taylor residual/Bernstein bounds with
Lmax.02108 are valid progress, not state/tolerance modification. Reusing the
full817152-point proof and supplementing only missing inverses is appropriate.

Limits: real inside inverse success is not a complex inverse/remainder
certificate; native bounds are conditional on the stated8ulp libm model;
frozen heat is the generating definition for THIS saved CONTROL, not a
general evolved-thermal restart. evaluate() must reject the opposite vertex
in its Cartesian chart. Do not replace a zero chart denominator with epsilon.

## Next single coherent bundle

Same384 sources, target68526, saved finite state, degree, precision and
whole physical tolerance. Implement one transactional ALL12 target transfer.
No old broken packet integration or old failed proof promotion.

### Vertex Duffy coordinates and AD2

At original vertex c, opposite a/b:

    lambda_a=s(1-t), lambda_b=st, lambda_c=1-s,
    E_ab=s^4 r_ab(t).

Implement value/first/second derivatives directly, cancelling Cartesian
1/s algebraically BEFORE interval evaluation. A single whole-triangle Duffy
chart would become singular at the other vertices: instead use FIXED SIX
subtriangles from each original vertex, adjacent edge midpoints and source
centroid. Each subtriangle's original vertex is its Duffy apex. Actually
verify positive lower bounds for the other extension denominators.
This changes coordinates only; SAME function/state, no basis swap/refit.

### Integration routes follow actual geometry

- Parts without the vertex singularity: physical fan + certified complex
  inverse as before, connected to the NEW canonical function.
- Parts containing a vertex: REFERENCE Duffy pullback. Do not apply an
  ordinary2D inverse function theorem to the degenerate s=0 Duffy map.

Same target-face functions in both routes:

    f_j(s,t)=n_j dot chi(X(s,t))-b_j <=0.

Canonical edges own position, tangent, AD2 and target intersection events
once; both sides reuse these records. A target plane through an actual
source vertex and a tangency require explicit event handling. No tiny-face
removal or inferred ownership. q.boundary may enforce exact zero bubble
value ONLY with provenance that the query is the canonical edge, never a
near-edge flag used to discard uncertainty.

Reference density has r=sqrt(X^2+Y^2), a separate origin branch. At a chart
attached to origin explicitly use r=s*R(t) with continuous chosen R; do not
assume complex sqrt(s*s) produces that branch automatically.

### Remainder and physical-function ledger

Compute Gauss/Cauchy remainders only in coordinates where analyticity is
actually established. C3 is not a high-order analytic remainder theorem.
Include Jacobians, original density, frozen thermal generator, gravity
potential branches and all denominators.

Positive packet must represent the SAME physical integral:

    dm=rho0(X)dX,
    dp=dm*v, dL=dm*(x cross v), dK=dm*|v|^2/2, dU=dm*e.

Native cache error is evaluation error relative to canonical integration.
Exact shared analytic traces have no old broken-state finite seam gap;
do not reimpose that old physical gap. Do retain conversion/evaluation error
and complete coverage/overlap checks. Keep rho0/J push-forward distinct
from rho0/(J+delta) EOS and preserve the immutable heat generator.

### Required completion for this call

1. All384 source classifications and actual domain coverage complete, no
   missing/duplicated regions.
2. Required ACTUAL integration domains, including118/246/254, have certified
   inverse or Duffy routes and remainders.
3. Sum ALL12 source/field/geometry/quadrature/summation errors against OLD
   whole target criteria; partial allocations cannot replace total budget.
4. Preserve original heat/mass/state SHA, positive packets/MPI ownership and
   caller on failure.
5. Publish target68526 ONLY if all the above pass.

Within the fixed work budget, any unproved term/region is preserved and
this call rejects. No degree/precision/cap/guard scan and no full old-origin
proof regeneration. Pass ->full targets conservative assembly and actual
weight law. One target pass is not compression/evolution/t6-10 stability.

## Driver implementation notes (inferences, no new state changes)

Exact source-vertex/target-plane coincidences may require factoring a
provably zero apex value before isolating roots, not dividing by small s.
The fixed six chart domains must still be checked in their COMPLEX Cauchy
neighbourhood: real denominator positivity alone does not prove that larger
neighbourhood is free of singularities. Current code has no AD2/Duffy or
whole-target implementation yet; begin with those API connections.
