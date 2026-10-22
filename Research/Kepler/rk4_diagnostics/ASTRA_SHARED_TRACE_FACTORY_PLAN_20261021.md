# Astra: prospective native shared-trace factory bundle

2026-10-21. End-of-whole-target bundle evaluation. Driver codes/runs; no
additional microreview is needed within this coherent bundle.

## Decision

End promotion of the old broken packet. Target68526 source118/246/254
certificates remain incomplete and canonical potential seam bound alone
exceeds unchanged tolerance. Reject/no-publish/caller preservation is the
correct result. No further old-packet integration or guard search.
This does not prove the physical solution is wrong/nonexistent.

## Native state and basis contract

Build a SEPARATE fresh representation control from the SAME original
continuous IC and SAME24-cell reference mesh. Do not average/project/refit
old coefficients, reset an existing physical trajectory, or use a future
orbit/gravity function as enrichment.

Use shared global vertex DOFs, shared canonical-edge DOFs and local interior
bubbles with STRUCTURALLY ZERO trace. Every shared DOF has a single global
ID and owner. Update it once and scatter it; do not reconcile independently
updated cell coefficients afterwards.

Preserve original nonpolynomial initial velocity with the decomposition

    v0 = I3(v0) + sum_edges E_e(r_e) + b_T
    r_e = v0|e - I3(v0)|e
    b_T = v0 - I3(v0) - sum_edges E_e(r_e).

Canonical edge r_e is defined ONCE. A fixed reference extension equals r_e
on that edge and zero on other edges. Interior b_T is zero on all edges.
Treat residual trace dependencies in the CANONICAL EDGE space, not by
tying residual coefficients over an entire connected component. Record
new dimension and relation to the old conforming space. Acceleration need
not match the rejected representation's acceleration.

The same conforming space represents x and v:

    x_h=sum N_a x_a, v_h=sum N_a v_a,
    M_ab=integral N_a N_b dm,
    M a=F_sigma+F_g+F_ext,
    rho e_dot=-P div(v).

Use the SAME stress and quadrature to connect Kdot and physical Udot.
Measure K=.5 v^T M v, L=integral x cross v dm and true mass-centred internal
angular momentum. No independent spin/kinetic reservoir, damping or heat
correction. Kepler Re=infinity and no generator steering.

Assembly uses sparse restriction/scatter. No production-wide dense Schur.
If local traction ledgers are required, recover through the existing hybrid
relation, without residual-fitting forces or heat.

## One bounded implementation/verification bundle

1. Native factory on same24-cell mesh. Preserve original rho/v/P and mass,
   physical U/K/L integrals. Explicitly report definition differences from
   the original legacy point-state.
2. Original INITIAL RHS and finite h=2e-4 CONTROL generated with newly
   computed shared acceleration. This is not an accepted timestep.
3. CPU/OpenMP/MPI/CUDA evaluate SAME state. Check canonical edge x/v trace,
   orientation reversal, zero bubble trace, Piola stress/work, energy,
   angular momentum and resolved rotation, positive J and thermal energy.
4. Native checkpoint roundtrip stores IDs/owners/epoch/shared coefficients
   and local bubbles. Duplicate/missing/epoch mismatch leaves output intact.
   Restart cannot regenerate independent per-cell trace coefficients.

Keep old initial curvature criteria max chi<=.002 and every cell mean
chi<=.001, plus existing conservation thresholds. Continuity alone cannot
replace these physical accuracy gates. No degree/mesh/tolerance scan to
turn a failure into a pass; preserve any defect witness.

## Next connection

On pass, verify the SAME whole target68526 once with the new packet. Exact
shared trace removes the stored representation's seam difference, but does
not automatically solve the three inverse/analytic-domain certificates.
Then full targets, screened weights, feasible cold reexpression and
transactional commit; finally native CUDA/MPI real t6/8/10 and the longest
specified GIZMO-duration evolution. Full long objective remains ACTIVE.
