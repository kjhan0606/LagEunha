# Astra: actual curved remap bundle

2026-10-21. One bundled review after fixed P3/h2 and actual finite-state
CPU/OpenMP/MPI/CUDA controls. Coding remains with the driver.

The initial spatial gate and finite stress/work/angular implementation pass.
Sampled J positivity and degree6 ray algebra do not prove injectivity or
complete overlap detection. Use the same preserved h=2e-4 control state.

## Selected sequence

1. Bound J from below over every reference triangle; check continuous shared
   boundaries, simple exterior boundaries and separation of components.
2. Degree6 root isolation in angular charts, including endpoints, close roots
   and tangencies. Include current-map coefficient/evaluation uncertainty.
3. Interval enclosures on radial/angular panels classify source and target
   inequalities. Subdivide uncertain panels and bound possible m,p,L,K,U,V.
   Adaptive radial quadrature agreement alone cannot exclude a small lens.
   Fail explicitly at the fixed subdivision cap; never drop a small overlap.
4. Integrate actual transported fields on a single common refinement. Keep
   fragments as the authoritative source field. No initial-state reset.
5. One finite solve M_i(W)+mu_i W_i=m_i0. Jacobian is diag(mu)+L with canonical
   shared-face c_ij=int rho_current/(2|qj-qi|) ds >=0. The ray expression is
   equivalent at simple roots; a tangency is not handled by a denominator
   floor. Carry residual and Jacobian geometric uncertainty into convergence.
6. One conditional cold velocity compression and one reverse remap. Preserve
   failures and fragments. No unbounded degree/fragment/refinement sweep.

## Compression contract

First solve mass-weighted L2 projection with Bv=0, cell P_i v=p_i*,
A_i v=L_i*, and v^T M_i v/2<=K_i*. This is a coupled convex feasibility
problem, not independent translation/rotation fits. Restore any kinetic
deficit only with an actual interior bubble satisfying zero trace, p,L and
M-orthogonality to v. Its direction must be the transported-field residual
projected into that nullspace. If unavailable, reject compression; do not
create noise, heat, a kinetic reservoir or an invented compensating flux.

Check cell m,p,L,K,U,E and transported velocity L2, gradients/divergence,
traces, physical pressure and thermal-rate changes. p,L,K conservation alone
is insufficient. Keep density/thermal fragments if no positive conservative
compression exists; never assign thermal energy from E-K after projection.

## Mandatory limits

Normalize absolute K integration uncertainty by cold U as well as K. Record
fragment/point counts, memory and time. Reverse remap is a transfer check,
not time reversibility or long-time stability.

The dynamical W0 individually-fixed-mass reduction remains unproved. For
general finite geometry M_i(0)!=m_i0. Record this actual difference; do not
delete mass fluxes or overwrite overlap masses. h0 identity and the screened
weighted law do not resolve that compatibility issue.

After this bundle, implement sparse/matrix-free hybrid and composed remap
using the fragment/compressed state, then density-wave/sliding/shocks and
native CUDA/MPI evolution through t6/8/10 and the longest GIZMO duration.
Keep the existing graded generator layout; no new IC or time-integrator
sweep in this bundle.
