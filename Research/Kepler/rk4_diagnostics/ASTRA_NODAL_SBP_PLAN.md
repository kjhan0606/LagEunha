# Astra revised decision: nodal SBP kinetic/thermal-compatible Euler

Implementation status: `NODAL_SBP_RESULTS.md` records the completed native
spatial screening. Identities hold, but actual Kepler divergence and the
advected density-wave Jacobian fail; this flux/operator is excluded from
runtime/long evolution pending the next bundled choice. The plan below is
retained as the historical decision and its explicit acceptance gates.

2026-10-21. One bundled review of `ENRICHED_COUPLED_RHS_RESULTS.md` and its
actual native CUDA/MPI/independent-reference evidence. Driver codes/runs.
The fixed-reference enriched Galerkin transport is excluded from evolution:
its projected Euler tangent breaks cold kinetic/thermal cancellation even
with exact initial reconstruction. Do not repeat weight/face debugging as a
substitute for this representation failure. Preserve the reusable native
screened solve, shared flux, angular ledger and old evolution path.

## Selected volume flux and state

Use actual positive-diagonal-mass **nodal conserved Euler state**, not enriched
coefficients reprojected through the failed dense Gram inverse each stage.
For L,R primitive states and a common direction n and symmetric mesh speed w,

    rho_bar=(rho_L+rho_R)/2, v_bar=(v_L+v_R)/2, P_bar=(P_L+P_R)/2,
    a=v_bar.n-w,
    f_rho=rho_bar*a,
    f_momentum=f_rho*v_bar+P_bar*n,
    f_E=.5*f_rho*(v_L.v_R)+P_bar/(gamma-1)*a
           +.5*(P_L*v_R+P_R*v_L).n.

This is a symmetric consistent smooth **volume** flux, not a shock entropy
stability claim. Momentum must use v_bar and kinetic flux must use v_L.v_R,
not a mean of squared velocities. With appropriate SBP/GCL and constant P,
the nodal pressure equation must reduce to Pdot=-gamma P div_h(v), separating
chain-rule error from velocity-divergence discretization error.

## Required spatial/metric structure

Construct on actual polygons:

    M=diag(positive nodal weights),
    Q_k+Q_k^T=R^T B_k R, Q_k*1=0.

Volume/face nodes, trace R, normals and measures belong to the same operator.
First construct positive diagonal-norm SBP; then flux differencing and same-node
gravity; then actual moving-node metrics/GCL; then native runtime/time step.
Nonpolynomial quadrature alone is not SBP evidence. Do not subtract the exact
Euler residual to make identities hold. Do not change IC or perform a new
degree sweep to hide failure.

Canonical full Euler HLLC/HLL fan may remain the boundary numerical flux,
connected by the proper SAT term to the split volume flux. Shock kinetic
dissipation stays inside total E and becomes physical Euler thermal energy.
No Pi/R reservoir, thermal repair or finite physical Re.

At each node use S_p=rho*g and S_E=v.S_p from identical states/precision;
their thermal contribution cancels. Check total potential-energy ledger
separately. Reuse screened mass/weight solve with the **new actual boundary
mass flux**. Internal flux differencing must not produce a cell mass source.
The fixed-fan Jacobian proof does not prove dynamic geometry existence.

## One next bundle, A through D

A. On the original failure cells and neighbors, prove SBP, row sums, metric
   closure, actual GCL, shared-face conservation/MPI, nodal kinetic/thermal
   identity and same-node gravity cancellation. Use the same exact initial
   state as the failed tangent test. Initial quadrature error must be disclosed.
B. Same b=(.013,-.007) perturbation: separately report numerical Pdot,
   -gamma P div_h(v), their residual, and div_h(v)-div(v). Small chain-rule
   error with large divergence error does not pass cold Kepler.
C. Actual discrete Jacobian for one nonuniform density-wave state at uniform
   P,v. If a continuum-absent growing mode occurs, preserve the failure and
   do not start long evolution or append arbitrary damping.
D. Only after A-C: native moving state one-step, MPI/restart, then same-path
   Sod/Sedov shock speeds, conservation, heating and positivity; then one
   refined Kepler path through t=.2, t6-10 and paper duration. Conservative
   remap and its errors are part of runtime if nodal relocation is needed;
   resampling at fresh quadrature nodes is not a substitute.

KEP/PEP or entropy properties do not establish local linear stability:
[Ranocha and Gassner](https://arxiv.org/abs/2009.13139). Moving-mesh entropy
arguments require their SBP/metric assumptions:
[Schnuecke et al.](https://arxiv.org/abs/1812.09093).

## Initial state and final scope

Record the fresh continuous-state versus legacy point-state energy difference
as different initial discretizations; add/subtract no heat. Future conservation
baseline is the explicitly adopted fresh initial integral. The periodic radial
exterior has real boundary angular transport; retain it and do not force it
stationary. Match GIZMO IC/boundaries for the final comparison independently.
Steering remains forbidden and Kepler Re=infinity. No nodal evolution is
implemented or certified by this plan.
