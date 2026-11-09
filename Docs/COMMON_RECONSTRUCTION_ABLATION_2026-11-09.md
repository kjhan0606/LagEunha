# Frozen reconstruction ablation — 2026-11-09

The previous goal turn made progress by completing the native candidate comparison, rejecting its accuracy regression, and preserving verified figures. This turn separates the scalar and velocity contributions before adding any new physical term. Production source and successful Kepler results are unchanged.

## Six operators on the same saved state

`common_reconstruction_factorial_20261109_v1` builds three scalar paths (original primitive hard-neighbor projector, continuous primitive face interval, current log-coordinate face interval) crossed with two velocity paths (original ordered-alpha, current tangent common-speed map). Primitive and log raw fits use the same-objective QR; the original hard projector remains an unchanged control. Original pressure-offset coefficients rho and log coefficients rho/P stay paired with their corresponding variables. No new evolution or tuned coefficient is introduced.

The original exact Riemann, wave impedance, direct force/volume/work and pressure-offset headers are identical. Five fixtures using the legacy branch reproduce the original portable-header RHS bitwise (`portable_parity/report.json`). The log/tangent branch reproduces the preceding Q audit wrapper bitwise. The tests are periodic 2D CPU saved-state evaluations, not a new CUDA/MPI/3D qualification.

Initial glass and the archived early failed glass state were reused directly. For completed glass control/candidate and regular control t=.1, independently saved polygons were converted to directed faces. Global-double clipping roundoff initially failed the native reciprocal-face check; these failed attempts are preserved. Recomputing the SAME identified power-plane intersections in long double changes geometry by at most2.72e-15 and gives native-volume discrepancies6.54–8.60e-16. No tolerance or physical state was changed. The candidate regular endpoint has three ambiguous tiny-edge owner assignments in this independent extraction; it is excluded rather than merged or declared passed. Full native shared faces are still needed for that fixture.

## Measured separation

At a fixed state the old ordered-alpha and new tangent velocity maps retain the same clipped normal jump. They change the pair's common velocity. Exact Riemann p-star and acoustic impedances therefore remain the same, while u-star changes. This explains identical direct force and different volume/work:

\[
\mathbf F=-A p_*\mathbf n,\qquad \dot V=A u_*,\qquad \dot E=-A p_*u_*.
\]

The glass endpoint pressure-star difference is below1.8e-15, normal-jump difference below1.6e-15, while the common-speed shift can reach0.40873 on the control state. This is an algebraic/frozen-state result; different volume histories can subsequently change forces.

| Same saved state | Scalar change: Fx/Fy relative L1 | Velocity change: Vdot relative L1 |
|---|---:|---:|
|glass initial|25.33% /31.37%|0|
|glass control t=.1|24.53% /23.69%|3.13%|
|glass candidate t=.1|23.29% /23.71%|4.19%|
|regular control t=.1|23.38% /23.38%|3.37%|

The initial acoustic heat ledger increases158.1553→193.4023; at the control glass endpoint1.42509→1.65106. This ledger is not identical to cell thermodynamic entropy. Primitive face intervals still produce reconstructed vacuum faces: the failed early state and both glass endpoints remain failures, regardless of smaller valid-face-only heat totals. There is no basis for simply reverting to that rejected path.

## Pressure-offset mechanism identified

`pressure_stages/report.json` follows actual initial face1055–1946. Both cell densities equal1 and both velocities are0, so velocity/density changes are excluded. Owner pressures are1e-6 and75.8833706. The original hot trace is39.0679644. Current log reconstruction produces75.8833510 before the interval map, and that map leaves it unchanged. Thus the near-owner hot trace is created by the fit, not clipping in this example.

For the hot owner the rho/P ratio makes the effective pressure displacement approximately(316002,440930), while the actual centroid separation is of order.01. The fitted log gradient is approximately(-8.0e-6,-3.54e-5). The hot-cell pressure is therefore nearly constant across its faces. Its partner trace and exact Riemann solution change the face p-star17.40943→33.81510, and the area-weighted acoustic ledger1.99544→5.40165. This is a directly reproduced mechanism, not proof that the original trace is the exact shock solution.

For equal density and two midpoint samples, the current log-offset ansatz gives the common reconstructed value

\[
P_f=P_L\left(P_R/P_L\right)^{P_R/(P_L+P_R)}.
\]

A finite linear-pressure interpolation instead gives(P_L+P_R)/2. At pressure ratio100 these are95.54 and50.5; at a weak ratio1.01 the relative difference is only1.24e-5. The differential identity grad(logP)=grad(P)/P does not make these finite reconstructions equivalent. This identifies an unnecessary strong-contrast distortion in the current hypothesis and the same large coefficients behind the earlier ill-conditioned pressure fit. It does not establish an entropy theorem or isolate every source of evolved Sedov error.

## Decision and next work

Keep the log/tangent combination unadopted. Retain QR and the inverse/audit convergence-contract fix as separately tested numerical corrections. The next scalar design must preserve the physical pressure/acceleration reconstruction while handling face admissibility and vanishing-neighbor continuity together. Do not add shock-specific floors, sensors, viscosity coefficients or silently restore the previously rejected primitive interval. Use the same saved failures, 1D reductions and full 2D tests to evaluate any derived replacement. No matched-accuracy GIZMO runtime or 3D claim follows from this audit.


## Subsequent admissibility checks

`COMMON_ADMISSIBILITY_ANALYSIS_2026-11-09.md` records exact-solver evidence against assuming a convex primitive non-vacuum blend and against hiding the pathological scalar trace with a velocity-only cap. Neither proposal is adopted.
