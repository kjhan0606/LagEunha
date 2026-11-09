# Admissibility proposals: exact-solver checks — 2026-11-09

Evidence: `common_admissibility_analysis_20261109_v1/{run.py,probe.c,report.json}`. These are local algebraic evaluations with the unchanged C exact Riemann response, not a new evolved hydrodynamic candidate. No new pressure floor, viscosity, sensor or limiter is adopted.

## A joint primitive blend is not a convex non-vacuum set

Let each primitive vary linearly from its owner state to its reconstructed state with a common factor theta. Positive density/pressure alone define an interval; existence of a shared material contact additionally requires

\[
 u_R-u_L < \frac{2(c_L+c_R)}{\gamma-1}.
\]

Take owner densities(1,100), both pressures1, candidate densities(20,30), both pressures5, and fixed normalized separation

\[
 \frac{(\gamma-1)(u_R-u_L)}{2\sqrt{\gamma}}=0.85.
\]

|theta|sqrt(P_L/rho_L)+sqrt(P_R/rho_R)|Existing C response|
|---|---:|---|
|0|1.1|material contact|
|0.5|0.7493569460|vacuum fan|
|1|0.9082482905|material contact|

The same response occurs for face normals along each dimension of the1D/2D/3D local API. Both endpoints are admissible while the midpoint is not. The density jump remains ordered and pressure jump is zero. Ordinary bisection cannot therefore claim the globally largest admissible factor on the basis of convexity or monotonicity. This disproves that particular shortcut, not every continuous reconstruction contraction. Physical owner-state vacuum cannot be repaired by contraction toward the owner.

## A velocity cap would conceal the pressure trace defect

The archived glass failure has positive reconstructed pressure9.094947e-19, despite a cold owner scale near1e-6. Its normal separation is4.367136e-6, while the exact non-vacuum bound is7.387148e-9: a factor591.18 smaller. Holding pressure/density fixed and capping separation preserves pair mean velocity and does not reverse expansion into compression. However it suppresses that expansion and changes pressure work; it is not a physical correction to the erroneous scalar trace.

An exact1D rotated/common-velocity-subtracted replay of this face is vacuum. Replacing its separation by the floating-point number immediately below the bound returns a material-contact status, but p-star is only1.60109e-98. Merely obtaining a successful status does not restore a useful thermodynamic state. This is a diagnostic calculation, not a selected numerical margin or a code modification.

Algebraic Galilean invariance and ideal affine-shear preservation do not remove this issue. Moreover the existing constrained pressure solver caches velocities and differentiates pressure response at fixed velocity. A sound-speed-dependent velocity map would require its pressure dependence in both the evaluated flux and Jacobian. Adding it only to the free RHS would break the common-operator contract.

## Decision

Do not implement either proposal merely to make the failed face pass. Preserve the pressure/volume/work exchange and investigate the scalar fit and tangential extrapolation that created the nearly zero face pressure. The earlier log-offset pressure distortion is also independently demonstrated in `COMMON_RECONSTRUCTION_ABLATION_2026-11-09.md`; positivity of that log candidate does not justify its degraded Sedov accuracy.

The new automated geometry-capacity correction is separate and has been adopted only after identical-state regressions. Reconstruction remains unresolved. Full2D validation, subsequent3D and a matched-accuracy <=2x GIZMO runtime certificate remain outstanding.
