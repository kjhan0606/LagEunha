# Initial Gresho spatial-force audit

Native C/CUDA K1, Syntax A10, same executable as the saved evolution.
16^2,32^2,64^2 dynamic and 32^2 frozen, all rc0. No evolution rerun.
Uniform rho=1, W=pi=0, quadrature-averaged initial pressure; velocity sampled
at the generators. Artifact: Research/MaterialVolume_20261105/gresho_initial_force_20261106_v1.

| N | radial acceleration L1 | L2 | mean residual .18<=r<.22 |
|---|---:|---:|---:|
|16|.08577057|.17172070|+.67712451|
|32|.03004544|.08997747|+.37675084|
|64|.00805862|.03461950|+.10376663|

32^2 frozen and dynamic acceleration diagnostics are identical. The rigid
core r<.15 has RMS residual 1.47e-14 (32) and 2.96e-14 (64). At16 the broad
stencil crosses the transition and core RMS is .00539. Native forces agree
with the centered difference of cell-average pressure to <2.8e-7 (mixed
precision paths). Thus a substantial transition-region error exists before
time integration or weight evolution; refinement reduces it. The data do
not prove uniform second-order convergence at derivative breaks.

The auxiliary weight drive also exists at initial equilibrium:
Q_i=sum_j A_ij/(2d_ij)*(P_i-P_j). In the rigid core, P=5+12.5r^2, square
cells give Q=-25h^2. Native results at all three resolutions reproduce this
to <3e-13 in Q/h^2; the constant-pressure exterior has Q=0. Uniform Q within
the core alone changes a local common weight offset and is not itself a
local volume change. Its difference from the exterior is not a global gauge
shift. This identifies an additional weight response to the stationary
pressure profile, but does not prove it alone causes the measured outward
ring drift; translation, geometry changes and weight response interact.

Do not attribute the whole failure to Wdot or time integration, and do not
add an empirical limiter/radial correction. Frozen W0 remains the successful
Kepler control. Dynamic weight inertia remains experimental. Spatial-force
accuracy and auxiliary equilibrium need independent assessment; no universal
1D/2D/3D or shock validation, full KH/RT result, or <=2x GIZMO claim follows.
