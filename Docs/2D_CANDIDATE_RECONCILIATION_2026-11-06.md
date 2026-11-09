# 2D candidate history reconciliation

The November 6 shear_weight_rate_2d diagnostic reused the psi-times-fitted-strain R candidate from material_volume_candidate_audit.py. Its smooth-shear improvement is valid for the measured frozen states but does not qualify the method. The candidate was already retired in UNIFIED_HYDRO_EQUILIBRIUM_RETIREMENT_2026-11-05.md: the uniform-pressure acoustic frequency lower ratio is approximately 36, pressure force accuracy worsens on deformed geometry, and the native CPU coefficient screen failed. Do not launch this candidate or optimize its coefficient solver based on the recent small diagnostics. Results are preserved with selection_status.json linking the original equilibrium evidence.

The previous two diagnostic turns were missing this history check. They must not be presented as a new viable development path.

The adopted Kepler baseline uses W=0, fixed mass, no steering, no internal spin, native CUDA variational pressure, and volume-based adiabatic thermodynamics. Its trajectory combines RK4 and KDK restarts; throughput is documented as 3.34 times GIZMO runtime, not within the goal of 2. Its successful cold-orbit thermodynamics does not establish shock entropy generation.

The 1D material-face PLM/SSPRK2 nine-case suite remains a separate validated comparison. Fixed-topology moving-polygon 2D pulses failed convexity. Variable-mass FV tests do not satisfy the chosen fixed-mass scope. No complete 2D KH/RT qualification exists for a unified new method.

Next development must use the existing native 2D Laguerre geometry and preserve fixed mass, while addressing pressure work, geometry consistency and shock entropy without reviving the retired strain operator. Establish the viable operator definition and its pressure/volume adjoint before launching the full 2D suite. 3D remains deferred until 2D completion.
