# Saved Q2 late-state radial accuracy

The existing aligned Q2 fixed-polygon snapshots at t=30 and t=48.4 were re-read by original particle ID against the same saved t=0 state. The t=48.4 state continues the preserved t=48.3753 predecessor through a wall impact; it is not a fresh simulation. IDs, masses and W0 match. The reproducible [reader](audit_kepler_q2_late_core_20261029.py) and [full receipt](KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json) include hashes and distributions. All measures below are mass-weighted; position error is normalized by max(r0, epsilon) against the exact softened circular trajectory.

| Original region | Cells | t30 absolute density change | t48.4 absolute density change | t30 normalized position error | t48.4 normalized position error |
|---|---:|---:|---:|---:|---:|
| 0.1<=r0<0.5 | 1,476 | 19.56% | 37.57% | 2.45% | 10.92% |
| 0.5<=r0<0.75 | 3,931 | 0.136% | 4.01% | 0.00334% | 0.0754% |
| 0.5<=r0<2 | 49,038 | 0.0367% | 0.616% | 0.000493% | 0.00924% |

There are **no generator sites with r0<0.1** in this Q2 initial state, so this audit cannot report that innermost cohort or compare its accuracy with Q1's sampled core. The 0.1–0.5 region is low-density background, not the paper's dense disk, but it is part of the numerical domain and can affect neighboring cells. Disk mass-mean absolute log-entropy change is 3.14e-5 at t30 and 4.39e-5 at t48.4; the latter's signed mean is -2.96e-5, not evidence of net numerical heating. Total energy change is -2.40e-8 and -3.51e-8 relative to initial magnitude; global conservation does not imply spatial accuracy.

The substantial t48.4 inner-disk density change and large inner-background orbit/density errors mean this saved aligned Q2 trajectory should not be treated as an accurate long-time reference solely because it crosses the wall event or the full-disk mean is smaller. The densities are per-cell and have not here been decomposed into moving-mesh geometry versus physical transport; this comparison alone cannot determine the cause or attribute the difference to the one wall collision. Nor does it invalidate the separate Q1 trajectory. Boundary dependence, core resolution, circular-reference and early-accurate long-time validation remain required before a GIZMO comparison.
