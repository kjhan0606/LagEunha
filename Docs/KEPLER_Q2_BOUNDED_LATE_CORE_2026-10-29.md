# Saved Q2 late-state radial accuracy

The existing aligned Q2 fixed-polygon snapshots at t=30 and t=48.4 were re-read by original particle ID against the same saved t=0 state. The t=48.4 state continues the preserved t=48.3753 predecessor through a wall impact; it is not a fresh simulation. IDs, masses and W0 match. The reproducible [reader](audit_kepler_q2_late_core_20261029.py) and [full receipt](KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json) include hashes and distributions. All measures below are mass-weighted; position error is normalized by max(r0, epsilon) against the exact softened circular trajectory.

| Original region | Cells | t30 absolute density change | t48.4 absolute density change | t30 normalized position error | t48.4 normalized position error |
|---|---:|---:|---:|---:|---:|
| 0.1<=r0<0.5 | 1,476 | 19.56% | 37.57% | 2.45% | 10.92% |
| 0.5<=r0<0.75 | 3,931 | 0.136% | 4.01% | 0.00334% | 0.0754% |
| 0.5<=r0<2 | 49,038 | 0.0367% | 0.616% | 0.000493% | 0.00924% |

There are **no generator sites with r0<0.1** in this Q2 initial state, so this audit cannot report that innermost cohort or compare its accuracy with Q1's sampled core. The 0.1–0.5 region is low-density background, not the paper's dense disk, but it is part of the numerical domain and can affect neighboring cells. Disk mass-mean absolute log-entropy change is 3.14e-5 at t30 and 4.39e-5 at t48.4; the latter's signed mean is -2.96e-5, not evidence of net numerical heating. Total energy change is -2.40e-8 and -3.51e-8 relative to initial magnitude; global conservation does not imply spatial accuracy.

The substantial t48.4 inner-disk density change and large inner-background orbit/density errors mean this saved aligned Q2 trajectory should not be treated as an accurate long-time reference solely because it crosses the wall event or the full-disk mean is smaller. The densities are per-cell and have not here been decomposed into moving-mesh geometry versus physical transport; this comparison alone cannot determine the cause or attribute the difference to the one wall collision. Nor does it invalidate the separate Q1 trajectory. Boundary dependence, core resolution, circular-reference and early-accurate long-time validation remain required before a GIZMO comparison.

## Circular-orbit geometric comparison

An independent all-cell bounded Voronoi reconstruction was subsequently performed on both saved states and on configurations formed by placing the original Q2 sites on exact softened circular orbits at the same times. The reconstructed actual cell areas agree with stored native areas to at worst 5.96e-8 relative. Two guard radii agree and each construction tiles the polygon; see [geometry receipt](KEPLER_Q2_BOUNDED_LATE_GEOMETRY_2026-10-29.json) and [reader](audit_kepler_q2_late_geometry_20261029.py).

| Region | t | Measured density change | Exact-circular geometry change | Actual density relative to circular geometry |
|---|---:|---:|---:|---:|
| Inner disk r0=.5–.75 | 30 | 0.1363% | 0.0407% | 0.1184% |
| Inner disk r0=.5–.75 | 48.4 | 4.0078% | 0.0412% | 4.0055% |
| Full disk r0=.5–2 | 48.4 | 0.6163% | 0.0225% | 0.6128% |

Entries are mass means of absolute per-cell fractional changes; the columns are ratios with different denominators and are **not additive**. The reference is a kinematic counterfactual, not a self-consistent evolved fluid solution. Nevertheless, pure circular-motion Voronoi breathing is much too small to account for the late inner-disk difference. This does not independently establish whether the deviation is caused by inner-background dynamics, boundary effects, force/work discretization, or another mechanism. It also does not quantify entropy production.
