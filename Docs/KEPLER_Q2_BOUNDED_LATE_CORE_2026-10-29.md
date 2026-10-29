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

## Inner-neighbor counterfactual

In the saved t48.4 state, 294 of 3,931 inner-disk cells share a Voronoi face with at least one generator originating at r0<.5. These cells carry 22.18% of the weighted absolute inner-disk density-change norm. Holding all disk and exterior positions fixed while replacing only r0<.5 generator positions by the corresponding circular-orbit positions changes the full inner-disk mean absolute density change from 4.007784% to 4.007802%: almost no net change in this norm. At t30 the analogous values are 0.136263% and 0.136242%.

This reconstruction verifies actual selected areas against native areas within 5.84e-8 relative, and verifies that all actual and counterfactual selected polygons lie inside the polygon's inscribed circle. See [receipt](KEPLER_Q2_BOUNDED_INNER_SUPPORT_2026-10-29.json) and [script](audit_kepler_q2_inner_support_20261029.py). The result argues against direct, instantaneous inner-background support geometry as the explanation of the whole late inner-disk norm. It does **not** exclude perturbations propagated earlier from that background, pressure forces or temporal evolution; replacing positions at a single instant is not a physical evolution experiment. Net changes in an absolute norm can also mask cell-level changes.

## Inner-disk position counterfactual

At t48.4 the inner-disk generator position error relative to exact circular orbits is a mass mean **2.831% of sqrt(initial cell area)**, with cell p99 **23.77%**. The 0.0754% position error quoted above is relative to orbital radius and therefore hides the scale relevant to neighboring cells. Retaining every non-inner-disk position, but replacing the 3,931 inner-disk generator positions by exact circular positions, reduces their mean absolute density change from 4.0078% to 0.04424%. Replacing all sites initially inside r=.75 yields 0.04139%; extending replacement to r=.85 yields 0.04124%. The actual-area reconstruction matches native stored volumes, and all selected reconstructed cells are interior to the wall. See [receipt](KEPLER_Q2_BOUNDED_ORBIT_SUPPORT_2026-10-29.json) and [script](audit_kepler_q2_disk_orbit_support_20261029.py).

This establishes that the saved inner-disk geometry is highly sensitive to the inner-disk generators' position deviations at this instant. It does not establish what caused those deviations: gravity timestep, pressure forces, earlier background/boundary effects or other mechanisms remain possible. The counterfactual changes positions without evolving velocity or energy, so it is not an alternate fluid solution or a numerical fix.

## Gravity-only continuation from the saved t30 state

The saved t30 native positions and velocities for all 3,931 inner-disk sites were separately evolved under only the same softened central gravity to the saved t=48.37534965270793 **preimpact** state. The old native log provides 11,579 printed timesteps, but each is rounded. Their sum exceeds the exact checkpoint interval by 7.08e-8; the independent replay rescales all printed steps by 0.9999999961474437 to match that exact interval. The result is therefore an approximate timestep replay, not an exact native-stage reconstruction. Source hashes and full details are in the [receipt](KEPLER_Q2_LATE_GRAVITY_CONTROL_2026-10-29.json) and [script](audit_kepler_q2_late_gravity_control_20261029.py).

The gravity-only endpoint differs from the actual native preimpact endpoint by a mass-mean position distance of **3.754e-4** (cell p99 **3.117e-3**). Repeating the independent gravity integration with two RK4 substeps for every printed step changes the result by mass mean **1.488e-10** (p99 **4.373e-10**). Thus the observed native-versus-gravity-only discrepancy is much larger than the measured one-versus-two-substep gravity sensitivity in this specific interval. The result does not give a rigorous gravity error bound and does not attribute all the discrepancy to pressure: the native calculation includes coupled hydrodynamics and the replay has rounded/approximately aligned timesteps. It also does not identify when the trajectory first departed or prove that the t30 state was accurate.
