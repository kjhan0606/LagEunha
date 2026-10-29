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

The [saved-state density comparison](KEPLER_Q2_INNER_DENSITY_T0_T30_T47P5_2026-10-29.png) shows t0, t30, and preimpact t47.5 using the same ±5% density-contrast scale relative to the initial mass-weighted disk mean. The nearest-generator images mask sites outside the original disk cohort; dark central area is the mask, not a vacuum claim. About 10.4% of visible t47.5 pixels saturate the scale; the image is not a cell-weighted error norm or a GIZMO figure comparison. [Settings and hashes](KEPLER_Q2_INNER_DENSITY_T0_T30_T47P5_2026-10-29.json), [rendering script](plot_kepler_q2_inner_density_history_20261029.py).

## Inner-neighbor counterfactual

In the saved t48.4 state, 294 of 3,931 inner-disk cells share a Voronoi face with at least one generator originating at r0<.5. These cells carry 22.18% of the weighted absolute inner-disk density-change norm. Holding all disk and exterior positions fixed while replacing only r0<.5 generator positions by the corresponding circular-orbit positions changes the full inner-disk mean absolute density change from 4.007784% to 4.007802%: almost no net change in this norm. At t30 the analogous values are 0.136263% and 0.136242%.

This reconstruction verifies actual selected areas against native areas within 5.84e-8 relative, and verifies that all actual and counterfactual selected polygons lie inside the polygon's inscribed circle. See [receipt](KEPLER_Q2_BOUNDED_INNER_SUPPORT_2026-10-29.json) and [script](audit_kepler_q2_inner_support_20261029.py). The result argues against direct, instantaneous inner-background support geometry as the explanation of the whole late inner-disk norm. It does **not** exclude perturbations propagated earlier from that background, pressure forces or temporal evolution; replacing positions at a single instant is not a physical evolution experiment. Net changes in an absolute norm can also mask cell-level changes.

## Inner-disk position counterfactual

At t48.4 the inner-disk generator position error relative to exact circular orbits is a mass mean **2.831% of sqrt(initial cell area)**, with cell p99 **23.77%**. The 0.0754% position error quoted above is relative to orbital radius and therefore hides the scale relevant to neighboring cells. Retaining every non-inner-disk position, but replacing the 3,931 inner-disk generator positions by exact circular positions, reduces their mean absolute density change from 4.0078% to 0.04424%. Replacing all sites initially inside r=.75 yields 0.04139%; extending replacement to r=.85 yields 0.04124%. The actual-area reconstruction matches native stored volumes, and all selected reconstructed cells are interior to the wall. See [receipt](KEPLER_Q2_BOUNDED_ORBIT_SUPPORT_2026-10-29.json) and [script](audit_kepler_q2_disk_orbit_support_20261029.py).

This establishes that the saved inner-disk geometry is highly sensitive to the inner-disk generators' position deviations at this instant. It does not establish what caused those deviations: gravity timestep, pressure forces, earlier background/boundary effects or other mechanisms remain possible. The counterfactual changes positions without evolving velocity or energy, so it is not an alternate fluid solution or a numerical fix.

## Gravity-only continuation from the saved t30 state

The saved t30 native positions and velocities for all 3,931 inner-disk sites were separately evolved under only the same softened central gravity to the saved t=48.37534965270793 **preimpact** state. The old native log provides 11,579 printed timesteps, but each is rounded. Their sum exceeds the exact checkpoint interval by 7.08e-8; the independent replay rescales all printed steps by 0.9999999961474437 to match that exact interval. The result is therefore an approximate timestep replay, not an exact native-stage reconstruction. Source hashes and full details are in the [receipt](KEPLER_Q2_LATE_GRAVITY_CONTROL_2026-10-29.json) and [script](audit_kepler_q2_late_gravity_control_20261029.py).

The gravity-only endpoint differs from the actual native preimpact endpoint by a mass-mean position distance of **3.754e-4** (cell p99 **3.117e-3**). Repeating the independent gravity integration with two RK4 substeps for every printed step changes the result by mass mean **1.488e-10** (p99 **4.373e-10**). Thus the observed native-versus-gravity-only discrepancy is much larger than the measured one-versus-two-substep gravity sensitivity in this specific interval. The result does not give a rigorous gravity error bound and does not attribute all the discrepancy to pressure: the native calculation includes coupled hydrodynamics and the replay has rounded/approximately aligned timesteps. It also does not identify when the trajectory first departed or prove that the t30 state was accurate.

## Instantaneous pressure-force reconstruction

An independent W0 variational face reconstruction was applied to the same t30 and t48.4 saved states for all sites originally at r<2. Maximum area discrepancies from stored native areas are below 6e-8. In the inner disk, the mass-mean magnitude of reconstructed pressure acceleration rises from **2.15e-7** at t30 to **7.26e-6** at t48.4, a factor of about 33.8; cell p99 rises from **1.43e-6** to **5.34e-5**. Its mass-mean acceleration ratio relative to central gravity rises from 7.74e-8 to 2.59e-6. The mean instantaneous absolute compressional-work rate `|Udot|/(U Omega)` changes from 0.00525 to 0.01612; it is not entropy production. The complete radial statistics and source hashes are in [the receipt](KEPLER_Q2_LATE_PRESSURE_2026-10-29.json) and [reader](audit_kepler_q2_late_pressure_20261029.py).

This is a reconstruction from two endpoints, not accumulated impulse or a causal attribution for the orbit discrepancy. The reconstruction method previously matched eight Q1 native RHS samples at t10, but no direct Q2 native RHS comparison was made here; actual Q2 forces and their time history require separate checks. The pressure acceleration is small relative to gravity at both endpoints and may vary significantly between them.

## Existing intermediate checkpoints

Seven saved checkpoints from the original t30→60 attempt fill in the interval without rerunning the simulation. Their exact times are read from the checkpoint footers; all retain the original IDs, fixed masses and W0. The [history reader](audit_kepler_q2_pressure_history_20261029.py) applies the same independent pressure reconstruction, compares actual areas to stored areas and records all hashes in [the receipt](KEPLER_Q2_PRESSURE_HISTORY_2026-10-29.json). The rows below give mass means for the original .5<=r0<.75 inner disk; displacement is normalized by sqrt(initial cell area).

| t | Absolute density change | Position error / cell scale | Pressure acceleration |
|---:|---:|---:|---:|
| 30.000 | 0.136% | 0.127% | 2.15e-7 |
| 32.502 | 0.202% | 0.177% | 3.37e-7 |
| 35.001 | 0.313% | 0.258% | 5.34e-7 |
| 37.500 | 0.489% | 0.391% | 8.35e-7 |
| 40.001 | 0.780% | 0.603% | 1.37e-6 |
| 42.501 | 1.229% | 0.943% | 2.13e-6 |
| 45.001 | 2.023% | 1.492% | 3.57e-6 |
| 47.500 | 3.301% | 2.385% | 5.92e-6 |
| 48.400 | 4.008% | 2.831% | 7.26e-6 |

These sampled quantities rise together throughout this interval, and the large change is already present at t47.5, before the single wall impact at t48.377. This excludes that specific final impact as the sole origin of the preimpact change. It does not exclude earlier wall effects or establish whether pressure perturbations drive the geometric displacement, result from it, or both. Reconstructed pressure forces still need direct native Q2 RHS validation, and snapshots 2.5 time units apart do not resolve all intermediate behavior.

## Nearest-generator pairs

For each original inner-disk site, select its nearest generator in the actual saved state and compare the distance to the **same pair's** distance when both original positions follow circular orbits. At t30 and t40 none of the 3,931 sites have a ratio below 0.8; at t47.5, **244 sites** do, and they contribute **33.4%** of the inner-disk weighted absolute density-change norm. At t48.4, **312 sites** do and contribute **39.7%**. The cell p1 distance ratio is 0.985, 0.914, 0.652, 0.585 at t30, 40, 47.5, 48.4, respectively. These counts are sites, not independent pairs; mutual nearest neighbors can be counted twice.

The t47.5 nearest-neighbor radial relative speed has mass mean -6.29e-5, while the same pairs on the circular reference have mean about 9.6e-18. At t30 the corresponding actual mean is -2.11e-6. Some local pair convergence is normal in shear, and selection of a nearest pair is not a direct viscous-force or shock measurement. Together these observations identify a growing close-neighbor component in this aligned Q2 case, not its cause or a universal pairing instability. They do not establish that AV supplies the force, or that internal cell spin will cure it. See [receipt](KEPLER_Q2_INNER_PAIRS_2026-10-29.json) and [reader](audit_kepler_q2_inner_pairs_20261029.py).

## Reconstructed pressure direction for close pairs

At t47.5003 all 244 close-neighbor **sites** identified above have a nearest peer on the same original radial ring, within 1e-10 in initial radius. Their separation is predominantly azimuthal: median absolute projection of the separation unit vector on the radial direction is 0.0232, and none exceed 1/sqrt(2). All 244 have negative relative speed along their separation line, with mass mean -5.95e-4. For the same 244 sites, the independently reconstructed **total** W0 pressure accelerations give mass-mean `(a_pressure,j-a_pressure,i) dot n_ij = +4.24e-5`; all 244 projected differences are positive. Positive means the pressure contribution to the second derivative of separation acts in the separating direction at that instant. At t48.4 all 312 close sites likewise approach while all reconstructed pressure differences have positive projection (mean +4.93e-5).

The reconstruction includes each site's pressure contributions from all faces; it does not isolate the shared face. Gravity and transverse relative motion also contribute to the full second derivative of pair separation. This is not a direct native Q2 RHS output, a time-integrated force, evidence about artificial-viscosity forces, or proof of the instability's cause. It argues against a claim that the reconstructed **net pressure acceleration at these snapshots** directly pulls these close pairs together. It also does not imply these results hold for the earlier HLLC/AV configurations or for other initial meshes. Full definitions and source hashes are in [the receipt](KEPLER_Q2_PAIR_PRESSURE_2026-10-29.json) and [reader](audit_kepler_q2_pair_pressure_20261029.py).

## Why approaching despite pressure response?

For separation d=|x_j-x_i|, the instantaneous identity is `d_ddot = n dot (a_j-a_i) + |v_rel - n(n dot v_rel)|^2/d`. For the same 244 close-neighbor sites at t47.5003, the mass-mean contributions reconstructed from the saved state are:

| Contribution to separation acceleration | Mass mean |
|---|---:|
| Central gravity projected along pair | -4.57455e-2 |
| Transverse relative-velocity curvature | +4.55791e-2 |
| Gravity plus curvature | -1.66417e-4 |
| Reconstructed total pressure-acceleration difference | +4.24079e-5 |
| Sum of these three terms | -1.24010e-4 |

The gravity and velocity-curvature terms nearly cancel, but their residual is negative. The positive pressure term offsets about a quarter of that residual at this saved instant. The modeled sum is negative for all 244 selected sites; all are already approaching. As a check, applying the gravity-plus-curvature expression to the same IDs on exact circular orbits yields a mass-mean absolute separation acceleration of 1.90e-16. The result describes the instantaneous stored state, not why the velocity or position perturbation appeared. Pressure may have helped generate the earlier velocity deviations even though its current projection resists them. The reconstructed pressure has not been checked against native Q2 RHS, and the sum does not claim to include any other native acceleration terms. See [receipt](KEPLER_Q2_PAIR_KINEMATICS_2026-10-29.json) and [script](audit_kepler_q2_pair_kinematics_20261029.py).

## Native CUDA check on three close pairs

The six cells of the three closest distinct pairs were then checked with an actual native CUDA K1 evaluation on Ada, using the untouched saved t47.50031653561342 checkpoint and the already-tested wall-event v2 executable. The latter is not the historical Q2 executable, so this directly validates this candidate on the stored state, not every historical execution path. The diagnostic log explicitly reports no accepted evolution. The independent polygon reconstruction matches the six native pressure accelerations to at worst 6.51e-18 absolute; its pressure-work comparison and geometry checks also pass. See [native audit](KEPLER_Q2_NATIVE_PAIR_RHS_2026-10-29/AUDIT.json), [input](KEPLER_Q2_NATIVE_PAIR_RHS_2026-10-29/INPUT.json), and [staging script](stage_kepler_q2_pair_rhs_20261029.py).

| IDs | Relative speed along pair | Native pressure relative acceleration along pair |
|---|---:|---:|
| 31088, 30249 | -1.47331e-3 | +1.26833e-4 |
| 29637, 23949 | -1.45954e-3 | +1.25256e-4 |
| 25319, 23161 | -1.06536e-3 | +9.03788e-5 |

The [projection receipt](KEPLER_Q2_NATIVE_PAIR_RHS_2026-10-29/PAIR_PROJECTION.json) confirms that these selected pairs approach while native total pressure accelerations oppose their approach at that instant. The six selected cells are an extreme sample, not a population-wide native RHS validation. The test neither establishes the perturbation's cause nor checks the time history, AV behavior in older configurations, intrinsic spin, or long-time stability.

## Radial separation and Kepler shear

For the same 244 close sites at t47.5003, each selected nearest peer originated on the same radial ring. Their present radial separation has mass-weighted RMS **1.7024e-4**. The actual angular-speed difference has RMS **1.15438e-3**; the difference predicted by the local softened circular frequency at each site's **current** radius has RMS **1.15755e-3**. The RMS residual is **3.26e-6**, and the uncentered weighted correlation is **0.99999978**; all selected differences have the same sign as the local circular-shear prediction. The t48.4 sample shows the same pattern. [Receipt](KEPLER_Q2_CLOSE_SHEAR_2026-10-29.json); [reader](audit_kepler_q2_close_shear_20261029.py).

This shows that the close-pair angular-speed difference at these snapshots is consistent with the Kepler shear induced by their small current radial offsets, rather than requiring a large noncircular angular velocity at that instant. The selection depends on close separation and mutual neighbors can be counted twice. This is an instantaneous correlation: it does not determine the origin of the radial offsets, prove their prior growth mechanism, or show that changing spin, weights or force law cures them. Any proposed remedy should track radial offsets, angular-speed differences and pair separation together through time, while preserving the no-steering requirement.

## Fixed-pair angular-momentum history

The same 244 late-selected sites and their t47.5 nearest peers were tracked backwards by ID through the existing t30–47.5 checkpoints, without reselecting neighbors. Their pairwise radial difference has mass-weighted RMS 6.12e-6 at t30 and 1.702e-4 at t47.5. Pairwise specific angular-momentum difference grows from 4.12e-6 to 1.154e-4. The circular-orbit relation between the actual angular-momentum difference and radial separation has an RMS radial residual 1.56e-6 at t47.5, much smaller than the 1.702e-4 radial difference. This is consistent with a growing angular-momentum split and associated radius split; it does not say which preceded the other.

At each checkpoint, total W0 pressure force was independently reconstructed and the difference of the two specific pressure torques was calculated. A trapezoid sum over the eight snapshots yields RMS integrated pair-torque difference **1.14775e-4** from t30 to t47.5, versus RMS measured specific angular-momentum difference change **1.11309e-4**; RMS of the per-pair residual is **3.74e-6**. The analytic central gravity is torque-free about the center. The approximate agreement supports pressure torque as a candidate mechanism, but the 2.5-time-unit sampling can alias orbital-scale variation; the quadrature and Q2 native RHS must be checked at finer temporal resolution before claiming a precise torque budget. The t47.5 native diagnostic validated only three selected pairs with the wall-event v2 executable. See [receipt](KEPLER_Q2_CLOSE_MOMENTUM_HISTORY_2026-10-29.json) and [reader](audit_kepler_q2_close_momentum_history_20261029.py).

For the same fixed late-selected sites, at every one of the eight sampled times the reconstructed torque difference has the **same sign** as the pair's specific angular-momentum difference at all 244 sites. The mass-weighted projection `mean(DeltaL*DeltaTau)/mean(DeltaL^2)` ranges from 0.183 to 0.213 per simulation-time unit; positive values mean this pressure contribution instantaneously increases the squared L split. The observed RMS split grows by about 28 times over the same 17.5 units, corresponding to an average logarithmic rate about 0.190 per unit. [Receipt](KEPLER_Q2_TORQUE_FEEDBACK_2026-10-29.json); [reader](audit_kepler_q2_torque_feedback_20261029.py).

This is stronger evidence of a pressure-torque feedback in this particular aligned Q2 trajectory than a mere correlation of the global norms, but it is still not a linear stability proof or a resolved time integral. The pairs were selected because they were close at the end, samples are sparse, individual native Q2 force agreement has only been checked on three pairs, and pressure force, geometry and velocity evolve together. Positive projected torque along the **angular-momentum split** does not contradict a repulsive pressure acceleration along the **instantaneous pair separation**; they are different projections. A remedy must break the coupled feedback without steering, while preserving the relevant conservation laws and testing other meshes and boundary conditions.

## Acoustic-scale comparison with the earlier shear-cell investigation

For those same 244 late-selected sites, the mass mean of `sqrt(gamma P/rho)/sqrt(cell area)` rises from 0.08325 at t30 to 0.09481 at t47.5. The instantaneous pair-torque projection divided by this mean is 2.17–2.27 across the eight saved times. This is an order-of-magnitude comparison, not a measured mode growth rate. The earlier investigation summarized in [the October 9 handover](HANDOVER_2026-10-09.md) reported about `1.1 c_s/d` for a high-pressure shear-cell mode; its definition of `d`, state, perturbation and growth estimator differ from the quantities above. The earlier detailed verification files referenced there are not available in this checkout, so their precise normalization cannot presently be reproduced. Comparable scales alone do not establish that Q2 exhibits the same instability, nor that acoustic scaling holds under controlled changes of pressure or resolution. See [receipt](KEPLER_Q2_FEEDBACK_SOUND_SCALE_2026-10-29.json) and [reader](audit_kepler_q2_feedback_sound_scale_20261029.py).

## Inner-ring angular-momentum spectrum

Independently of the late-selected pairs, all 3,931 original inner-disk sites in 11 initially uniform rings were sorted by their initial azimuth. For each ring, a cyclic discrete Fourier transform was applied to the per-site change in specific angular momentum after subtracting that ring's mean change. Summing the resulting powers with each ring's total mass as its weight, the fraction in `0.2 <= m/N < 0.3` rises from **44.3% at t30** to **68.7% at t47.5**; the fraction in `0.4 <= m/N <= 0.5` is **0.135% at t47.5**. At that endpoint each ring's largest individual mode has `m/N` between approximately 0.20 and 0.28. The ring-index DFT uses the closed angular coordinate; it does **not** apply periodic spatial boundary conditions to this bounded Q2 calculation. See [receipt](KEPLER_Q2_RING_SPECTRUM_2026-10-29.json) and [reader](audit_kepler_q2_ring_spectrum_20261029.py).

The earlier handover's `m≈0.4N` referred to tangential displacement in different initial conditions, while this is the angular-momentum change on Q2 material rings. Therefore their difference does not rule out a shared underlying mechanism. These eight sparse snapshots establish neither a mode's temporal growth rate nor a causal relationship; the Fourier amplitudes, ring population and radial domain are also different from the previous 244-site torque sample. A matched displacement-and-momentum mode analysis and controlled pressure/resolution comparison are needed before designing a correction around a particular wavelength.

The matched tangential-displacement analysis now uses the same 11 rings and snapshots. It measures each site's angular offset from its same-time analytic softened circular orbit, multiplies by its original radius, and removes the per-ring mean. At t47.5, **63.3%** of the ring-mass-weighted displacement power is in `0.2 <= m/N < 0.3`, and **0.037%** is in `0.4 <= m/N <= 0.5`. The largest individual displacement mode has the same integer `m` as the angular-momentum residual in all 11 rings at this checkpoint. At t30 the two largest modes agree in 5 of 11 rings. See the [receipt](KEPLER_Q2_RING_DISPLACEMENT_2026-10-29.json) and [reader](audit_kepler_q2_ring_displacement_20261029.py).

This is a comparison of spectra, not a demonstrated common growth rate or pressure-driven causal mechanism. The old `m≈0.4N` result used a different experiment; the individual Q2 peak modes are integer-valued and their normalization differs across rings. Both the ring-index transforms and analytic circular reference are post-processing of a fixed-polygon run, not a periodic physical boundary. No simulation was restarted for this comparison.

## Pressure-torque projection onto the full-ring modes

The independently reconstructed W0 pressure force was projected onto the
mean-subtracted specific-angular-momentum error of **all** original sites in
the same 11 inner rings, rather than selecting close pairs. For the band
`0.2 <= m/N < 0.3`, the ring-mass-weighted projection
`sum Re(conj(L_k) tau_k) / sum |L_k|^2` is **0.190–0.205** per unit time
across the eight snapshots from t30 to t47.5. The band power increases by a
factor of **894.3** over that interval; the corresponding average logarithmic
amplitude rate is **0.1942** per unit time. Reconstructed cell areas match
stored volumes to relative error below 5.95e-8. See the [receipt](KEPLER_Q2_MODE_TORQUE_2026-10-29.json) and [reader](audit_kepler_q2_mode_torque_20261029.py).

The projection is the instantaneous contribution of reconstructed pressure
torque to the amplitude's logarithmic rate, not physical energy injection.
Its agreement with the interval-averaged rate is evidence consistent with
pressure-torque feedback across the whole inner ring population, beyond the
late-selected close pairs. The snapshots are sparse, the force is reconstructed
and only three extreme pairs have a native RHS cross-check, and the coupled
trajectory includes other possible terms. This does not prove causality or
show a remedy. A pressure/resolution-controlled native comparison and a
conservative operator that changes this feedback without steering remain
necessary. Weak-power bands in the receipt can have large, noisy normalized
rates and should not be ranked by rate alone.

## Frozen-state normal-traction counterfactual

On the same eight saved geometries, replace the reconstructed W0 force in
post-processing by a single equal-and-opposite normal traction per interior
face, using the arithmetic mean of its two stored pressures as face pressure.
No native trajectory or snapshot is modified. In the main `0.2 <= m/N < 0.3`
band its torque projection is **0.1901–0.2047** per time unit; at t47.5 it
is **0.204713**, versus **0.204999** from the existing variational-force
reconstruction. The eight-snapshot difference is at most 0.000286 per time
unit. See [receipt](KEPLER_Q2_NORMAL_TRACTION_2026-10-29.json) and
[reader](audit_kepler_q2_normal_traction_20261029.py).

Thus, on these frozen states, replacing the existing force by this particular
mean-pressure normal traction does not materially alter the projection driving
the main angular-momentum error. This is a negative screening result, not a
stability test: an actual Riemann pressure, evolved spin, mass-centroid
closure, energy/pressure work and future trajectories are not represented.
Equal-and-opposite force gives pair momentum balance but is not by itself a
complete energy- and spin-conservative solver. The outcome does not exclude a
different traction or evolved affine/spin model; it argues against spending a
long GPU run on this simplified substitution alone.

## Candidate normal-pressure torque about cell centroids

On the t30 and t47.5 saved states, each of the 3,931 original inner-disk
polygons was independently reconstructed. Assuming uniform cell density,
its area centroid and polar mass moment define a provisional mass centroid
`C` and inertia `I`. For the same arithmetic-mean shared normal traction used
above, add `(face centroid - C) cross force` over all faces of each cell. The
sum of absolute **net** cell torques divided by the sum of absolute individual
face torques is **0.00090 at t30** and **0.0241 at t47.5**: most face torques
cancel, although the residual grows. The mass-mean of the instantaneous
scale `|tau|/(I Omega^2)` is **9.95e-8** and **4.14e-6** respectively;
the t47.5 cell p99 is **3.33e-5**. For faces between selected cells, the
orbital-plus-intrinsic pair torque cancellation closes to relative error
4.01e-16 at t47.5. Reconstructed areas agree with stored volumes to better
than 5.88e-8. See [receipt](KEPLER_Q2_PRESSURE_SPIN_2026-10-29.json) and
[reader](audit_kepler_q2_pressure_spin_20261029.py).

This is an instantaneous hypothetical traction calculation, not an evolved
spin measurement. The mass centroid and inertia are not native evolved state;
uniform density is an assumption. `tau/I` alone omits `I_dot`, affine strain,
spin transport, pressure work and the choice of Riemann state. Its small value
relative to `Omega^2` does not establish that rotation is irrelevant over many
orbits or that a spin model would fail; equally, nonzero torque does not
establish that adding spin cures the measured error. An integrated conservative
spin/affine operator remains necessary to answer that question.

## Face-normal velocity diagnostic for the fixed-mass design

For interior faces whose two sites originated in the inner disk, the exact
instantaneous Voronoi geometric normal speed was compared with the arithmetic
mean of the two saved-gradient linear velocity reconstructions at the face
centroid. At t30, the face-length-weighted mean absolute difference is
**0.286 times the mean local sound speed**; at t47.5 it is **0.395**. At the
latter time, **6.51% of the included face length** has a difference greater
than that sound-speed scale. There are 11,078 selected faces at each saved
time. The face-speed formula differentiates the perpendicular-bisector
equation and includes the off-center face correction. See [receipt](KEPLER_Q2_FACE_SLIP_2026-10-29.json)
and [reader](audit_kepler_q2_face_slip_20261029.py).

This compares geometric speed with a simple reconstruction average: it is
not the actual Riemann contact speed or measured native mass flux, and saved
velocity gradients can come from a nearby RK stage. It cannot establish that
the existing fixed-mass solver violates its own equations. For an eventual
spin/ALE design it does show why assuming that generator-defined face speed
automatically equals reconstructed fluid speed is unsafe. The W0 fixed-mass
limit must explicitly define its interface motion and flux closure, including
energy and moments, while keeping the no-steering condition.
