# RT128 reflected-wall independent reference — 2026-11-09

The native128x256 trajectory uses reflecting physical walls at y=0,1 and evolves every owned cell. The published GIZMO RT picture instead fixes particles initially at y<.1/y>.9. These different boundary conditions cannot establish a matched quantitative comparison.

## Reused reference implementation

The existing `Research/MaterialVolume_20261105/gizmo_rt_mirror_20261106_v1` source already implements a mirrored periodic reference: physical half0<y<1 has gravity−.5; reflected half1<y<2 has gravity+.5. Equal scalar states and reflected vertical velocity represent the reflecting planes through symmetry. Existing restart timeline overflow protection is retained. This is not asserted to be discrete operator equivalence with native ghost construction; mirror symmetry is measured.

The new isolated `common_rt128_gizmo_mirror_20261109_v1` copies that source. Only Config.sh changes: BOX_LONG_X=.5, with BOX_LONG_Y=2 and BoxSize=1 unchanged. No hydro, gravity, limiter or native production code changes. Physical128x256=32768 particles plus32768 mirrors are evolved. Physical input position, velocity, mass and specific internal energy match the existing current-native IC bitwise. A saved mapping relates GIZMO positive IDs to native IDs. Gamma is nativefloat1.4. No pressure-normalizing U adjustment.

The copied source builds successfully; binary SHA256 f7c1e4acbee4fb7c2c3821d9a3ecc85b9490960caaeb67c2facaed8268bc7f52. Build commands/logs and the complete initial/source manifest are retained. Reference CPU24,25 with2 OpenMP threads on LagEunha; no GPU allocation or changes to ongoing native runs. Mirror particles double computed count, so these timings cannot qualify the <=2x matched-accuracy performance requirement.

## Preflight and continuation

The actual t=.002 preflight completed in16 steps and7.9552s. All masses remain exact and density/U positive. Mirror position and velocity residuals are1.11e−15 and1.28e−13; maximum relative density asymmetry7.50e−14. Physical particles stay strictly between the walls. The initial GIZMO kernel-density initialization produces a mass-weighted relative pressure discrepancy0.0008319663 (0.0831966%) from native, maximum0.0008337206. This is recorded, not silently corrected or described as identical initialized pressure. Initial diagnostic-snapshot velocities include a possible first half kick; parity is certified against input HDF5 rather than that snapshot.

The detached runner continues only from verified successful restarts to actual targets1.2,2.4,4. Existing stop-condition/timeline conventions use TimeMax half a maximum step below subsequent targets, then explicitly require actual final snapshot time to match within1e−10. An actual-time mismatch or mirror/positivity failure stops the batch without restarting or changing physics. No scheduled endpoint snapshot is used; final forced output follows the second kick. Mass, finite positive state, physical-half occupancy, mirrored position/velocity/rho/U residuals and actual endpoint time are checked at every stop. Mirror numerical tolerances1e−8 are symmetry checks, not physical-accuracy criteria.

The output plan allows at most8 snapshots, estimated12MiB each, plus512MiB restart-file allowance. The preflight output, including restart state, is56,194,944bytes. Free space before launch149,076,251,246,592bytes. Effective parameters/environment/command/lineage are saved per segment. Initial preflight checkpoints are reused, not discarded and rerun. Full reference morphology is pending; no accuracy or speed certificate follows from this preflight.

## Saved native boundary-cohort measurement

`common_rt128_20261109_v1/boundary_cohorts_t15_v2/report.json` follows initial material IDs through existing snapshots, without new hydro. At t=1.5, each of the y<.1/y>.9 cohorts contains3328 cells. The lower cohort has displacement RMS0.000908926 and speed RMS0.0135637; the upper cohort has displacement RMS0.000127125 and speed RMS0.00509711. Maximum displacements are0.00161981 and0.000423099. These are nonzero, unlike fixed-coordinate layers. This confirms that the published frozen-layer setup cannot be substituted for this reflecting-wall trajectory. It does not prove that the active-wall solution is wrong, or that the boundary difference causes any interior morphology error. The initial script used the wrong result field name and stopped before writing a report; its source remains in boundary_cohorts_t15. The v2 uses and verifies saved_end=last_accepted_time=target.

No3D evolution or second-paper drafting is started. Current 1D/2D/3D method consistency and matched-accuracy performance remain incomplete.


## Same-time saved-state comparison pipeline

The CPU-only comparison worker531641 waits for verified native/GIZMO endpoints1.2,2.4,4 and never starts or restarts hydro. The actual live process and command were verified after launch. A12-hour finite wait checks each dependency supervisor; missing handles or failed results stop the analysis. Every comparison requires matched actual time, exact material masses/ID mapping, initial-checkpoint hashes and the recorded final GIZMO snapshot hash. Native geometry is independently rebuilt and certified with existing rt_power_polygons.reconstruct.

At the completed t=.002 preflight, the rendered actual-native-polygon versus GIZMO-nearest-density image was visually inspected: axes touch, shared.8..2.8 color range, no clipped cells, titles inside the panels, density reconstruction types explicitly labelled. These are initial interface profiles, not evolved RT mushroom results. Native density range.99999687–2.00031403, GIZMO.99916071–1.99863731. No density-error norm is formed from the different renderings. The original preflight analysis source is preserved next to its results; later endpoint code additionally records the energy sum.

Common point-mass observables are computed at native generator/GIZMO particle positions, without interpolation: mass and x-momentum use cos(m*pi*y)*exp(-2*pi*i*k*x/.5), vertical momentum uses sin(m*pi*y) with k=0..4,m=0..8. Full, nonzero-x and seeded k=2 mode subsets are retained with numerator and reference norms. They are point-mass statistics, not finite-volume Fourier integrals or exact-solution errors. At t=.002 the seeded vertical-momentum relative difference is1.35385e−5; horizontal momentum is initially very small, so its relative difference2.24205% must be read with the absolute difference9.44049e−8 per total mass.

Independent physical-half point-potential energy sums use sum m*(u+|v|²/2+.5*y). At t=.002 the native relative change is−1.12816e−16, consistent with its existing ledger−1.14272e−16; GIZMO is+1.66469e−13. No morphology/late-time accuracy or speed claim follows. Source/evidence: common_rt128_gizmo_mirror_20261109_v1/compare.py, comparison/preflight/{run.py,report.json,energy.json,density.png}, comparison/launcher.json and status.json.


## Matched t=1.2 reference endpoint checked

The automatic comparison produced `common_rt128_gizmo_mirror_20261109_v1/comparison/t12/report.json` and density.png/pdf; the figure was visually inspected. Both physical domains contain32768cells/particles (128×256); GIZMO additionally evolves a mirror half,65536computed particles. At t=1.2, the primary RT wave shapes are close. Native/GIZMO vertical-velocity RMS=.03397297/.03382905 (difference about.43%). Point-mass low-mode differences are2.02% for y momentum,5.84% for x momentum; seed-k2 y-momentum difference1.46%. These are defined spectral moments, not image density error. Initial kernel-pressure mismatch=.0832% persists as a comparison caveat. Native polygons and nearest GIZMO density samples use different display representations; neither image is an exact solution. Total energy relative changes−4.69e-15/native and+4.90e-11/GIZMO use particle-point gravitational potential.

Actual native616677 and GIZMO613490 processes were confirmed live advancing toward t=2.4; comparison follower531641 is alive. These simulations use the frozen current model, **not the new experimental face-local scalar candidate**. Shared host timing and GIZMO's additional mirror particles do not establish a matched-cost speed ratio.

## t=2.4 mirror certificate failure and separate diagnostics

The original GIZMO t24 trajectory exited0 at exactly2.4. Mass/ID, finite and
positive states, physical-half containment all pass. The predeclared1e-8
mirror certificate fails: position max1.0538184e-9 passes, velocity max
2.1569889e-8, density relative max1.6703223e-8, specific-U relative max
1.2115922e-7 fail. Evolution supervisor512343 and comparison531641 correctly
terminated; both and child613490 were confirmed absent. Original result and
failed statuses remain unchanged. No t4 was launched by that failed batch.

Read-only sorted-ID checks confirm all65536 unique IDs and correct i/i+32768
pairing. The snapshot itself is Peano-reordered, not ID-sorted. Thermal
asymmetry is localized near the interface:20/32768 pairs exceed1e-8,6 exceed
1e-7; its mass-mean relative error is5.714e-11. Velocity mass-RMS asymmetry is
5.441e-10. This is small compared with current native/reference differences;
it does not establish whether asymmetry growth is roundoff or a BC defect.

A separate `mirror_failure_t24_v1/t24` saved-state comparison visibly labels
GIZMO diagnostic-only/failed mirror gate. Actual-cell native density and
nearest GIZMO density samples use the same0.8..2.8 colors. Both show the RT
bubble/spike pattern, with visible interface differences. These two display
reconstructions do not define a density-error norm. At2.4, density ranges are
[.913591,2.120586] and[.914986,2.119435]; vy RMS .0996402 versus.0972865.
Same point-mass seed-kx2 moments differ7.23%(mass),8.13%(x momentum),5.19%
(y momentum). Point-potential total energy residuals are-5.01e-15 versus
7.37e-11. Neither method is an exact reference; the failed mirror certificate
is not waived by these small energy residuals. The figure was visually checked.

After investigating that terminal failure, a deliberately separate unchanged-
solver diagnostic continuation was launched in
`continuation_after_mirror_failure_v1` to4, from the actual terminal2.4 restart
pair. Source run.c writes the restart at final exit; runtime first sync is
19201 at2.400125. Binary, IC, EOS, gravity, timestep and original1e-8 mirror
threshold are unchanged. The continuation records both core invariants and
strict mirror result separately; it never represents completion as a certified
reflection reference. It tests how the measured asymmetry grows rather than
silently accepting it. New directory and lineage hashes preserve the failure.
CPU24,25; maximum2 snapshots(~25MB), restart allowance512MiB,149.10TB free.
No native RT or cylinder process was interrupted. No3D evolution was started.
