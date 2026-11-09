# Current common-method smooth-wave qualification — 2026-11-09

The previous goal turn made progress by building/starting a reflected-wall GIZMO RT reference from the current128 native IC. This turn prepares a missing smooth-wave convergence check for the current common2D volume/pressure-work binary. Existing1D smooth-wave evidence and old frozen-acoustic/entropy-star diagnostics do not establish this current2D whole-period result. No3D evolution is started.

## Physical and numerical contract

Periodic unit square, gamma=(double)(float)(5/3), rho0=1, P0=1/gamma, c0=1, zero external gravity. For integer wavevector k=(1,0) or(1,1), theta=2*pi*(k dot x - |k|*t), initial rho=1+epsilon*sin(theta), P=1/gamma+epsilon*sin(theta), v=epsilon*sin(theta)*k/|k|. Each case advances one period1/|k|. Six baseline cases use16²,32²,64² andepsilon1e−5. Two diagonal32² controls halve epsilon and the timestep cap separately.

Density, pressure and momentum are integrated analytically over each initial square cell. In particular, rho*v contains epsilon*sin(theta)+epsilon²*sin²(theta). Cell masses may differ initially and are then fixed individually. The native velocity is cell momentum/mass. U is the physical thermal integral, not a thermalized subcell kinetic variance. This omission from total represented kinetic energy is O(epsilon²); the reference is a linearized wave, not an exact finite-amplitude nonlinear solution. The half-amplitude control tests contamination at that order.

The existing native CSV converter/checkpoint codec is reused in a frozen local copy; only its auxiliary gamma metadata changes from double1.4 to nativefloat5/3. The regular production IC writer and hydro sources are unchanged. Every written x/y/v/U/W/m/V/ID field roundtrips bitwise. Independent12-point tensor quadrature agrees with analytic rectangle sine/sine-squared integrals within2e−14. Independent actual-polygon triangle8/16 integration of all eight pristine states passes: normalized density discrepancy≤1.74e−13, pressure discrepancy zero at stored precision, momentum discrepancy≤2.26e−16. Polygon volumes agree within1.12e−16 relative. The32² half-dt control has a byte-identical initial checkpoint to baseline.

## Existing runtime and no test-specific physics

Frozen executable SHA25611a0234acc882f6f2fc405b0053c2e5210262863be665f88cc5ae91b821ddf4e is identical to current Sedov/Gresho/RT. The existing periodic restart route is selected with EUNHA_IC=gresho; this label chooses the periodic driver, while the loaded checkpoint defines a sound wave rather than a vortex. No Gresho force law is present in that material path. Existing SSPRK2, pressure/density/velocity reconstruction, dynamic weight solve, endpoint-volume tolerance1e−11, C=.4, original acceptance/retry rules and disabled LU/memo options are retained. Baseline maximum dt=.1/N; one control halves it. No sound-wave limiter or damping coefficient is added.

The finite supervisor waits for the actual glass Sedov worker to finish and its native GPU handle to exit, then uses the same third Ada GPU slot on CPU12,13. It never kills/restarts existing cylinder/RT/Gresho/glass trajectories or allocates a fourth concurrent GPU evolution. Initial and endpoint checkpoints are retained, with one requested midpoint: at most3 per case including input,24 total,26,729,184bytes; free space at preparation149,066,654,679,040bytes. Effective params and environments are saved per case. A failed invariant stops the batch rather than changing numerical tolerances. The six-hour dependency wait checks the live supervisor.

## Assessment and interpretation

Independent actual weighted polygons are rebuilt from each saved final state using the existing geometry diagnostic, checking positive convex cells and volume agreement. Triangle8/16 quadrature integrates the linear reference on those cells; density, pressure and momentum integral L1 differences are normalized by epsilon times total volume. The analysis also records transverse velocity and the material-ID cell-mean entropy proxy. These are not new evolution variables or entropy repairs.

Only after actual completion will log2 error ratios across resolutions be reported. The resolution sweep refines dt with h, so it is coupled space/time convergence. Use half-dt and half-amplitude controls plus quadrature estimates before assigning order. Completion/conservation do not automatically establish second-order accuracy, nonlinear shock/shear morphology, full1D/2D/3D equivalence, MPI scalability or <=2x GIZMO timing.

Evidence root: Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/common_acoustic2d_20261109_v1. prepare.py, convert.c, compile receipts, plan.json, per-case ic_audit.json/initial_assessment, run.py and assess.py are retained. At this document creation the initial audits are complete; time evolution is queued, not completed.


## Startup guard error and isolated recovery

The first k10_n16 launch in `common_acoustic2d_20261109_v1` exited99 before any accepted step. Its requested one-period stop1.0 exceeded inherited `GFS_KDK_MAX_WINDOW=0.4`. This variable is the legacy driver's run-duration guard, not the timestep or integrator. Both failed supervisor525970 and child548387 were confirmed absent; no accepted trajectory was restarted.

The isolated `common_acoustic2d_20261109_retry1` preserves all original files and changes only this environment value to1.0, within the source-supported limit2. The common binary, ICs, CFL, timestep caps, quadrature, volume tolerance and all physical operators are unchanged. The original failed attempt remains intact. Recovery provenance is in `recovery.json`; the same24-checkpoint bound and26.73MB estimate apply, with149.05TB free before launch.

The16² axial case has now completed one full period:161 accepted steps,0 retries, fixed mass, positive U/V, energy drift3.35e−17, max volume mismatch1.95e−15. Actual-cell density L1 divided by the initial wave amplitude is.0972817, not a9.7% error of the background density. The32² case is running. Refinement, diagonal propagation, half-amplitude and half-timestep results remain pending; no convergence order is inferred from the16² case alone.


## Completed one-period direction and control assessment

All8 trajectories reached their original full-period endpoints with the unchanged common binary. Fixed masses, positive states, unit-box volume, geometric-volume tolerance and global energy checks passed. These are separate from the accuracy assessment below. The recovery changed only the run-duration guard, never the time step or hydro method.

| Case | Density L1 / epsilon | Momentum L1 / epsilon | Fitted density amplitude / initial | Phase error (radians) |
|---|---:|---:|---:|---:|
| k10_n16 | 0.0972817259 | 0.0969885826 | 0.848177639 | -0.0251742863 |
| k10_n32 | 0.0400325587 | 0.0405479659 | 0.969613255 | -0.0146544226 |
| k10_n64 | 0.0116253939 | 0.0116979793 | 0.993879101 | -0.00445052945 |
| k11_n16 | 0.156192493 | 0.157007361 | 0.741764159 | -0.0196557417 |
| k11_n32 | 0.0570260192 | 0.0581285284 | 0.935239129 | -0.0146903579 |
| k11_n64 | 0.0204963671 | 0.0205663348 | 0.987000886 | -0.00461147917 |
| k11_n32_half_amplitude | 0.0570267819 | 0.0581286245 | 0.935238738 | -0.014690706 |
| k11_n32_half_dt | 0.0569226308 | 0.057961214 | 0.935242465 | -0.0140949542 |

The density refinement orders are1.281 and1.784 along k=(1,0), versus1.454 and1.476 along k=(1,1). The latter does not demonstrate second-order accuracy on the tested16/32/64 sequence. At fixed N the diagonal wavelength is shorter by sqrt(2), so these errors are not a comparison at equal cells per wavelength. Both return after one full period; no incomplete propagation is used.

For diagonal32, halving dt changes normalized density L1 by−0.1813% and momentum L1 by−0.2878%; halving amplitude changes them by+.001338% and+.0001654%. Quadrature8→16 changes wave cell averages by approximately machine precision. Thus neither the tested time cap nor nonlinear-amplitude contamination explains the dominant error. This points toward spatial reconstruction/mesh coupling, but does not prove a particular limiter or geometry operation is responsible. The phase diagnostic still changes4.05% under half-dt, so time discretization is not absent.

The volume-weighted sine/cosine regression uses the exact sine and cosine averaged over the actual saved power cells, not generator-point values. At64, density amplitude ratios are.993879 axial and.987001 diagonal, phase errors−.00445053 and−.00461148rad, with residual RMS normalized by epsilon .0146096 and.0246540. Errors are measured relative to the initial perturbation amplitude, NOT the background density. Neither profile agreement nor conservation alone is a physical accuracy certificate.

Artifacts: `interpretation_v1/report.json`, `analyze_convergence.py`, and visually inspected `figures_v2/convergence.png/.pdf`. The original figure had overlapping automatic minor-tick labels; v2 changes only tick formatting and matches legend/body font size. The reference N^-2 line is a guide, not a fitted order.

Source inspection confirms `gfs_material_scalar_project2d.h`, `gfs_material_gradient_project.h` and `gfs_material_project_shared2d.h` match the frozen common source hashes. Scalar gradients enforce neighbor extrema and positive face values, and the same prepared states feed C/CUDA. This identifies a candidate path for a saved-state reconstruction audit; it is not evidence to disable limits or introduce a smooth-wave-specific correction. Full2D qualification and the ≤2x matched-accuracy GIZMO runtime target remain open.


## Native scalar-reconstruction audit of the saved initial waves

A separate read-only C shared library calls the frozen `gfs_project_one_scalar` and `mg_project` functions directly. It records the existing least-squares gradient before projection and the projected gradient, without changing either function. The diagnostic consumes the saved exact cell-average initial waves and their regular periodic face graph. It compares reconstructed face values with the analytical FACE-AVERAGED sine, not with a point value. No evolution, production binary, limiter or CUDA state is modified. Sources, compilation command, library and hashes are in `scalar_initial_audit/`.

| Case | Density face L1/epsilon before projection | After projection | Cells with changed gradient |
|---|---:|---:|---:|
| k10_n16 | 0.00439899784 | 0.0137089888 | 256/256 |
| k10_n32 | 0.00104274715 | 0.00325852851 | 1024/1024 |
| k10_n64 | 0.000256932272 | 0.000791462161 | 4096/4096 |
| k11_n16 | 0.00861111306 | 0.0236119225 | 192/256 |
| k11_n32 | 0.00207272217 | 0.00607701088 | 896/1024 |
| k11_n64 | 0.000513048958 | 0.00153028811 | 3840/4096 |

The scalar neighbor bounds act over much of the smooth initial wave, increasing face-value L1 by roughly3x. The strict positive-face margin is inactive for all these states; disabling positivity would not address this observed increment. Crucially, the projected INITIAL face-value errors still converge at about2nd order: axial2.073 then2.042, diagonal1.958 then1.990. Thus these results do not establish that initial scalar projection explains the1.48 evolved diagonal order, and do not justify turning the limiter off for a wave-specific passing result.

Diagnostic consistency checks give an exactly zero gradient for constant positive density/pressure. Joint rigid rotation of the graph and velocity state through32 angles rotates all measured gradients with maximum absolute difference4.07e−20. This is NOT the same as propagating a differently directed wave on a fixed grid. `consistency_checks.json` also records reconstruction/jump refinement ratios. Source headers match the frozen common copy; the diagnostic uses the same face-measure least-squares and pressure-offset definitions.

Next causal separation remains the velocity reconstruction, the assembled semidiscrete force/volume response, and evolved mesh geometry. A full native RHS audit must verify its exact compile-time path before comparing it with production; this scalar-only check is not relabelled as full hydro verification. No new dissipation coefficient, smoothness sensor or test-specific correction is introduced.


## Full free-RHS and evolved-geometry separation

The CPU reference now executes the frozen common trace build, face-normal velocity preparation, shared scalar projection and `gfs_material_rhs_maps2d`. Relevant physical macros were checked against `common_noh_volume_20261109_v1/native/extract_command.json`; the current portable build preserves those archive members. The CPU branch is explicitly present beside the CUDA branch in frozen `gfs_material_trace_probe.inc`. This audit does not execute the native volume inverse, time integration or GPU path, and is not a bitwise production-parity claim.

For exact initial regular cells, the integrated analytical pressure force and material volume rate are compared to one RHS at16/32/64/128. Axial force refinement orders are2.062,2.023,2.006; diagonal2.101,2.033,2.009. Volume-rate orders approach2 in both directions (last pair1.961 and1.962). The diagonal initial-state velocity limiter is inactive at all these resolutions. Global force, total volume rate and compatible total-energy rate cancel to roundoff. The128 case is a new instantaneous analytic state, not a completed128 trajectory. Sources, build commands, raw arrays and metrics are in `rhs_initial_audit/`.

Next the same analytical wave was integrated over the ACTUAL final saved power cells, replacing evolved thermodynamic/velocity values only in diagnostic memory. Power-plane support, reciprocal faces, geometric closure and areas were independently checked. Native resolved-face roundoff handling (256 eps times generator separation and the counterpart rule) was copied exactly after the first unfiltered diagnostic failed. No production geometry or tolerance changed.

For axial16/32/64, the manufactured-field RHS errors are essentially unchanged by the final geometry. For diagonal16/32, force L1 rises from.03555/.008288 to.04322/.01733, and volume-rate L1 from.08219/.02439 to.10828/.03756. The independently reconstructed64 geometry missed the native strict reciprocal-face-length tolerance on THREE pairs, with differences9.99e−16 to1.413e−15. This is a postprocessor representation failure; the production simulation completed normally. Failed reports and geometry are preserved.

For that64 diagnostic, both choices of one already reconstructed face as the shared lower/upper-owner representation were tested, with exact reciprocal coordinate transforms and NO tolerance changes. Maximum geometric change1.413e−15, volume error≤4.12e−14. Both pass the unchanged C operator and give force L1.004400626/.004400627, volume-rate L1.011907052 in both cases. Original regular64 values were.002024701 and.006487668. Sensitivity is far smaller than the increase; these diagnostic representations are not relabelled as native saved geometry. Results are in `rhs_evolved_geometry_shared/report.json`.

This implicates the reconstruction's response to mesh topology/geometry, without proving the complete origin of the time-integrated error. It motivated the controlled vanishing-face experiment documented in `COMMON_RECONSTRUCTION_CONTINUITY_2026-11-09.md`.
