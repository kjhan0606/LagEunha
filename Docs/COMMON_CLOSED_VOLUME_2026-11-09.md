# Closed-domain endpoint-volume qualification, 2026-11-09

The preceding goal turn made progress by integrating a common read-only stage-volume audit. This turn builds a private extension of the existing endpoint inverse to periodic domains and reflected stationary walls. The running high-resolution Cylinder executable, canonical sources, IC and pressure-work model are unchanged. No 3D evolution is started.

The dimension-independent algebra wrapper chooses the minimum immutable material ID in each positive-edge connected component as the constant-weight gauge. Only the auxiliary linear system replaces that row with a zero weight increment. It then checks every ORIGINAL raw-face equation, including each omitted physical row, at the existing 1e-10 relative residual. It never subtracts an RHS mean, modifies target volumes, changes masses, or relaxes nonlinear volume tolerance. Fourteen actual C/MKL tests pass with reuse off/on, covering disconnected compatible components, incompatible components whose global RHS sums cancel, zero RHS, negative coefficients, directed nonsymmetric compatible equations and isolated cells. This is algebra qualification, not fluid validation.

Private native adapter resolves existing periodic aliases and RT reflection ownership, computes the same A/(2d) coefficient and samples the previously integrated stage audit. Nonlinear tolerance remains 1e-11, maximum iterations eight. A nonperiodic RT coordinate uses its direct y difference; dividing by a zero periodic length is excluded.

## Native qualification: not accepted

The same existing checkpoints and four-step nominal intervals as the preceding audit were used. All executables reach their requested endpoint, but reaching it does not qualify this candidate:

| Test | Accepted steps | Retries | Linear calls / accepted | Wall seconds |
| --- | ---: | ---: | ---: | ---: |
| kh | 4 | 0 | 0 / 0 | 1.094 |
| rt | 88 | 387 | 0 / 0 | 66.184 |
| sedov | 79 | 377 | 380 / 3 | 51.682 |

KH already meets the nonlinear volume threshold, so its four successful steps do not exercise the new solver. RT and Sedov mostly shrink timesteps until no inverse is needed. This is a failed performance/implementation gate, not a successful common method or justification to relax tolerances. RT has no successful inverse call; Sedov has only three successful calls out of 380. Their final volume discrepancy near 1e-11 therefore cannot establish an effective inverse method.

A diagnostic executable identifies RT's first early return: a mirrored edge with a different owner has length 8.67e-19 at separation about 0.011. The existing pressure/RHS adapter excludes edges below its 256*DBL_EPSILON*distance coordinate-roundoff threshold before enforcing self-wall ownership; the new raw-CSR inverse checks ownership first. This contract mismatch must be fixed by sharing the existing resolved-face construction, not by inventing a test-specific cutoff. The diagnostic log is `rt_owner_diagnostic/log` and original states/results are retained.

Sedov reaches the closed linear wrapper, where the ORIGINAL residual test rejects most corrections, ranging up to 0.001394 relative error. The present records do not yet distinguish target-sum floating-point incompatibility from reciprocal face/Jacobian inconsistency. Next capture the failed original matrix, unmodified RHS, geometric target/actual volumes and correction, and localize residuals to the gauge row versus retained rows. Do not call this a roundoff-only problem without that evidence. Do not repair targets or claim acceptance from the auxiliary pinned matrix alone.

The common source remains at the previously qualified audit and anchored solver. The closed inverse is private and not adopted. Mathematical compatibility, the actual resolved-face contract and current full-row residual come before longer physics runs. General 2D, distributed MPI, 3D and matched-accuracy GIZMO runtime remain incomplete. Second-paper drafting remains after the first astronomy draft.

Evidence root: `Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/common_closed_volume_20261109_v1`. Includes private source/build commands, 14 algebra cases, exact input lineage, environments, output policy, actual run logs, failure diagnostic and `results.json` / `summary.json`.

## Shared resolved faces and finite-precision closed compatibility

The first failed Sedov linear system was captured before any correction was accepted (`sedov_first_system/closed_system.bin`, input-volume CSV and `analysis.json`). All 6,144 directed edges have reciprocal partners. The largest coefficient mismatch is 4.44e-16. The original equation residual is 1.789e-18 at the gauge row and at most 2.50e-25 elsewhere. The unmodified RHS sum is 1.789e-18; its infinity norm is 1.721e-9. Thus a relative-to-correction residual of 1.04e-9 fails the old 1e-10 gate despite an absolute volume discrepancy at roundoff scale. This statement concerns this captured system, not every historical failure.

The pressure adapter's existing resolved-face construction was extracted verbatim into `gfs_material_resolved_faces_native.inc`. Both pressure work and the closed volume Jacobian now use it. It retains the existing coordinate-roundoff cutoff, paired-face rule and reflection ownership checks; no new cutoff or boundary model was introduced. The cylinder retains its existing raw-face route. With inversion disabled, KH, RT and Sedov reproduce all 13 stored fields bitwise over the original four-step intervals. With the original algebra-only acceptance unchanged, RT now reaches the linear solver but all 387 attempts fail; Sedov remains at 3 successes in 380 calls. Sharing faces fixes an adapter mismatch, but does not itself solve closed-system compatibility.

The first RT system after this refactor was also captured (`rt_resolved_first_system`). Its original residual is 2.378e-18 at the gauge row and at most 4.50e-28 elsewhere; the RHS sum is -2.378e-18 and its infinity norm is 4.895e-12. Again all directed edges are paired and coefficient asymmetry is at machine precision. The apparently large relative residual, 4.86e-7, divides a tiny volume-sum discrepancy by a still smaller correction scale. Neither target volumes nor RHS were projected or adjusted in these measurements.

### Inexact Newton certificate in the existing geometric accuracy units

The private `physical_native` candidate keeps the auxiliary pinned-system algebra certificate at 1e-10. It then recomputes every ORIGINAL row residual in long double and requires

    max_i |(J delta_W)_i - (V_target_i - V_actual_i)| / V_target_i <= 1e-11.

This deliberately replaces the *closed full-row* residual divided by `max(abs(rhs))` when physical target volumes are available. It does not retain that old full-row criterion as an additional mandatory gate. The original relative algebra residual is still logged. The 1e-11 value is exactly the pre-existing nonlinear geometric-volume criterion, passed from the native adapter, not a new physical parameter or test-dependent coefficient. The updated weights must then regenerate the actual GPU geometry and pass that same nonlinear criterion before the timestep is accepted. The eight-iteration limit, target volumes, masses, pressure work, IC and timestep proposal are unchanged. The algebra-only entry point without target-volume scales retains its original strict certificate.

This is a change to numerical convergence certification, not a new dissipation model, a thermal correction, proof of exact material-face advection, or proof of improved shock morphology. An incompatible target larger than the requested geometric accuracy must still fail. Forty-two C/MKL cases pass, including the original 14, gauge changes, volume-unit scalings from 1e-12 to 1e12, rejected large incompatibilities, disconnected component defects whose global sums cancel, and invalid target volumes. Both direct and preconditioner modes are exercised.

In the first four-step native comparison, RT and Sedov each complete four steps at the original dt=0.000125 with zero retries and four accepted linear corrections. The maximum *actual rebuilt* relative volume errors are 7.41e-14 and 4.02e-12 respectively. Cell masses remain bitwise fixed and internal energies and volumes remain positive. KH needs no inverse in this interval and remains bitwise identical. All three inverse-disabled paths remain bitwise identical to the established control. These are coupling/solver tests, not completed KH/RT/Sedov physics validation or a GIZMO runtime comparison.

Evidence: `resolved_native`, `resolved_results`, `physical_native`, `physical_cases`, `physical_results` and the two captured systems, under the evidence root above. The original failed candidate, immutable fine-cylinder executable and existing results are preserved. An 80-nominal-step paired native comparison is the next qualification. Canonical integration is not implied by these private artifacts.

## Extended qualification and canonical opt-in integration

The 80-step matched checkpoint comparisons are complete. All cases use the same existing dt proposal, zero retries, bitwise fixed masses and positive internal energy/volume. Actual rebuilt geometry is checked after correction at each accepted endpoint.

| Test | Final time | Accepted linear solves | Maximum volume error | Inverse-disabled error | Wall-time ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| kh | 0.014 | 5 | 8.578e-12 | 8.259e-11 | 1.088 |
| rt | 3.011 | 81 | 2.605e-12 | 1.711e-05 | 1.263 |
| sedov | 0.11 | 83 | 9.217e-12 | 1.003e-05 | 1.245 |

These wall-time ratios compare the same native implementation with its endpoint inverse disabled/enabled on a shared GPU; they are not GIZMO comparisons, a scalability result, or a statistically isolated speed benchmark. KH now exercises five actual inverse corrections. RT and Sedov exercise 81 and 83, with no rejected linear corrections. No physical error or morphology claim is inferred solely from volume consistency.

Five scoped files have now been integrated into canonical source, retaining the existing opt-in `GFS_MATERIAL_VOLUME_COORDINATE=endpoint` activation. Default execution, fixed-solid cylinder behavior and the running high-resolution binary are not changed. Backups, exact patch and before/after hashes are in `canonical_before`, `canonical_closed.patch` and `canonical_integration.json`. The canonical build was made after synchronizing its source snapshot to canonical headers, excluding unrelated inherited radial-image experiments. Main C and native adapter objects were rebuilt; no new GPU kernel was introduced.

The rebuilt canonical executable reproduces the qualified private endpoint-enabled results and established inverse-disabled results bitwise in all 13 stored fields for KH, RT and Sedov over the four-step checks. A cylinder checkpoint with the original anchored endpoint inverse also remains bitwise identical in all 13 fields. General distributed MPI, 3D, long-time morphology and matched-accuracy GIZMO timing remain unqualified. Next assess predictor-stage as well as endpoint consistency, then longer 2D physical comparisons; do not substitute these small coupling checks for full KH/RT/Sedov accuracy tests.

## Predictor-stage necessity, temporal refinement and next full-history comparison

The previous goal turn integrated and built the common closed-domain inverse (progress). This turn tests its already available `stages` option rather than inventing a new pressure/heat model. KH, RT and Sedov each complete the same 80-step intervals without retries. Both the predictor and final rebuilt geometry satisfy the unchanged 1e-11 volume criterion. Actual linear solves are 160, 320 and 320 respectively. Shared-host wall times are about 1.47, 1.66 and 1.67 times the endpoint-only runs; these are local method-cost comparisons, not GIZMO performance qualification.

Eight additional native runs halve and quarter dt for endpoint/stages in RT and Sedov, retaining identical input states and final times (RT3.001->3.011, Sedov.1->.11). All finish with zero retries, fixed masses and positive U/V. Some runs have one negligible final exact-stop remainder step (161 or321 instead of160 or320); logs retain actual dt. In Sedov, the difference between endpoint and stages decreases by about four per halving: density relative L1 1.280e-6,3.180e-7,7.929e-8; pressure 2.893e-6,7.166e-7,1.785e-7. RT coarse steps exhibit a different local history, but the two finer levels' cross-method differences decrease at about second order. Do not describe this as universal second-order time convergence: within each method, observed Sedov temporal self-convergence is about1.05--1.09 in L1, while RT varies by field and mode (about0.73--1.66). Nonsmooth reconstruction/topology effects have not been isolated here.

Endpoint-only remains the economical candidate for the next whole-history comparison. This is a scoped choice: extra predictor inversions did not demonstrably remove the remaining field errors in these intervals, and cost more. It does not prove they are unnecessary in every flow, or that the whole physical method is qualified. Both options and all negative evidence remain available; no default or high-resolution Cylinder trajectory was changed.

### Signed excess work versus physical entropy

Saved first-RHS data already show negative `H=Udot+P*Vdot` at the selected Sedov/RT cells before time integration: SedovID330 H/U=-.676592 per unit time, RTID4727 H/U=-.00452201. Their cell-average entropy-proxy decreases persist under dt refinement: Sedov minimum delta log(P/rho^gamma) remains near-.00479 (about-.48%), RT near-4.98e-5. Neither endpoint nor predictor volume correction removes this behavior. The mean mass-weighted proxy increment is positive in these intervals.

This does NOT establish a negative physical viscosity or require an entropy repair. The current direct-flux header explicitly defines H as signed excess work, including reversible reconstructed subcell work; it is not the earlier adjoint model's nonnegative production term. A logarithm of cell-averaged quantities is not the cell integral of material entropy. The existing November7 history already discusses this distinction and rejected or limited entropy-star attempts (`MATERIAL_RT_PORT_AUDIT_2026-11-06.md`, ordered-direct thermodynamics and later entropy-failure sections). Do not revive those patches from this scalar diagnostic. Reference/spatial accuracy must decide whether the remaining proxy behavior is harmful. The initial v1 analysis used nominal parameter gamma; `stage_assessment_v2` uses the actual gamma recorded by the executable and is authoritative. Both versions are preserved.

Evidence: `stages_extended/{results,summary}.json`, `time_refinement/{results,analysis}.json`, `stage_assessment_v2/summary.json`, `stage_assessment/first_rhs_heat.json`. No new physics or acceptance threshold was introduced this turn.

### Full KH128 history now executing

`full_kh128` uses the frozen canonical executable c3554c696ec0b2044ac34ba20f54c4dfe7c920bf142eb6fd3e67a91f799a9ffb, 128x256=32768 cells, from the existing pristine dynamic-weight IC to t1.35 through checkpointed .3,.6,.9,1.2,1.35 segments. This is a new-method whole-history comparison, not a restart or replacement of an existing trajectory. Input x/v/m/specificU match the existing GIZMO128 IC exactly. GIZMO's initialized kernel density differs as previously documented; equal physical accuracy is not assumed. Only endpoint volume inversion is enabled, with previous-LU preconditioning off for this reference trajectory. All other retained IC, reconstruction, EOS and CFL settings come from the prior high-resolution KH configuration. The effective params/environment and exact source/checkpoint hashes are in each segment directory.

The finite supervisor advances only from each verified successful saved checkpoint and stops on a failed invariant or endpoint. Five final-output checkpoints are planned, approximately228MB including copied input checkpoints, plus diagnostic logs/field extracts. CPU12,13 and the shared Ada are used; this is not an exclusive GIZMO timing experiment. The high-resolution fixed-wall Cylinder binary/process is unchanged. Do not count a launched calculation as a physics pass. Use actual `status.json`, process handles and per-segment results when reporting.

## Closed-domain reuse qualification

The existing previous-LU right preconditioner was exercised in twelve ABBA native runs (off/on/on/off for each of KH, RT and Sedov) over the same 80-step intervals. Symbolic reuse is ON in both sides of this comparison, so the comparison isolates previous-factor preconditioning. CPU14/15 and the shared Ada were used; the long Cylinder and KH128 trajectories were not modified. The current pinned matrix supplies all GMRES products and its strict1e-10 residual certificate; the closed wrapper separately verifies EVERY original equation at the unchanged geometric-volume accuracy. No physical model, stopping tolerance or kernel changed.

Each repeated off pair and on pair is bitwise reproducible, and on/off accepted times and dt agree exactly. Masses remain exact, U/V positive, final geometry within1e-11. KH rho/U/V/velocity are bitwise identical on/off. RT maximum field differences normalized by the field's maximum magnitude are below1.77e-12; Sedov below5.10e-15. Preconditioning succeeds in4/4 KH attempts,78/80 RT and79/82 Sedov; unsuccessful attempts fall back to current factorization, rather than changing the equation or accepting a loose solution.

Mean ABBA linear-work ratios on/off are .870,.473,.842; whole-run ratios are1.024,.917,1.021. The smaller cases show no total-time benefit, so the option is not enabled globally. RT's measured whole-run reduction is about8%, not53% (53% refers only to inverse linear work). These shared-hardware measurements do not certify the GIZMO<=2x target. The continuing full KH128 reference deliberately retains its recorded immutable off configuration. Current KH128 profiling shows symbolic analysis is substantial; a representative high-resolution saved-state comparison is needed before choosing settings for subsequent production trajectories.

Artifacts: `closed_reuse/{results,summary}.json` and per-run IC/binary/environment hashes. The short temporal tests are not presented as a full KH/RT physics pass. A separate planar Sod reduction comparison, using the existing64x16 IC to its originalt=.1 endpoint, is now running on the same canonical method to check transverse invariance and exact cell-average accuracy before broader dimensional claims.


## Completed planar Sod reduction and matched comparison mask

The canonical 64x16 planar Sod comparison has completed t=0 to0.1 in500 accepted steps, zero retries, with endpoint inversion off/on. All masses remain fixed and U/V positive. Maximum accepted geometric-volume discrepancy falls from2.4531e-7 to1.2327e-14. The endpoint case makes500 successful linear corrections. Maximum transverse speed is1.56e-12 (off) and1.55e-12 (endpoint); transverse x-row spread is below3.14e-13. Reconstructed planar cell areas agree with the saved volumes within2.31e-12. These establish planar symmetry and the geometry coupling, not general2D shock accuracy or identity with the discrete1D solver.

The exact Sod solution uses the executable gamma1.399999976158142, integrated over each actual power-cell interval with64/128-point quadrature. The original whole-cell mask inside x=[.25,.75] selected31 cells without inversion and32 with it: a boundary moved across.25 by a tiny amount. Therefore the original per-case error difference must not be reported as an improvement. The authoritative cross-method comparison uses the same31 material IDs, covering length0.484375241878 and0.484375241881 respectively. At quadrature order128, relative density L1 errors are0.0301483385540 and0.0301484236956; pressure errors are0.0356845194210 and0.0356846026762. The approximately8e-8 differences are much smaller than the64/128 quadrature changes (about2e-5 in density and6e-6 in pressure). Accuracy is unchanged at this reference precision while geometric-volume consistency improves.

Evidence: `planar_sod/results.json`, `analysis.json`, and especially `matched_cohort.json`; scripts and profiles are retained beside them. The1D Riemann reference is an accuracy reference, not proof of algebraic dimensional equivalence. No physical patch, limiter, initial condition or tolerance was changed for this test.


## Evolved KH128 reuse and adaptive-step qualification

Six reverse-order native comparisons (off, symbolic, precondition, precondition, symbolic, off) advanced the same existing32768-cell KH checkpoint fromt1.2 to1.202. This saved state comes from the earlier full trajectory, not the currently running new-method reference. Canonical binary c3554c696ec0b2044ac34ba20f54c4dfe7c920bf142eb6fd3e67a91f799a9ffb and all physics/CFL/tolerances are identical. All runs reach the target in24 accepted steps with24 accepted inverse corrections, fixed masses, positive U/V and volume error below1.91e-12. The preconditioned pair is bitwise reproducible; off and symbolic repeats differ at floating-point level, quantified in assessment.json. No claim of bitwise reproducibility is made for those pairs. Preconditioning changes rho/U/V/velocity by at most1.007e-13 relative to each field maximum. Accepted times agree exactly; dt differences up to8.958e-18 remain and are recorded rather than called bitwise identical.

All six runs have10 existing directed stage2 CFL retries. The initial harness had incorrectly required zero retries and stopped AFTER the first completed baseline. Inspection showed all10 are ordinary stage2 CFL-bound checks, with no failed inverse. The completed baseline was retained and the five remaining comparisons were continued; no timestep/solver tolerance or acceptance rule was changed. Retry count and cause are now compared explicitly. Original harness, stop diagnosis and continuation are preserved.

Every inverse call changes the sparse pattern, so symbolic-only reuse yields no cached pattern hits. Mean whole-process times are14.9465s (off),15.0992s (symbolic) and13.6110s (previous-LU right preconditioner). Preconditioning succeeds in23/23 attempts in EACH repeat, reusing the first factorization; the current matrix still supplies all products and strict1e-10 raw-residual certification, followed by the full original closed-volume and actual rebuilt-geometry checks. Its linear-work ratio is0.4046 and whole-process ratio0.9106. These are shared-hardware short comparisons; the linear-work reduction is not the total speedup, and none proves GIZMO<=2x at matched accuracy.

This supports using the existing opt-in preconditioner for subsequent comparable evolved-KH computations. It does not justify a universal default (the earlier small KH/Sedov runs did not improve total time). The already running full KH reference and fixed-wall cylinder remain immutable. No new physics, geometric clipping tolerance, pressure limiter or thermal repair was introduced. Evidence: `kh128_saved_reuse/plan.json`, `results.json`, `summary.json`, `assessment.json`, `harness_stop.json`, plus each exact log and final state.


## First KH128 whole-history checkpoint

The new endpoint-volume trajectory reachedt=.3 from its pristine IC in3137 accepted steps and137 retries. Fixed masses, positive U/V and maximum geometric-volume error9.7224e-12 pass the coupling gates. The supervisor continued from that verified checkpoint towardt=.6 without replacing the binary or configuration. Compared with the earlier128x256 trajectory at the samet=.3 and identical material IDs, relative density L1 difference is3.66284e-5 and mass-weighted transverse velocity RMS differs by-0.0159572%. Total energy agrees with the initial value at the reported double-precision sum. These are early whole-history method differences, not reference error, full morphology validation or proof of improvement. Different timestep histories are retained. Evidence: `full_kh128/native_t03/result.json` and `full_kh128/assessment_t03/comparison.json`. Fullt1.35 and matched-accuracy GIZMO assessment remain outstanding.


## Noh exterior follow-up

The previously recorded Noh exclusion is now superseded by the opt-in common anchored volume adapter, qualified with unchanged existing exterior physics. See `COMMON_NOH_VOLUME_2026-11-09.md` and its source/result manifest. The original closed periodic/reflected and fixed-solid cylinder regression states remain bitwise unchanged. Full Noh64 evolution is launched, not yet accepted as a morphology pass.


## KH128 checkpoints t=.6 and .9: whole-history method differences

The same initial13-field arrays are bitwise identical. Old and new trajectories retain their different accepted timestep histories. Both completed endpoints preserve exact individual masses and total volume2, with total energy at rounding level. Saved periodic power-cell owners were independently certified at256x512 and512x1024 raster points; figures use cell averages without smoothing. These are method differences, not exact-reference errors.

At t=0.6, same-ID mass-weighted relative density L1 difference is 3.4832329e-05; common-space density L1 difference is 6.6961551e-05 on512x1024. Transverse RMS changes by 0.00019494223 relative, while the common-space transverse-velocity relative L2 difference is 0.020590112. The latter includes shifted cell boundaries and is normalized to a small perturbation. Do not replace it by the much smaller global RMS difference.

At t=0.9, same-ID mass-weighted relative density L1 difference is 3.2631625e-05; common-space density L1 difference is 9.8346767e-05 on512x1024. Transverse RMS changes by -3.3767834e-05 relative, while the common-space transverse-velocity relative L2 difference is 0.011019572. The latter includes shifted cell boundaries and is normalized to a small perturbation. Do not replace it by the much smaller global RMS difference.

The density morphology is nearly unchanged at these early/intermediate times; this does not establish final nonlinear accuracy. The new trajectory completed.9 and continues to1.2/1.35. A finite CPU-only analysis_pipeline waits for actual successful endpoints, checks supervisor liveness, and renders the same previous-method comparisons. The existing GIZMO128 final snapshot is verified to have time1.35 and32,768 cells, but no new-method final comparison is claimed before that endpoint exists.

## Final KH/GIZMO comparison preflight

The finite CPU postprocessor now waits for the current common-method KH128 t1.35 endpoint. It reuses the archived spatial Fourier-moment definitions and synchronized GIZMO final restart state, rather than comparing a native fully integrated velocity with a GIZMO drift/first-kick snapshot velocity. The old128 comparison is reproduced within1e-12 as a preflight. Exact input x/v/m/u equality is rechecked; GIZMO's own density initialization produces a mass-weighted pressure difference .0005813957953 from native P=10. This caveat is retained instead of relabeling the initialized states identical.

At completion the diagnostic will compare previous/common/GIZMO low mass/momentum modes and render density. Native panels use certified actual power-cell ownership; GIZMO uses nearest stored density samples explicitly labelled as such, not an asserted MFM Voronoi mesh. No density-error norm is derived from unlike display reconstructions. These are method differences, not exact-solution errors; shared-hardware runtimes are not certified as matched-accuracy speed. Output and live status are in full_kh128/gizmo_analysis_pipeline and gizmo_comparison_t135. No new fluid run was launched by this processor.

The current full KH128 history has now reached1.2 and continued to1.35 (nativePID239693 observed live). The .9→1.2 segment took3191 accepted steps and maximum volume mismatch4.264e-12. Automated same-space previous/current density comparison gives relative L1 .000231023 at256×512 and .000226853 at512×1024, approximately .0227%. The transverse-velocity spatial L2 difference is .0221153, approximately2.21%; it must not be replaced by the much smaller difference in a global RMS statistic. The inspected side-by-side density fields retain the same large-scale interface waves; developed vortex roll-up is not established by this image. This is close agreement with the previous method, not exact-reference accuracy or the pending final GIZMO comparison.

## Full KH128 completion and GIZMO comparison — 2026-11-09

`full_kh128` has completed all five segments to exactly t=1.35: 14,409 accepted steps, 1,649 retries. Every segment preserves each mass bitwise, positive U/V, and the unchanged 1e-11 actual rebuilt-volume criterion (worst 9.7224e-12 over the history). Final total energy relative drift is -4.53e-16 in the GIZMO-comparison accumulation. The final previous-method comparison and GIZMO worker both completed successfully; their figures were visually inspected.

Final point-mass low-mode relative L2 differences (mass, x-momentum, y-momentum) between common native and synchronized GIZMO128 are 1.4097%, 1.0308%, 15.7705%. The previous native/GIZMO128 y-momentum difference was 15.7866%. Between new/previous native, corresponding differences are .03627%, .01614%, .27316%. These modal norms must not replace same-space local-field norms or be called exact-solution errors. The interface remains a wave at this time, not fully developed vortex roll-up. The new geometry consistency alone has not removed the remaining KH/GIZMO growth difference.

Evidence: `full_kh128/gizmo_comparison_t135/comparison.json` and `density.png`, plus `assessment_t135/comparison.json`. Both native panels use certified actual power-cell ownership. GIZMO's panel is explicitly labelled nearest stored density samples, not actual MFM control volumes; no density error norm is computed between those unlike renderers. All input x/v/m/specificU agree, but the previously documented initialized pressure difference (~.0581%) remains. This is not an accuracy pass or matched-accuracy speed certification.

Full-history step/cost evidence and the new inverse/C=.8 qualification are documented in `COMMON_TIMESTEP_COST_2026-11-09.md`. No old successful Kepler or earlier method trajectory was replaced.
