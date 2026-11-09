# Whole-history KH timestep and cost audit — 2026-11-09

The new endpoint-volume KH128 history completed t=0..1.35. Its 14,409 accepted steps and 1,649 retries are compared with 2,048 GIZMO synchronization intervals. All 32,768 particles are active at every one of GIZMO's 2,049 logged evaluations (including initialization). Individual stepping did not reduce the active population in this particular reference calculation.

Evidence and reproducible read-only analysis: `Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/common_timestep_cost_20261109_v1/{audit.py,report.json,accepted_steps.npz,profile_groups.json}`. Actual completed logs, effective environments, frozen source and binary hashes are recorded. No trajectory or timestep was changed for the audit.

| Quantity | Common endpoint inverse, KH128 | Archived GIZMO MFM128 |
| --- | ---: | ---: |
| Cells | 32,768 | 32,768 |
| Final time | 1.35 | 1.35 |
| Evolution steps / sync intervals | 14,409 | 2,048 |
| Median accepted dt | 9.56308e-5 | 6.59180e-4 |
| Process wall seconds | 6,620.783 | 1,250.637 |
| Process seconds / accepted interval | 0.45949 | 0.61066 |

The raw wall ratio 5.29393 decomposes into step-count ratio 7.03564 times mean interval-cost ratio .752444. The native GPU is shared on LagEunha; the reference GIZMO CPU run was on syntax. These are descriptive costs, NOT a same-hardware benchmark or matched-accuracy certification. An accepted SSPRK2 interval and a GIZMO synchronization interval are not identical computational work. The GIZMO cumulative CPU-log clock is 1,254.44s, slightly different from the archived process wall clock; the table consistently uses process wall clocks.

Disjoint native timing groups, excluding each segment's first accepted interval, assign 44.29% to the hydro adapter/face calculations, 27.04% to the endpoint inverse linear solve wrapper, 16.87% to GPU tessellation, and 11.80% to remaining geometry/transfer/overheads. Nested factorization timers are excluded from the sum. This does not support attributing the whole slowdown to Voronoi construction.

The actual material CFL uses

`dt = C min_i [V_i / sum_f(A_f max(c_i + |v_i.n-w_f|, c_j + |v_j.n-w_f|))]`.

The boundary normal speed includes generator velocities, face-centroid offset and `(Wdot_i-Wdot_j)/(2d_ij)`. It is a relative-speed expression, not laboratory advection. On a uniform, stationary Cartesian grid its limit is `C h/(2 d c)` in dimension d. For this 2D KH initial state C=.2 gives a quiescent-limit estimate 9.5683193e-5, close to the actual first bound 9.5682750e-5; the small difference includes the actual perturbed initial velocities/state. The geometric bound is at least 14.52 times the selected bound throughout the completed history; damping CFL is disabled. The acoustic face bound, not geometry or the .001 maximum step, therefore supplies the selected local restriction. Endpoint truncation and retry/predictor history can further shorten accepted steps.

The archived GIZMO MFM source uses effective default CourantFac=.4 and particle size divided by half its signal speed. The numeric Courant coefficients are not directly comparable. The current standalone 1D material driver uses `.3*width/max_face_relative_signal`, which also differs from the multidimensional face-sum normalization. This arithmetic is not a 3D runtime test or proof of an admissible larger step.

## Existing history and next decision

These CFL formulas and enlarged-step studies were already documented in `MATERIAL_RT_PORT_AUDIT_2026-11-06.md` (Geometric CFL normalization audit; same-start late KH controls; current-source full KH timing). They are NOT new discoveries. Prior C=.4/.8 KH, Sedov and shear controls and the full KH64 C=.8/GIZMO raw ratio 1.43548 remain valid within their documented scope; matched physical accuracy was not established. Do not repeat those old studies unchanged.

The new evidence here is the complete higher-resolution history with the new endpoint-volume inverse and its all-step cost partition. This coupling must be qualified at the existing efficient C=.8 setting before adoption. A separate `common_closed_volume_20261109_v1/full_kh64_cfl08` run now uses the same 8192-cell IC, C=.8, dt cap .001 and physical operator as the archived efficient candidate, with the new common endpoint inverse. It freezes the same canonical binary used by the just-completed KH128 reference. LU/memo acceleration remains off so the new comparison is not mixed with a second performance change. All original volume, positive-state, mass and energy checks remain enforced. The old and new CFL normalizations, coefficients and existing live reference jobs are unchanged.

The finite supervisor records effective params, lineage and process handles per segment, then a separate CPU-only worker compares saved t=1.35 states with the old C=.8 trajectory, current KH128 and synchronized GIZMO endpoints. This is a new-coupling accuracy gate, not a claim that increasing CFL is universally safe. Only 2D evolves; 3D and the full goal remain incomplete.


## New inverse with existing efficient C=.8: completed full history

The KH64x128 candidate completed all five targets through1.35, 1795 accepted steps and451 retries. Every mass stays fixed; U/V stay positive; maximum target/geometric volume mismatch9.9448e-12, maximum absolute relative total-energy drift1.814e-15. Total observed process wall404.833s on the shared Ada is not comparable as a controlled speed benchmark to the older154.551s native or107.665s GIZMO runs on a different load. No <=2x claim follows.

Final certified power-cell spatial comparison to the previous C=.8 full trajectory gives density relativeL1 .000811430 (about.0811%) and transverse-velocity relativeL2 .049642 (about4.96%) at512x1024 sampling. The coarser raster gives.000791704/.049659, showing the latter norm is not a display-only effect. Low y-momentum mode difference is.7151%; vyRMS differs-.2393%. Do not replace the larger local-field difference with the small RMS statistic. The visually inspected density plots retain the same large-scale interface shape.

Compared with the current128 C=.2 result, coarse C=.8 low y-momentum modes differ17.426%. Because both resolution and timestep differ, this is not a pure spatial convergence study. GIZMO64 growth disagreement remains about60.13% in vyRMS; the new inverse alone has not resolved it. C=.8 remains a candidate pending current-operator temporal controls and shock gates; no global setting changed. Evidence: `common_closed_volume_20261109_v1/full_kh64_cfl08/assessment_t135` and `assessment_summary.json`. Both finite hydro and analysis workers exited successfully. Original live cylinder/RT histories remain running.


Current-operator same-start temporal controls, polygon-intersection norms, and full Sedov timestep comparisons are now complete; see `COMMON_TEMPORAL_ACCURACY_2026-11-09.md`. These retain candidate status and do not change any production CFL.


## KH128 C=.8 complete cost audit

The same current endpoint-volume implementation completed t=1.35 in 3609 accepted steps and 815 retries, compared with14409 accepted steps and1649 retries at C=.2. Candidate process wall is 2017.170s versus6620.783s for C=.2 and 1250.637s in the archived GIZMO history. Observed candidate/GIZMO ratio is 1.61291; native accepted steps/GIZMO synchronization intervals is 1.76221.

These are shared-Ada versus archived-syntax-CPU measurements. The ratio below2 is NOT a same-hardware, same-load or matched-accuracy speed certificate. It must not be promoted to completion of the performance goal. There are no inactive-particle savings in the2048 GIZMO intervals.

Candidate disjoint profile fractions are40.78% hydro adapter,31.05% endpoint inverse linear wrapper,16.42% tessellation and11.75% remainder. The geometric bound exceeds the selected face bound by at least3.795 throughout this history. All scripts assert endpoint completion, binary identity, step/time sums, and profile accounting before writing the report.

Final spatial temporal differences and GIZMO mode comparison are recorded in `COMMON_TEMPORAL_ACCURACY_2026-11-09.md`. The final density L1 difference from C=.2 is0.04268%, transverse velocity L2 is2.9011%, and low y-momentum mode difference is0.31942%. The GIZMO discrepancy remains about15.8% in those modes. Efficient candidate status is supported; a global default and general accuracy/performance pass are not.

Evidence: `common_closed_volume_20261109_v1/full_kh128_cfl08/cost_audit/{audit.py,report.json,accepted_steps.npz,profile_groups.json}`. No live calculation was restarted or modified for this audit.


## Current-method exclusive H200 benchmark submitted

Job456949 (`lag-common-kh128-bench`) depends afterok on current portability gate456829. Existing old64-cell benchmark454560 and exact-reuse gate456838 are unchanged. This new job is necessary because the older benchmark lacks the current endpoint-volume coupling and uses a different resolution; its timing cannot certify the current method.

The benchmark freezes current portable binary11a0234…ddf4e, GIZMO binary, exact existing128×256 inputs, settings and source hashes. Native/GIZMO initial IDs, positions, velocities, masses and specific internal energies are bitwise identical. ABBA order is current native, GIZMO, GIZMO, current native, all through t=1.35, on the same exclusive H200 node and the same two allocated CPU cores. Native uses one H200 in addition to those CPU cores; GIZMO uses CPU only. No claim of equal hardware acceleration is made. Native retains its qualified C=.8, cap.001, endpoint tolerance1e−11, physical operator and five-segment history. GIZMO parameters match its archived128 reference; no viscosity or pressure adjustment is made.

Current H200 six-case validation must pass before this job can start; runtime rechecks that report and executable hash. A GIZMO t=1e−6 executable/IC/output preflight precedes costly repetitions and is excluded from timing. Native elapsed time sums solver processes including initialization/I/O in allfive segments; GIZMO times one solver process per history. Postprocessing/copies are excluded and total batch wall is separately recorded. ABBA spread and worst observed native/GIZMO ratio are reported; accuracy remains a separate requirement.

Output budget: at most30 native checkpoints including copied inputs,9 GIZMO snapshots including preflight, about711.4MB estimated checkpoint/snapshot storage and a4GiB total allowance. Frozen executables total18.0MB; preparation free space149.3TB. Native dump interval10 exceeds the end time within each segment; GIZMO stores initial and final endpoints. AST and Bash checks passed; the GIZMO endpoint reader was tested on the existing actual128 snapshot, including exact mass/ID and positivity checks. No GPU calculation was launched on the login node.

Evidence: `common_current_h200_benchmark_20261109_v1/{plan.json,prepare.py,run.py,run_native.py,benchmark.sbatch,preflight.json,slurm_submission.json}`. No measured speed result exists yet. The full dimensional/MPI/physical-accuracy goal remains active.
