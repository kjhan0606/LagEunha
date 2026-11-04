# Native Kepler exact20pi assessment, 2026-11-04

This package evaluates the existing Syntax job449508. It launches no solver,
changes no trajectory, and adds no physical operator. All work is under
/gpfs. The exact endpoint and final whole-disk assessment are now complete.
Job449508 is COMPLETED0:0 at20pi, with227343 accepted steps fromt0.
See final20pi/final_review.json and final20pi/goal_completion_audit.json:
the practical20pi stability objective is supported under the user's latest
GIZMO-relative criterion, with the scope limitations below retained.

The active user objective is Kepler disk stability through20pi<63, with
self-consistency, simplicity and improved evolution efficiency. The latest
user instruction accepts errors comparable to GIZMO at the corresponding
time. The historical cohort1% and custom wall budgets are retained as
diagnostics, not independent rejection criteria. No new t120/600 or vacuum-
boundary development is made a prerequisite for this bounded objective.

## Authoritative inputs and assessment

- Main run: /gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/native_u64_kdk_individual_20261104_v1/kdk_block2_wallkdk_v1/run_wall_lts_kdk_t4030_to20pi_v2
- Original reference: /gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/native_u64_w0_initial_20261101_v1/initial696.bin
- Acceptance context: /home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs/KEPLER_GIZMO_COMPARISON_CRITERIA_2026-11-04.md
- Evaluator: assess_endpoint.py; call only after the native completion marker.

The evaluator verifies the exact footer/receipt time, native state integrity,
original IDs and cell masses, zero W/pi, U64/V64 EOS, polygon area/containment,
and original-t0/continuation error metrics. It scans only accepted prefixes
of the actual RK4→globalKDK→block2 lineage. It reports global and original-
disk angular momentum, energy, density and phase/radial errors. The original
disk has49,185 IDs selected with inclusive.5<=R0<=2; the published annular
L1 uses current.8<R<1.2 with equal-element weighting. These are not mixed.

The density figure uses actual mass/current Voronoi volume, no smoothing,
centered[-2,2] axes, linear scale[.5,2], jet and fixed color ticks. Large
changes are also reported numerically: original-disk density ratios,
quantiles, and mass/count tails at10%,50%,100%. These are descriptive
distributions, not new failure thresholds.

## Evidence completed before the endpoint

- baseline_source_build_integrity.json binds22 unchanged source/build files
  to the baseline manifests and records the actual H200 GPU identity.
- interim_t61/lineage_ledgers.json covers222,158 accepted steps t0→61.00045
  with zero recorded invalid/clipping/floor/dual-energy injections.
- interim_t61/EOS_geometry.json: pressure-EOS relative max5.91e-8;
  total polygon area error-7.37e-15; all sites inside the polygon.
- interim_t61/density_t0_t61.png shows the preserved ring and edge noise.
- matched_time_throughput_55_to60.json: current production path is5.24x
  faster than the preceding RK4/circle path; GIZMO is3.34x faster than the
  current path in the same physical interval. Different boundaries and
  operators prevent an isolated integrator-speed attribution.

## One bundled Astra stage evaluation

The existing astra_next_kepler_stage agent reviewed the updated criteria,
current native code, saved states, figures and timing evidence. Its final
assessment was conditional on the still-pending exact endpoint:

1. Finishing the existing exact20pi endpoint audit is sufficient for this
   stage; new spin/weight operators or additional long targets are not
   required. Its earlier recommendation to reject solely from the old1%
   and wall budgets is withdrawn under the latest user instruction.
2. Check the sum of cell volumes and expose the density tail hidden by
   color saturation. Both checks are included in assess_endpoint.py.
3. The adiabatic-volume map is consistent with fixed-mass reversible PdV;
   it is not an arbitrary thermal-energy deletion. Total finite-step
   energy error must still be measured, rather than claimed exactly zero.
4. Block2 uses bin-dependent pressure-kick time weights. Instantaneous
   pair-force/torque cancellation does not imply exact finite-step
   momentum conservation. Report measured drift and unlogged pressure-wall
   torque; do not claim exact conservative multi-bin LTS.
5. The existing fixed-time GIZMO velocity graph is a valid magnitude
   comparison, while a density-error ratio at the same time is not known.
   A satisfactory whole-ring endpoint supports the bounded stability claim
   without inventing a missing density reference or rerunning GIZMO.

This is a stage criteria/source review, not a claim that Astra inspected a
future endpoint. Final numerical and visual inspection must use actual files.

## Scope and limitations to retain in the final result

The trajectory uses native RK4 through40.000091, globalKDK through40.2,
and then block2 KDK to20pi. It is not a pure-KDK trajectory fromt0. The
exact adiabatic-volume map starts at40.02 and preserves inherited entropy
errors. Configured Gamma=float32(1.666667)=1.6666669845581055.

The tested branch is fixed-mass W=0 Voronoi, no orbital steering, no
explicit viscosity (Re=infinity), and intrinsic cell spin inactive. The
CUDA build defines GFS_VARPRESS_W0: the active momentum force uses the
geometric variational pressure gradient and its corresponding pressure work
(exam_gpu.cu1246–1270,1290–1315). The generic params.dat av_mode=5 label
does not make that compiled pressure branch HLLC momentum flux. The thermal
map is the corresponding reversible fixed-mass adiabatic-volume law. This
qualifies the smooth cold Kepler configuration, not shock accuracy for every
other build/mode of the general hydro code.

The weight/inertia/thermal native layout and GPU geometry retain a Laguerre
extension path; evolving weights and spin are not claimed validated here.
Earlier spin/strain and dynamic-weight evidence is retained rather than
repeated solely to chase errors below the GIZMO target.

The current run uses one MPI rank, eight OpenMP threads and one H200.
Block2 offers periods1/2; geometry is evaluated globally at fine events.
Deep-bin activity estimates are hypothetical, and multi-rank scaling is
unverified. The5.24x improvement is measured for actual production paths,
not an assertion that complete individual-step or multi-GPU work is done.

The reflective polygon and polar Q1 particle layout differ from the paper
reference. The annular velocity comparison does not establish universal
method superiority. KDK pressure-wall torque and cumulative affected-cell
IDs are not fully logged; no fictitious zero replaces those missing data.
Core/background error is separately visible in the original-radius metrics.

Final evidence is in final20pi/summary.json, final20pi/lineage_ledgers.json,
final20pi/disk_morphology_metrics.json, final20pi/EOS_geometry.json and
the final density image. Root and Astra inspected the actual final image
and numerical distribution. final20pi/final_review.json supersedes the
pending manual-review placeholder preserved in the immutable summary.json.
The separate GIZMO449138 reference job remains running under prior
authorization and was not changed.
