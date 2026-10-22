# H200-first single paired transition dispatch — 2026-10-22

Goal ACTIVE, unachieved: t6/8/10 and longest normalized GIZMO stability.
Prepared 2026-10-22T00:03:21.099420+00:00.

## Current live handles and accepted state

Original positive-rule producer launcher2166274 / ranks2166305–2166308
continues on Ada, unchanged and not restarted. Committed prefix at this report:
278136/1138432triangles (24.43%), all compressed;
498437632original ->60850479retained points (8.19x within this prefix).
Projected total from start remains approximately5.2–5.4h. Only FULL ready
completion can release the next physical step.

Current dependency is **material_full_positive_dispatch_v1**, PID2351718,
state=verified_live_original_producer_wait. Its actual ps handle was revalidated; do not infer
liveness from state.json alone. It polls the original producer60s apart.

The Ada-only wait process2296727 was replaced childless before any physical
run. The dispatcher's first waiting revision2347915 was corrected before ANY
submission to handle squeue terminal/purged-job replies through sacct of the
SAME job. Both obsolete handles are gone. Replacement records and previous
source/PID/state are preserved. No GPU preprocessor/trajectory was killed.

Current accepted trajectory is still original-rule step4 at
t=8.2355880701376864e-5. The prepared actual next high k1 remains reusable;
its force and mass solve are not replayed. No new accepted time is claimed.

## Identical physical sources; architecture-specific builds

material_full_positive_transition_v2_ada_mpi4 and
material_full_positive_transition_v2_h200_mpi4 have **68identical physical
source files**, including the original native C entry and full operator.
See dispatch/identical_physical_sources.json for hashes. Ada uses sm89,
H200 sm90; CUDA/C++ were rebuilt together. Both build statuses are
BUILD_PASS_NO_PHYSICAL_EXECUTION. No runtime H200 success is claimed yet.

The first high k1 acceleration is the retained actual Ada result, also when
H200 runs the new comparison. Missing high source summaries then run on the
selected device without force/CG. This first-stage comparison includes any
hardware numerical difference rather than claiming an all-H200 regenerated
high force. Future high stages2–4 are evaluated on that same chosen device.
Do not mix Ada timings into the H200 performance report.

## Resource observation and actual dispatch policy

Read-only Slurm observations: normal MaxTRESPU gpu8; h200 partition QoS
part_h200 has GrpTRESgpu6; syn104 has8H200,6allocated, CPUAlloc46/64.
The current running list shows5untyped GPU allocations in h200 and1typed H200
in the dedicated sjjlee_hgpu partition. Total node allocation is not identical
to a proven part_h200 QoS counter, so neither immediate availability nor
exhaustion is inferred from those physical counts alone. No other job/QoS was
changed. The actual scheduler will decide AFTER full preprocessing completion.

After completion and ready/hash checks, the new dispatcher submits exactly
one job.sbatch with:

    --partition=h200 --nodes=1 --ntasks=4 --cpus-per-task=4
    --gres=gpu:H200:1 --mem=128G --time=02:10:00

No GPU is reserved while the CPU-heavy preprocessing remains incomplete.
The Slurm CUDA_VISIBLE_DEVICES allocation is preserved. Another user's other
GPU on syn104 being busy is not a reason to reject our allocated device.

If the job is PENDING after600s, scancel --state=PENDING is requested. Ada is
used only after authoritative terminal cancellation AND absence of a physical
claim. A job which won the start race is observed; it is not killed as a
pending fallback. Unknown/transient status continues observation of the same
job rather than starting or resubmitting another calculation.

H200 and Ada both acquire a shared atomic O_EXCL
physical_transition_claim.json BEFORE entering the actual native runner.
The claim is never automatically removed. A started/failed/partially accepted
physical run is preserved without an automatic second execution. Ambiguous
submission or missing accounting also cannot authorize a duplicate. Ada uses
MPI4/OMP8 on its one card; H200 uses MPI4/OMP4 on its one allocated card.
The purpose is one actual same-stage RK4 transition; long production remains
subject to the already documented accuracy/cost/map guards.

## Verification performed this turn

- Native CUDA/C++ compile/link for H200: exit0. Same68source hashes as Ada.
- Shell syntax and Python AST parsed; source/binary/runner hashes checked.
- Already cancelled owned job433103 was queried READ-ONLY: squeue returns1
  with Invalid job id; sacct returns433103|CANCELLED by10396|0:0. The current
  exact accounting format/parser detects its terminal state. No new probe job,
  cancellation, force evaluation or accepted interval replay occurred.
- The original producer and current dependency were confirmed live by ps.

## Scientific diagnostics and strategic review

See RESULT_SHEAR_AND_COMPRESSION_DIAGNOSTICS_20261022.md. New saved-data-only
all71080cell analysis2.196s separates uniform bulk compression from the
positive-rule Jensen residual; at step4 the latter is only about1e-16 and
unresolved beyond numerical arithmetic. The exact material-map shear and
near-origin full-gradient dt reference were derived separately without
prescribing motion. These are not evidence that t>6 heating is solved.

One bundled follow-up was sent to the existing Astra planning agent for the
conditional continuous enriched-map guard, shear lifetime and minimum future
compression diagnostics. It must not alter/poll current runs or pre-judge the
unexecuted paired cost step. Expected memo:
ASTRA_MATERIAL_SHEAR_LIFETIME_AND_GUARD_20261022.md. No repeated routine reviews.

## Next authoritative evidence

Read dispatch/state.json and the actual live handles. On actual completion,
backend H200 maps to material_full_positive_transition_v2_h200_mpi4/run;
backend Ada maps to material_full_positive_transition_v2_ada_mpi4/run.
Inspect its contract.json, all four actual stage comparisons, endpoint,
physical finite-step gates, v3 restart roundtrip and selected_rule_cost.json.
Only then adopt its step5accepted_r[0-3].bin. If a wrapper fails after native
success, audit/adopt the saved successful state; do not replay it just to
repair a wrapper. Preserve the original full scope and previous Astra cost
abort decision. Global injectivity, enriched spectral CFL, long material mesh
lifetime, conservative repartition and paper normalization remain unfinished.

## Astra memo received and latest live check

ASTRA_MATERIAL_SHEAR_LIFETIME_AND_GUARD_20261022.md was received and read in
full. Its conditional guard/lifetime/output contract is adopted for the
post-paired-step decision. It changes no current run, rule, dt, state or
existing accuracy/cost gate, and instructs no extra guard simulation now.
No further routine Astra call was made.

At 2026-10-22T00:04:28.539399+00:00, original producer
launcher/four ranks and current dependency PID2351718 were verified live by ps.
Committed prefix 281592/1138432 (24.74%), 8.19x within prefix.
No new H200 job or physical transition claim exists: preprocessing is still
in progress and no GPU is being held for the downstream step.
