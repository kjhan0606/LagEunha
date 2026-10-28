# Dedicated Ada execution of the native seven-case comparison

The previously authorized LagEunha RTX 5000 Ada server was idle when checked:
no GPU compute processes, 0% utilization, 2 MiB used, and ample host memory.
Its shared GPFS paths, CUDA 13 runtime, driver 580.82.07 and Intel MPI/Fortran
dependencies were verified. This dedicated server has no working Slurm
configuration. The comparison therefore uses the existing local-MPI execution
method previously used for the Ada tessellation tests, with two OpenMP threads
and one numerical process at a time. No simulation runs on grammar.

The new output directory is
`/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/vertex_halo_pair_t2654918_ada_20261029`.
Existing parent checkpoints and binaries are linked, not overwritten.
Estimated new output is about 1.8 GB: seven endpoints and twenty full native
face/force dumps. No new Slurm jobs were submitted, and existing Slurm jobs,
the main halo40 trajectory and its t=3 post-checkpoint hook were not modified.

## Numerical comparison contract

The same accepted t=2.6549183918680868 checkpoint and exact frozen binaries
are used as in the pending A100 comparison. Cases remain qualifier_old,
qualifier_new, control40, vertex40, vertex20, vertex20_reuse and
vertex20_reuse_local_sync. Physical parameters, target times, flags, complete
face/force dumps, native 15-check gates, exact ten-field endpoint comparison,
step traces and four RK4 geometry/operator comparisons are unchanged.

Only the local MPI bootstrap and execution-context validation differ. The
adapter records the actual hostname, GPU UUID, benchmark UUID and supervisor
PID. It does not invent a Slurm allocation number. It checks all requested
flags, the same actual host/GPU/session, and every recorded supervisor. All
cases use GPU UUID `GPU-42b3367f-1df7-2b2d-7961-7781c87b022d` on LagEunha.
Five synthetic context tests reject wrong flags, missing records, mixed
hardware/session, unrecorded supervisors and fabricated Slurm IDs, while
accepting an honestly recorded continuation. They are not numerical tests.

## Packaging correction without replay

The first two native qualifier calculations both exited 0 and passed their
15-check endpoint audits. The comparison then stopped because the stager had
omitted a root-level transitive Python import,
`audit_kepler_polar_long_gate_20261027.py`. The original stdout, stderr,
exit status, adapter, manifest and staging receipt are retained in
`packaging_recovery_20261029/`; numerical outputs were not replayed.

The missing reader was added and an import preflight added to the stager.
The qualifier comparison was rerun on its existing snapshots and passed:
all ten physical fields, printed step trace, tessellation/face counts and
native conservation/injection gates agree exactly. The original supervisor
was confirmed terminated; the continuation records a new actual PID while
retaining the same benchmark UUID and actual GPU. Only the five previously
unstarted numerical cases continue. Physics and geometry thresholds were not
relaxed. Initial status 1 remains recorded; the continuation has a distinct
exit-status file and logs.

## Observing this finite run

The continuation is executed through a live SSH command under a 2400-second
timeout. It is a finite comparison batch, not a login-node polling daemon.
The case wrappers refuse to overwrite existing numerical logs or endpoints.
Actual numerical evidence is in per-case logs and endpoint gates,
`qualifier.json`, and eventually `comparison.json`. Execution metadata is in
`manual_execution.json` and `manual_continuation.json`; overall continuation
status is `manual_resume_exit_status.txt`. An initial script exit or a live
SSH handle alone does not establish successful numerical comparison.

At this record the qualifier has passed and control40 is executing. Halo
certificate/reuse/local-sync correctness, net runtime gain, production
promotion, t>6 and the GIZMO long endpoint remain unproven. Ada results do not
automatically prove A100 throughput or MPI rank scaling.

Artifacts: [stager](stage_kepler_vertex_halo_ada_20261029.py),
[manual context adapter](audit_kepler_vertex_halo_manual_20261029.py),
[context tests](test_kepler_vertex_manual_context_20261029.py),
[stage receipt](KEPLER_VERTEX_HALO_ADA_STAGE_2026-10-29.json),
[packaging recovery](KEPLER_VERTEX_HALO_ADA_RECOVERY_2026-10-29.json),
[recovery implementation](recover_kepler_vertex_halo_ada_packaging_20261029.py).
