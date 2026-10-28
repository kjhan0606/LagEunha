# Native Kepler execution recovery — 2026-10-28

## Execution path

Default Codex shell execution fails before starting Bash because its configured
`codex-linux-sandbox` executable is missing. Reviewed `exec_command` calls with
`sandbox_permissions="require_escalated"`, `/bin/bash`, `login=false`, and the
explicit LagForce working directory succeed on `grammar`. The full shared
context was read again after recovering access. Simulations still run through
Slurm on `syntax`; no global context or system configuration was modified.

## Recovered results, 22:15–22:29 KST

| Job | State | Authoritative evidence |
|---|---|---|
| 438626 | COMPLETED | Exact t=2.2, all 15 checks pass; geometry-CFL 0.05 crossed the former t=2.145 positivity failure. |
| 438092 | FAILED | t=2.6549183918680868, step 23008; four halo-clipped cells, zero other geometry-invalid counts. |
| 438416 | FAILED | Same time, step and halo-only guard as 438092. |
| 438533 | FAILED | t=2.7548867117105607; four halo-clipped cells, zero other geometry-invalid counts. |
| 438957 | RUNNING | A100 PCIe syn103; resumes accepted t=2.654918 state with halo 40 and has crossed the failed stage. |
| 438959 | RUNNING | Two-A100 MPI comparison; syn102, confirmed at 22:35 KST. |

The `normal` QoS caps total GPU allocation per user at eight; the MPI job
initially waited on that cap and started once capacity became available.
The H200 partition additionally has a
six-GPU group cap. Original H200 MPI job 438394 remains queued; other jobs
were preserved.

The [t=2.2 receipt](KEPLER_Q1_AREA_CFL_GEOM005_T2P2_LONG_GATE_AUDIT_2026-10-28.json)
reports disk density error 5.24809e-4, entropy-log error 1.71842e-4, orbit
position error 8.63593e-8 and global relative energy change −4.46291e-8.
Fixed mass, W=0, positive internal energy and volume, and zero recorded
floor/dual-energy/source-floor injection pass. This is a finite-interval result.

## Accepted restart and automatic continuation

Source inspection confirms `khdiag.step023008.start_geometry.dat` contains
the accepted RK4 step-start state before the failed geometry stage. The
438416 file has SHA-256
`c34a846dedaf2e01a8e5548de2aa091f98db53f761424ce6cd931107ff0fc4a4`.
All 147,456 cells and native metadata are complete, and the existing 15-check
gate passes against the original t=0 state.

The 438092 and 438416 accepted states have identical positions, mass, volume,
internal energy, density, pressure, W and y velocity. Five x velocities differ
by at most 1.110223e-16, at most two float64 ULPs. The
[staging receipt](KEPLER_Q1_HALO40_RECOVERY_STAGE_2026-10-28.json) records this
explicitly. At the accepted restart, disk density error is 2.62246e-3,
entropy-log error 1.97615e-4, position error 1.25185e-7 and global relative
energy change −2.79093e-8.

Recovery retains the same clip-swap native CUDA binary, SHA-256
`6718551061d70e24e09f0fe000a35b5cf1ac8c1a15b579f1e97e6d907b42baa0`.
Only `GFS_KEPLER_HALO_CELLS=28→40` changes. Original RK4, force/energy
bookkeeping, fixed masses, W=0, no steering and inviscid settings are retained.
Transported intrinsic cell spin is still absent.

Run root:
`/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/halo40_from_t2654918_a100_20261028`.
Job 438957 targets exact t=3. Later endpoints are t=4, 5, 6, 6.5, 7, 8, 9
and 10. Each job submits its successor only after the native 15 checks and
nine additional checks for three initial-radius disk bins pass: mean
absolute density error <1%, mean absolute entropy-log error <0.01 and mean
position error <0.01. Every segment preserves 0.5-unit checkpoints and its
exact endpoint, approximately 1.5 GB total through t=10. Each requests one
A100, two CPU threads, 32 GB and 48 h. Preparation is reproducible with
[the staging script](stage_kepler_q1_halo40_recovery_20261028.py).

The 22:29 shared-log read had passed t=2.65695 and 25 steps beyond the former
failure. Halo-clipped counts, nonpositive internal energies and floor
injections remained zero. The t=3 and t=6.5 gates remain pending.

At 22:35 KST, Slurm confirmed both 438957 and 438959 RUNNING. Recovery
438957 reached t=2.66689, global step 23157, 149 accepted steps beyond the
former failure. Its restart-relative energy ledger is −5.806e-12, with zero
nonpositive internal energies, floor/dual-energy/source-floor injection and
halo-clipped cells. The [live observation receipt](KEPLER_Q1_HALO40_LIVE_CHECK_2026-10-28.json)
preserves the exact diagnostic lines; it is not an endpoint audit.

## Independent MPI comparison

The [MPI receipt](KEPLER_Q1_MPI_A100_PAIR_STAGE_2026-10-28.json) records the
original native binary, parameters, restart and audit hashes. Job 438959
runs one and two ranks sequentially on the same two-A100 allocation, two
OpenMP threads per rank, to exact t=1.55. Both branches receive native and
radial audits. The final comparison records field differences and measured
wall-time speedup against limits fixed before execution. The original H200
comparison is preserved. Preparation is reproducible with
[the MPI stager](stage_kepler_q1_mpi_a100_pair_20261028.py).

No MPI speedup or rank-consistency result is yet available. Accurate t>6
evolution, an efficient CUDA+MPI path, a consistent intrinsic-spin treatment
where warranted and the published GIZMO long endpoint remain unverified.

## User-requested queue cleanup, evening update

The user requested removal of unnecessary submitted work. The
[cleanup record](SLURM_CLEANUP_2026-10-28_225639/README.md) supersedes earlier
queue notes: 34 pending Kepler jobs were archived and cancelled, including
the failed-parent q2 long chain, old H200 MPI/KDK variants, and the optional
cap16 probe and its follow-up audit. All 34 cancellations are confirmed with
zero elapsed allocation time. Three running jobs remain: native recovery
438957, A100 MPI comparison 438959 and inner-resolution control 437772.
Its immediate continuation 437779 is the sole pending Kepler job. Existing
sources, snapshots, receipts and scripts are preserved; other projects were
excluded. The active native recovery still advances via its guarded
successive segments.


## 2026-10-29 — measured MPI result and late-checkpoint acceleration gate

The completed 438959 comparison supersedes its earlier pending status: both
one- and two-rank runs reached exact t=1.55 and passed the native and pair
checks. Wall times were 893.109 and 821.519 s, respectively (1.08714x
speedup, 54.36% parallel efficiency). Adding GPUs alone is insufficient.
See `KEPLER_Q1_MPI_A100_PAIR_RESULT_2026-10-28.json` for field differences.

The single focused padding-reuse benchmark was revised before any allocation
started. Unstarted job 439001 was replaced by 439029, using the same late
accepted t=2.6549183918680868 checkpoint, native candidate executable and
halo=40. The seven cases include the canonical-executable qualifier, one-
and two-rank reuse off/on pairs, and one extra single-rank case enabling the
existing `GFS_GPU_LOCAL_PERIODIC_SYNC=1` device alias refresh. Only that
extra case changes the alias flag. The accepted checkpoint has all 147456
IDs in identity order, satisfying the existing alias kernel's indexing
precondition; runtime image/state checks and fallback remain in the kernel.
The late-state runtime result is not available yet. Initial t=.001 parity
and wrapper speedup are historical feasibility evidence, not promotion.

The revised auditor requires exact field bytes, identical printed accepted
step traces and per-stage tessellation/force counts for each same-rank pair,
distinct devices per MPI rank, existing native physics gates, and zero
negative-energy/floor/dual-energy/source-floor events in every ledger. The
extra device-alias case must actually use that path without host fallback.
Nine targeted acceptance/rejection tests pass. See
`KEPLER_PADDING_REUSE_LOCAL_SYNC_STAGE_2026-10-29.json` for scripts, hashes
and job replacement provenance. Seven endpoint snapshots require about
595 MB; the extension adds one approximately 85 MB snapshot and no extra
queued job. The first archive command incorrectly used `scontrol show
batch_script`; it stopped before changing the auditor or cancelling the
job. Corrected `scontrol write batch_script` preserved the queued script
and the replacement completed. Original long job 438957 was untouched.

No acceleration path has been promoted. t>6, intrinsic-spin dynamics and
GIZMO's t=120/t=600 density-field comparison remain unverified; retain the
paper's longer approximate 1000-orbit prose target for explicit reference-
radius interpretation, as required by the benchmark comparison contract.
