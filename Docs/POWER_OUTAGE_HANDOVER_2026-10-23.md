# Kepler power outage handover — 2026-10-23

## User instruction and stop

Stop this project's work at 10:50 Asia/Seoul, commit and push before the announced
11:00 server shutdown, then cease work. No operating-system shutdown is performed.
Recorded at: 2026-10-23T10:50:19.177690+09:00.
Only owned processes/jobs bound to `/gpfs/kjhan/LagForce/Kepler_` are selected.
No new simulation, restart, requeue or automatic continuation is authorized after
this stop until the user explicitly resumes. Preserve partial files and all raw data.

Stop observations (empty lists mean no matching work was running):

```json
{
  "local": {
    "before": [],
    "events": [],
    "after": [],
    "all_project_workers_stopped": true
  },
  "slurm": {
    "before": [],
    "events": [],
    "after": [],
    "all_project_jobs_stopped": true
  }
}
```

## Latest native continuation, superseding earlier pending/resource-wait reports

Job 435137 completed two new accepted native steps (6 and 7), with independent
actual v3 restart readbacks. Latest accepted time is
`0.00019779036422122955`; total accepted steps
`7`.
The prepared genuine step-6 selected k1 was consumed exactly once; the old full
shadow and accepted intervals were not repeated. Do not reuse that consumed k1.

Latest accepted rank checkpoints remain on shared GPFS:

```json
[
  {
    "path": "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/material_cell_reference_continuation_v1_h200_mpi4/run/continuous_step7/accepted_r0.bin",
    "bytes": 86580269
  },
  {
    "path": "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/material_cell_reference_continuation_v1_h200_mpi4/run/continuous_step7/accepted_r1.bin",
    "bytes": 86580269
  },
  {
    "path": "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/material_cell_reference_continuation_v1_h200_mpi4/run/continuous_step7/accepted_r2.bin",
    "bytes": 86580269
  },
  {
    "path": "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/material_cell_reference_continuation_v1_h200_mpi4/run/continuous_step7/accepted_r3.bin",
    "bytes": 86580269
  }
]
```

Complete step-7 wall time including readback was
`38.4340076` seconds.
Long release: `False`.
Stop reason: `measured_complete_step_exceeds_24h_t1_or_7day_t10_budget`.
The measured cost gate failed; t>6–10 Kepler stability has NOT been demonstrated.
Global map injectivity guard is also incomplete. Historical step-5/pending files
retain their observation times and do not override this actual continuation result.

Before resuming, inspect latest checkpoint/source/operator bindings and any
geometry precomputation binding. Do not silently replay accepted intervals or
apply old-state geometry to the accepted step-7 state. Preserve the consumed
`material_step6_reference_continuation.claim` and its terminal record.

## Code and data preservation

Authoritative research root: `/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics`.
Original parallel solver source: `/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src`.
Isolated publication worktree: `/gpfs/kjhan/LagForce/Kepler_publish_20261023` (Git commands run on syntax).
Remote: `git@github.com:kjhan0606/LagEunha.git`, intended branch `master`.

The commit contains source, small JSON/log evidence, exact frozen source closures,
and this handover. Native binaries, large quadrature/geometry/factor data and restart
checkpoints stay on GPFS; they are not deleted or uploaded as Git objects.
`Research/Kepler/SNAPSHOT_SOURCE_MANIFEST.json` binds copied source bytes.
`SNAPSHOT_CONTENTS.sha256` indexes the archived files.

The final commit/push hash and verification receipt are written AFTER publication
to `/gpfs/kjhan/LagForce/PowerOutage_20261023/publish_receipt.json`; this document does not claim a future push
has already succeeded. Read that receipt and `/gpfs/kjhan/LagForce/PowerOutage_20261023/result.json` before resuming.
No new physical validation was launched by the shutdown procedure.
