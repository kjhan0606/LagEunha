# Kepler Slurm cleanup — 2026-10-28

The user explicitly requested removal of unnecessary submitted work. All
34 selected jobs are confirmed `CANCELLED by 10396`, each with elapsed time
`00:00:00`; none consumed an allocation. No selected cancellation target
remains in the queue. Source, run files, checkpoints and audit results are
preserved. Jobs from other projects were excluded by ownership and project
path checks; job 436751's home-directory workdir was resolved by its actual
LagForce batch-script paths before inclusion.

## Retained work

| Job | Verified state | Purpose |
|---|---|---|
| 438957 | RUNNING | Current native 384² W=0 long-time trajectory, recovered at t=2.654918; validated segments progress toward t=6.5 and t=10. |
| 438959 | RUNNING | Matched 1-rank/2-rank A100 correctness and throughput comparison. |
| 437772 | RUNNING | Independent 420² inner-resolution control. |
| 437779 | PENDING | Continuation of that resolution control toward t=6.5, with restart-input checks. |

## Removed pending work

| Reason | Jobs |
|---|---:|
| q2 long chain descending from failed parent 437250 | 18 |
| t=6.5 RK4/KDK pair whose parent 438092 failed, and its successor | 3 |
| Superseded q1 trajectory and optional H200 execution, integrator or scaling variants | 7 |
| Deferred q2 MPI/OpenMP variants | 4 |
| Unstarted capacity probe 438420 and follow-up CPU audit 438970 | 2 |
| **Total** | **34** |

The optional experiments remain reproducible from their preserved files.
Reconsider them after the active comparison and physics gates identify a
specific need. Keep immediate prerequisites in the queue; the current 384²
trajectory submits only its next segment after its acceptance gates pass.
The original long-time Kepler, Laguerre-extension, spin-consistency and
CUDA+MPI objectives remain open.

## Evidence

- [Decisions and original scheduler records](before.json).
- [All original cancelled batch scripts](batch_scripts).
- [Cancellation requests](cancellation_requests.json).
- [Authoritative post-cleanup accounting and queue](after.json).
- [Machine-readable summary](summary.json).
- [Exact cleanup script used](cleanup_script_executed.py).
