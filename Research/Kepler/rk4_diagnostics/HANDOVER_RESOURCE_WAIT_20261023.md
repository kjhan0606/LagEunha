# Resource wait handover — 2026-10-23

The long-time Kepler goal is unachieved. Three consecutive goal turns encountered the same external execution-resource block. The last turn before this audit was a verified wait on the specific live Slurm job435137, not a failed or missing process.

## Authoritative current state

- Job435137 remains PENDING/Resources with no allocated node and zero native runtime. The scheduler predicts2026-10-24T04:44:42; this is not an allocation or a guaranteed start time.
- Its original H2001/MPI4/OMP8/CPU32/192G request, frozen source/binary and runner remain intact. No cancellation or resubmission was done.
- Ada has live geometry probe processes using19,270MiB of32,760MiB. Samples418197,482095,482150 were confirmed by ps. Other-thread processes were not stopped.
- The approved resource decision preserves the existing H200 request. The H100 test-only estimate was later; no H100 replacement was submitted. Ada lifetime v2 built successfully but was not executed.
- The shared step6 claim and the H200/lifetime output run directories do not yet exist. The actual accepted state remains step5 at t=0.00011263377801805954. Its genuine saved next k1 remains available.

## Resume from evidence

Read the SAME job435137 through squeue/sacct first. If it has begun, inspect its live handle/logs rather than launching a replacement. If it has ended, inspect native exit, reference/cache comparison, per-stage full-M/cell gates, accepted checkpoint and restart/cost results. The queued script can run independently of the Codex goal status when Slurm assigns resources.

Do not infer a stopped job from an observation timeout, a lock file or an old status file. Do not clear the shared claim automatically. If the job accepted step6 before a later failure, preserve and resume that accepted state rather than repeating step5/k1. Do not run competing Ada/H200 continuations.

The whole-step cost, map guard, t6/8/10 stability and longest GIZMO comparison remain unproven. No numerical tolerance, rule, arithmetic, tile or thermal context was relaxed during this resource wait.
