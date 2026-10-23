# Scheduled power outage — 2026-10-23

This directory is the Kepler research root. The original source and Git publication worktree are separate; confirm their identity before changing them.

The user instructed that this project's work must stop at **2026-10-23 10:50 Asia/Seoul**, followed by commit and push, before the announced server shutdown at11:00.

- At or after the deadline, do not launch, restart, requeue or automatically resume simulations until the user explicitly resumes work.
- Preserve completed checkpoints, partial outputs, immutable sources and restart bindings. Do not delete calculation data as part of stopping.
- Stop only this project's jobs and worker processes; other projects are outside this instruction.
- One-time control files and completion receipts are in `/gpfs/kjhan/LagForce/PowerOutage_20261023`.
- Read `POWER_OUTAGE_STOP_20261023.json` and the final handover/status before any post-outage continuation. A previously pending Slurm job may already have completed; use current scheduler/output evidence.

The power outage does not establish the scientific goal. Long-time t6/8/10 and the longest GIZMO comparison remain subject to their existing physical, conservation, restart and cost checks.
