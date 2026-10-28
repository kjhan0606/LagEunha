# Q1 bounded-domain experiment, 2026-10-29

The failed periodic Q1 t=3 state showed exterior sites entering disk support. This experiment changes the boundary model, using the native fixed 1024-gon at R=2.75 and retaining original t=0 generators with r<2.65. It starts a new physical IC at t=0; it does not modify the contaminated t=3 state or discard any previous checkpoints/results. The model's ultimate equivalence to the GIZMO long-time test still requires boundary studies.

## Initial-state evidence

`INITIAL_AUDIT.json` pins the canonical full Q1 t=0 source and the new 90,352-site native restart. All 48,680 original disk IDs, positions, velocities, masses and thermal energies are exact. Every complete core/disk record byte is preserved. Only the exterior is initialized once with its new bounded volume, original Hopkins density profile and uniform pressure; its masses then remain fixed. W=0, Re=infinity, no steering and no intrinsic spin are retained.

Two independent guard-radius constructions agree to 1.18e-10 in relative cell area. Their total area differs from the exact polygon area by 9.95e-14. Original disk float volumes agree with the independent bounded construction to 3.03e-8; the minimum conservative wall clearance is 0.1001557. The full initial restart is about 52 MB and remains on GPFS. The native restart argument is 1 solely because the legacy argument 0 generates a fresh IC; the stored physical time is exactly zero and this label represents no elapsed evolution.

`NATIVE_PREFLIGHT.json` checks the actual GPU first RHS for all owned sites and an accepted step to t=1e-9. All native volumes match the independent mesh to 5.72e-8. Eight selected core/disk/exterior cells pass independent force/work checks; the largest acceleration discrepancy is 3.64e-27. Masses/IDs remain fixed, states remain positive and no floor/dual/entropy injection occurs. This is an initialization check, not a short-time substitute for the requested long evolution.

One independent polygon contained a 2.20e-16 edge from co-circular degeneracy. Its direction is numerically unresolved. The audit reports its length and a conservative acceleration/work contribution bound, adds that bound to the error budget and leaves every native face unchanged. The acceleration bound is 3.48e-18, below the predeclared 1e-12 independent-force threshold. No physical gate was loosened and the completed native calculations were not rerun.

## Active continuation

Slurm job **439333**, `KepQ1BoundEvents`, was confirmed RUNNING on `syn104`, with `--partition=h200 --gres=gpu:H200:1`, eight CPUs/OpenMP threads and 32 GiB. Its first native accepted steps were observed advancing with zero injections. The exact source/binary and stage inputs are pinned in `CHAIN.json`; `job.slurm` records the submitted resource request.

One allocation executes the dependent sequence `.825 -> 2.5 -> 3 -> 6.5 -> 10`, reading each preceding accepted checkpoint. Only five exact endpoints are requested, about 260 MB total, plus small logs/audits. Periodic snapshots and face/tessellation dumps are disabled. The driver refuses existing output rather than replaying a stage. `.825` and `2.5` are diagnostic comparison checkpoints; invalid state, IDs, mass, injection or geometry stops immediately. Exact t=3 and later must pass all original 15 physics thresholds and unchanged 9 radial thresholds; the area invariant correctly uses polygon area instead of periodic area 36. Extra ID/population, endpoint-domain and rollback checks are retained. A failed gate stops the chain without altering thresholds.

The actual root is `/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/chain_h200_v2`; inspect live Slurm state, `chain_state.json` and the active native `log`. A stale state file alone is not evidence of a live job. If interrupted, inspect the last complete native checkpoint and resume from it; do not rerun t=0. The existing independent job 437779 was left running and pending job 439029 was untouched.

No Q1 t=3/6.5/10 result is asserted by this record. After those gates, boundary dependence, runtime cost and GIZMO t=120/600 and the interpreted approximately 1000-orbit target remain required. Preserve the full long-time objective.

## Diagnostic-setting correction during the same process

The initial launch mistakenly retained `GFS_RK4_AUDIT_STEP_MIN=0` and `MAX=1000000` after removing its selected-ID/stage filters. This enabled all-particle TSV emission at every RK stage. The 0.9-second CPU summary did not represent the full observed time per step. The error was found by comparing scheduler elapsed time with actual accepted progress; the generated TSV is auxiliary data, not a restart checkpoint.

The same MPI PID 3947439 was retained. Within its own Slurm allocation, the debugger verified the exact native binary/process identity and changed only `gfs_rk4_stage_diagnostic.audit_min` from 0 to -1. No physical field, simulation time, timestep, executable or checkpoint changed. Parent/solver processes were continued, and accepted native progress resumed. Subsequent stage scripts now unset both range variables; their old scripts/hashes and the original chain configuration remain preserved. `DIAGNOSTIC_CORRECTION.json` and the debugger transcript record the change. This is a corrected runtime diagnostic setting, not a restart from t=0.

Large redundant all-particle TSV files from this mistake may be reduced to bounded header/tail samples with a size/provenance receipt once their writes have stopped. Every native initial/final checkpoint, selected-cell preflight TSV, physical audit and native log is retained. The effective chain configuration records the corrected diagnostic settings.
