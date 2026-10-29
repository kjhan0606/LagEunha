# One dependent accepted-restart bundle after t=10

Slurm job **439375** (`KepQ1RestartDt`) requests one H200, eight CPUs/OpenMP threads, 32 GiB and two hours. Its actual dependency is **afterok:439333**, with automatic cancellation if that dependency becomes invalid. The parent native chain is unchanged. This is one finite computational allocation, not a polling daemon.

The worker additionally requires the exact accepted t=10 endpoint, all original 28 physical checks, matching parent job/parameters/executable, complete native IDs and fixed masses, and matching snapshot hashes. It refuses an incomplete parent or an earlier green checkpoint. The local preflight correctly reported NOT READY and created no native execution claim or case directory. Slurm `--test-only` accepted the request; its hypothetical job/start prediction was not an actual submitted run or an allocation guarantee. The actual job was confirmed PENDING for Dependency with StartTime unknown.

## Finite sequence

1. Four original/cached/cached/original runs start from the same accepted t=10 state on H200, over approximately eight recent parent timesteps. Full restart record/footer parity, exact printed numeric ledgers, state/physics checks and at least 1.05x ABBA median speedup qualify the cache on this hardware. If the completed qualification report fails, controls retain the original v2 binary. A native failure or incomplete qualification stops with data preserved.
2. From that same original t=10 state, three native controls advance exactly **0.1 physical time** with original CFL, a half cap and a quarter cap based on the minimum of the recent 64 accepted parent intervals, excluding its shortened endpoint. These runs use the same selected executable, parameters and hardware. No trajectory is restored to an analytic orbit, and no mass/weight/thermal repair is applied.
3. Each control must pass the unchanged physical checkpoint gates. The comparison requires genuine measured timestep separation, reports core/disk endpoint differences and incremental entropy/energy changes, and retains absolute original-IC metrics. It makes no formal convergence-order or inherited-error-repair claim and does not automatically promote a t120/t600 run.

The interval remains fixed. An allocation-cost guard stops before the three controls if the estimated work would exceed the remaining allocation; it does not silently shorten the experiment. The driver rejects existing output/claims instead of replaying t=0 or duplicating a partial execution. An interrupted worker requires an explicit audit of its real Slurm state and retained outputs before any resume.

## Output and safeguards

Four hardware-qualification endpoints plus three timestep-control endpoints total about **365 MB** of native snapshots, plus small maps/logs/audits. The accepted parent is reused by symlinks. Periodic snapshots, full-particle RK-stage TSV and face/tessellation dumps are disabled. The effective environment is checked immediately before each MPI launch, including target, cache mode, requested timestep cap and diagnostic bounds.

`GUARD_TESTS.json` records fault-injection checks covering the previous full-TSV configuration error, incompatible physics, wrong target/cap/cache, incomplete parent, and use of a real green t3 checkpoint in place of t10. `AUDITOR_REGRESSION.json` records reanalysis of existing completed GPU data after making the comparator accept a dynamic native restart filename. No gas evolution was rerun for those checks.

`bundle.json`, `bundle.sha256`, `native_template.sh`, `params.dat` and `job.slurm` preserve the actual submitted contract. Driver/auditor source is in the project-level `run_kepler_q1_restart_dt_bundle_20261029.py`, `audit_kepler_q1_restart_dt_bundle_20261029.py`, and `audit_kepler_const_config_cache_20261029.py` files; frozen runtime copies and their hashes are in `/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/post_t10_restart_dt_h200_v1`.

This is prospective work. Parent t6.5/t10 and this bundle are not asserted complete. Once the parent finishes, check actual scheduler allocation rather than a stale state file or preflight forecast. If resource waiting is long, assess the already authorized Ada alternative after verifying and cancelling only the unstarted duplicate; never start the same bundle twice. No fallback is launched by this script.

## Diagnostic-only revision after t6.5

`AUDIT_REVISION_T65.json` records the actual PENDING/Dependency check before the auditor update. The original bundle/source are retained; `bundle.before_t65_audit.json` and its checksum preserve the previous contract. The new auditor additionally reports core-adjacent/exterior distributions, cell p99/max density/entropy/orbit errors and incremental regional energy/angular momentum. Its orbit normalization divides by max(original radius, Plummer softening), including the exactly central generator. No native duration, force, physics gate or executable changed. `REGIONAL_DISTRIBUTIONS_T65.json` in the bounded Q1 evidence folder validates these new metrics using existing states; no common-restart timestep trajectories have run yet.
