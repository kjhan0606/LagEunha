# Kepler research snapshot — 2026-10-23

This archive preserves the existing CUDA + MPI + OpenMP work and the research leading to the current full-space cell-condensation/reference-cache candidate. The long-time Kepler stability goal remains unachieved. Read [the current handover](../../Docs/HANDOVER_2026-10-23.md) first.

## Layout

- `../../Exam/`: modified original solver and new physics headers.
- `stream_candidates/`: original GPU-stream, historical t8 and radial-resolution scripts, reports and small data summaries.
- `rk4_diagnostics/`: development source, native test programs, build/run wrappers, analyses, plans, reports and small mathematical data bundles. Historical alternatives and rejected candidates are retained; presence in this directory is not endorsement or evidence of a successful run.
- `frozen/`: exact source/include closures for cell structure, executed factor revision2, native readback, successful actual RHS and capacity-rejected reference continuation. Use these closures to identify the exact source of a reported result. Development headers can differ from these snapshots.
- `evidence/`: small original manifests, native result JSON, status/exit records and integrity logs. Older result fields marked pending refer to their own point in time, not to the newest state. Original manifests intentionally retain GPFS absolute paths.
- `SNAPSHOT_SOURCE_MANIFEST.json`: every copied file's original path, repository path, size and SHA256. Generated documentation and this index are not original calculation inputs.

## Established versus pending

The latest accepted full-space state is step5 at t=0.00011263377801805954. The next actual k1 was computed and saved; no new step was accepted by the cell-condensation comparison. Full original-M residuals and all 71,080 source-cell energy/angular gates passed for that RHS. The measured solver comparison is not a full RK4 speed measurement.

The reference-cache continuation built for Ada and H200 but exited2 on Ada before cache generation because the original20% free-memory requirement was not met. The old1024-triangle metadata estimate is invalid; native code uses4096. See `evidence/material_cell_reference_continuation_capacity_failure_v1/corrected_memory_plan.json`. The approved buffer-lifetime change is a next implementation task, not completed work.

## Build and run provenance

The frozen sources and archived `run_env.sh`/runner scripts record the actual Intel MPI, GNU C++ and CUDA environment. GPU builds used sm89 and sm90 with FMA disabled; host builds retained binary80 arithmetic with fast math and contraction disabled. This snapshot does not change those arithmetic contracts.

The archived wrappers contain original absolute run paths and are historical records, not relocatable one-command launchers. Rebuild and verify source/input bindings in a new run directory before use. Never overwrite accepted outputs or consume the saved next k1 in two competing continuations.

Large quadrature, basis, mass, factor, source-state and restart binaries remain at:

```text
/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics
```

Those data and their original checksum manifests are required to reproduce the native calculations. They are deliberately excluded from Git, together with executables, object files, Fortran module build artifacts and large per-cell tables. Existing tracked repository files outside this snapshot are preserved.
