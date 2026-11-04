# Native common FV--ALE checkpoint

This archive accompanies `Docs/UNIFIED_HYDRO_FV_NATIVE_2026-11-05.md`.
It contains a CPU/OpenMP prototype using the original native Voro library,
not full production MPI/CUDA evolution or a new accepted Kepler result.
`manifest.json` records lossless source and stored-file SHA256 hashes.
Native executables remain in the authoritative GPFS run root; build/source
manifests identify their exact binaries and inputs.

## Evidence groups

- `fv_native_build_v2` and `fv_native_coupling_v1`: initial implementation
  using historical geometry policies. The two multi-dimensional Sod
  failures are retained; their incomplete endpoints are not passes.
- `fv_native_build_v3` and `fv_native_coupling_v2`: corrected explicit
  true-power geometry, ten complete native contact/Sod coupling cases.
- `fv_native_sod_large_v1`: complete 2D64/128 and 3D16 Sod to t=0.2,
  with exact-reference actual-cell comparisons.
- `fv_power_geometry_gate_v1`: retained original raw derivative gate
  failure at the degenerate extruded fixture.
- `fv_power_geometry_gate_v2`: fixed three-increment diagnostic and two
  Richardson pairs; the original raw failure remains recorded.
- `fv_sod_reference_v1`: native exact Riemann integration source/build
  provenance. This is an analytic reference, not another PDE evolution.
- `fv_push_default_geometry_regression_v1`: native parent/current legacy
  W=0 geometry checks, both passing with byte-identical output.

Files named `running.json` retain historical launch PIDs. They do not
indicate a job currently running. Case `result.json` and `process.json`
give terminal status. Binary `.fv.gz` checkpoints, source, logs and text
data are gzip-compressed without rewriting their contents.

## Replay from the published source tree

Use a new, previously nonexistent run directory under the authorized GPFS
project root. An Intel oneAPI environment providing `mpiicx` is required
to reproduce the recorded build configuration. From the repository root:

```sh
set -eu
REPLAY_ROOT=/gpfs/kjhan/LagForce/FVNative_replay_new
mkdir "$REPLAY_ROOT"
export TMPDIR="$REPLAY_ROOT"
mpiicx -O2 -qopenmp -std=gnu11 -DXYZDBL -DGOTPM -Wall -Wextra \
  Exam/gfs_fv_native.c Exam/gfs_fv_geometry.c \
  Exam/HydroCore/fv_native_evolution.c Voro/voro.c -lm \
  -o "$REPLAY_ROOT/fv_native"
python Exam/HydroCore/run_fv_native_suite.py \
  "$REPLAY_ROOT/coupling" "$REPLAY_ROOT/fv_native" --phase coupling
python Exam/HydroCore/run_fv_native_suite.py \
  "$REPLAY_ROOT/sod_large" "$REPLAY_ROOT/fv_native" --phase sod_large
mpiicx -O2 -qopenmp -std=gnu11 -DXYZDBL -DGOTPM \
  Exam/HydroCore/exact_sod_cell_reference.c -lm \
  -o "$REPLAY_ROOT/sod_reference"
python Exam/HydroCore/assess_fv_sod.py \
  "$REPLAY_ROOT/sod_large" "$REPLAY_ROOT/sod_reference"
```

The scheduler requests eight OpenMP threads per case and runs up to two
cases at once. Its current binding does not guarantee disjoint cores.
Recorded times are not matched or isolated GIZMO benchmarks. Gresho/cold
shear, the remaining new-operator test suite, temporal GCL convergence,
restart replay, full native GPU/MPI and <=2x GIZMO remain outstanding.
The manuscript remains on its separately verified production baseline.
