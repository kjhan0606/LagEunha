# Kepler wall-validation OpenMP candidate

Base: `ded504a6624d16623a5fd05e4dbf69943fc0fe83` (`parallelize independent RK4 particle loops`).

## Change

`GFS_WALL_PARALLEL_VALIDATION=1` enables an OpenMP scan in `Exam/gfs_wall_events.inc:gfs_wall_valid`. Each thread checks one trial cell and atomically marks its accepted global-ID slot. A reduction collects invalid-state flags. If the scan finds any invalid ID, duplicate, mass change, or invalid state, the code increments the generation and reruns the original serial validator so the first-failure reason and rejection path remain deterministic. With the environment flag absent, validation follows the existing serial path.

The optimization is limited to the one-rank W0 wall-event path. It does not remove any validation, change the physical RHS, make the event transaction MPI-safe, or establish an OpenMP speedup.

## Current evidence

- Host object compile with Intel MPI C compiler and `-qopenmp`: passed on 2026-10-30.
- Linked GPU executable: `/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/wall_validation_omp_build_20261030_v1/eunha2_wall_validation_omp.exe`.
- Executable SHA256: `0d7595d03699a6b904d7cdcb6136212ef141ae22456aab3f742e3aff92b13817`. Link receipt and dependency check are in the same build directory; `ldd` found no unresolved shared libraries.
- The build reused the exact-hash-matched `exam_gpu.o` and `exam_gpu_extract.o` artifacts from the existing OpenMP build. Source comparison confirms neither GPU source file changed between the base commit and this candidate.
- The executable contains the new environment-string branch and links `libiomp5` plus `libcudart.so.13`.
- Two pre-existing compile warnings remain: the `ptrdiff_t` format warning at `Exam/exam.c:8274` and the pointer-type warning at `Exam/exam.c:11007`. The Fortran link driver also warns that it ignores the C-only `-fno-strict-aliasing` flag.
- No runtime or state-parity test has been made for this candidate.

## Required qualification

1. From the same accepted t=14 checkpoint, compare the serial validator (`ded504a`) and this candidate through the same exact endpoint. Keep identical `OMP_NUM_THREADS`, wall-event checks, and all physics settings.
2. Require bitwise equality of every stored particle field after sorting by global ID, successful boundary/rollback/conservation gates, and a measured reduction in median and p95 wall time per accepted RK4 step.
3. Exercise an invalid-state or duplicate-ID rejection fixture and confirm the serial fallback reports the same first-failure reason. Do not enable this path for the long run before these gates pass.
