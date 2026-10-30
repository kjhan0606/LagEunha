# Kepler wall-validation OpenMP candidate

Base: `ded504a6624d16623a5fd05e4dbf69943fc0fe83` (`parallelize independent RK4 particle loops`).

## Change

`GFS_WALL_PARALLEL_VALIDATION=1` enables an OpenMP scan in `Exam/gfs_wall_events.inc:gfs_wall_valid`. Each thread checks one trial cell and atomically marks its accepted global-ID slot. A reduction collects invalid-state flags. If the scan finds any invalid ID, duplicate, mass change, or invalid state, the code increments the generation and reruns the original serial validator so the first-failure reason and rejection path remain deterministic. With the environment flag absent, validation follows the existing serial path.

The optimization is limited to the one-rank W0 wall-event path. It does not remove any validation, change the physical RHS, make the event transaction MPI-safe, or establish an OpenMP speedup.

## Current evidence

- Host object compile with Intel MPI C compiler and `-qopenmp`: passed on 2026-10-30.
- The object references the Intel OpenMP runtime and contains the new environment-string branch.
- Only two existing warnings were emitted: the `ptrdiff_t` format warning at `Exam/exam.c:8274` and the pointer-type warning at `Exam/exam.c:11007`.
- No linked GPU executable or runtime test has been made for this candidate.

## Required qualification

1. Link on a node with the matching CUDA 13.0.2 runtime, reusing only GPU objects whose source hashes match the baseline build.
2. From the same accepted t=14 checkpoint, compare the original validator and this candidate through the same exact endpoint. Keep identical `OMP_NUM_THREADS`, wall-event checks, and all physics settings.
3. Require bitwise equality of every stored particle field after sorting by global ID, successful boundary/rollback/conservation gates, and a measured reduction in median and p95 wall time per accepted RK4 step.
4. Exercise an invalid-state or duplicate-ID rejection fixture and confirm the serial fallback reports the same first-failure reason. Do not enable this path for the long run before these gates pass.
