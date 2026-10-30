# Kepler wall-validation OpenMP candidate

Base: `ded504a6624d16623a5fd05e4dbf69943fc0fe83` (`parallelize independent RK4 particle loops`).
Candidate code: `dd8800f493be1148b91d5c0b7b70c6c11a3ca449`.

## Change

`GFS_WALL_PARALLEL_VALIDATION=1` enables an OpenMP scan in `Exam/gfs_wall_events.inc:gfs_wall_valid`. The actual scan loop is shared with `Exam/gfs_wall_validation_scan.h`, which the unit harness compiles directly. Each thread checks one trial cell and atomically marks its accepted global-ID slot. A reduction collects invalid-state flags. If the scan finds any invalid ID, duplicate, mass change, or invalid state, the code increments the generation and reruns the original serial validator so the first-failure reason and rejection path remain deterministic. With the environment flag absent, validation follows the existing serial path.

The optimization is limited to the one-rank W0 wall-event path. It does not remove any validation, change the physical RHS, make the event transaction MPI-safe, or establish an OpenMP speedup.

## Current evidence

- Host object compile with Intel MPI C compiler and `-qopenmp`: passed on 2026-10-30.
- OpenMP scan harness: `Exam/tests/run_gfs_wall_validation_scan.sh` passed with 1, 4, and 8 threads. Cases cover unique shuffled IDs, duplicates, out-of-range and missing IDs, changed mass, invalid thermal and RHS state, nonzero W0 weight, and a valid new generation after rejected cases.
- The harness tests the shared scan loop, but not the RK4 caller's serial fallback reason or a full simulation run.
- A CPU-only synthetic scan microbenchmark at 90,352 cells measured 2.86x speedup for increasing IDs and 2.81x for reversed IDs, but only 0.45x for randomly shuffled IDs. The exact t=0, t=10, and t=14 native global-ID sequences have now been compared and are identical: 99.6248% of adjacent IDs increment by one and 99.7477% differ by at most 64. A source audit counts 23 wall-validation calls in the active RK4 routine (8 with RHS checks, 15 without). Applying the RHS-enabled ordered-ID microbenchmark to all 23 gives a nominal 1.89% per-step saving; the production mix may save less. At the observed 10.49-day continuation projection, this corresponds to 4.75 wall hours before the cost of qualification. The two one-unit A/B variants are estimated at 0.86 simulation hours. This derived estimate justifies the same-binary H200 A/B after the t=14→17.8 physics gate if the queue permits without delaying the reference chain; it is not end-to-end evidence. Detailed inputs and arithmetic are in `/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/wall_validation_same_binary_ab_t14_t15_20261030_v1/cpu_microbench_20261030.json`, `restart_id_order_t0_t10_t14_20261030.json`, and `wall_scan_callsite_estimate_20261030.json`.
- Linked GPU executable: `/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/wall_validation_omp_build_20261030_v2/eunha2_wall_validation_omp.exe`.
- Executable SHA256: `594dd70b39830ee953e38ae38ad99d5fe0d26ba08e3df9af483dcefcbcd2d07c`. Link receipt and dependency check are in the same build directory; `ldd` found no unresolved shared libraries.
- The build reused the exact-hash-matched `exam_gpu.o` and `exam_gpu_extract.o` artifacts from the existing OpenMP build. Source comparison confirms neither GPU source file changed between the base commit and this candidate.
- The executable contains the new environment-string branch and links `libiomp5` plus `libcudart.so.13`.
- Two pre-existing compile warnings remain: the `ptrdiff_t` format warning at `Exam/exam.c:8274` and the pointer-type warning at `Exam/exam.c:11007`. The Fortran link driver also warns that it ignores the C-only `-fno-strict-aliasing` flag.
- No runtime or state-parity test has been made for this candidate.

## Required qualification

1. Compare the base executable and candidate's default serial path from the same accepted t=14 checkpoint through exact t=15. Require all physics gates and bitwise equality of every saved field after matching global IDs.
2. After the t=14→17.8 physics gate, compare this candidate with GFS_WALL_PARALLEL_VALIDATION unset and set to 1, using the same executable, t=14 checkpoint, parameters, environment, and H200 allocation, provided the queue does not materially delay the reference chain. The input ordering is already recorded in `/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/wall_validation_same_binary_ab_t14_t15_20261030_v1/restart_id_order_t0_t10_t14_20261030.json`. Require bitwise endpoint equality and passing physics gates, then measure median and p95 RK4 wall time. The candidate qualifies as a speed optimization only if this end-to-end measurement shows a benefit.
3. Exercise an invalid-state or duplicate-ID rejection fixture and confirm the serial fallback reports the same first-failure reason. Do not enable this path for the long run before these gates pass.


## Continuation verification (2026-10-30)

Reran Exam/tests/run_gfs_wall_validation_scan.sh with OMP_NUM_THREADS=1, 4, and 8. All three invocations returned PASS: OpenMP wall-validation scan cases. This revalidates the shared scan kernel and its rejection fixtures; it does not establish integrated RK4 parity or an end-to-end speedup. Those still require an H200 run.
