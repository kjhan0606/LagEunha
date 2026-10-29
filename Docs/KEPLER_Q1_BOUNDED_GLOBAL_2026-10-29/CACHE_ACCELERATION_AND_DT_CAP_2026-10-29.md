# Q1 cache and timestep-cap follow-up, 2026-10-29

## Completed t=10 restart controls

Job 439375 restarted each case from the same exact t=10 checkpoint and advanced to t=10.1. The coarse, half-cap and quarter-cap cases all completed the original physical gates. Frozen source/restart hashes and native full-record checks passed; no masses, W=0 state, or stored restart fields changed across cache A/B comparisons except the runtime next-pointer field.

| Case | Accepted steps | Median h | Incremental core entropy absolute mean | Incremental disk entropy absolute mean | Incremental core density absolute mean |
|---|---:|---:|---:|---:|---:|
| Coarse | 258 | 3.976e-4 | 2.96e-4 | 3.01e-6 | 0.37775 |
| Half | 2,101 | 4.760e-5 | 4.23e-6 | 6.03e-5 | 0.37776 |
| Quarter | 4,202 | 2.380e-5 | 5.45e-6 | 1.076e-4 | 0.37776 |

The inherited core density discrepancy is essentially unchanged by these short continuations. Entropy increments do not decrease monotonically with the smaller step: the core increment falls after halving, while the disk increment rises. This is timestep-sensitivity evidence over only 0.1 time units, not a convergence proof or an explanation of the accumulated core error. See post_t10_restart_dt_h200_v1/controls/comparison.json and qualification/audit.json under /gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029.

## Isolated process-constant configuration cache

The source change is limited to Exam/exam.c and Exam/gfs_thermal64.h. With GFS_NATIVE_CONST_CONFIG=1, the native process reads the Kepler center/softening and thermal-mode environment once instead of repeatedly calling getenv in the hot path. The default remains off and retains dynamic reads. When enabled, GFS_NATIVE_CONST_CONFIG, GFS_THERMAL64, EUNHA_IC, and the Kepler center/softening environment values must remain fixed for the process.

Two ABBA restart benchmarks passed exact native record/footer parity apart from the runtime next pointer:

- H200 at t=10: 1.974x median wall speedup; baseline/cache medians 0.6627/0.3357 s per accepted step.
- The bounded t=.825 comparison: 1.515x median wall speedup; baseline/cache medians 0.7587/0.5010 s per accepted step.

Both compare short restart intervals and do not establish MPI scaling or long-time performance. The updated source also compiled into Exam/libexam.a with CUDA 13.0.2 and Intel oneAPI 2025.3 using two build threads. The compiler emitted pre-existing unused-variable and format warnings outside the changed lines; no compile errors occurred.

## eta=.06 original-IC branch

The frozen gravity-only RK4 screen selected eta=.06 as the largest tested cap that kept all original core sites within its preregistered orbit-position budget through t=10. The corresponding global ceiling is 7.6420e-5 and the gravity-only screen required 130,858 steps to t=10. This is not a hydrodynamic accuracy result.

A new native branch was prepared from the preserved original t=0 IC using the cache executable (SHA256 7364a05c09eeb393d4e8528f01a82efb7f6d9c6fba9b394d2d1e62e6bd78345f), eta=.06, W=0, fixed mass and the existing bounded polygon model. It saves exact accepted checkpoints at .825, 3, 6.5 and 10, with no periodic particle/tessellation dumps. Its resource request is one H200 GPU, one MPI rank, eight OpenMP threads and 32 GiB. Four checkpoints are about 208 MB total.

Slurm job 440745 was submitted to partition h200 with --gres=gpu:H200:1 and is now RUNNING on syn104. At the last log read, chain_state.json showed stage to0p825 and the native log had reached step 157 at t=0.0119215; no checkpoint gate has completed yet. The preserved run package is /gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/dteta006_cache_t0_to10_h200_v1. The next scientific decision must use its actual native checkpoints and the original t=0 core/disk metrics. This branch does not yet establish t=10 accuracy, t=120/t=600 stability, or the effect of spin/variable weights.
