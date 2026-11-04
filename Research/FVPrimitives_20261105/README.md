# Shared native dimensional FV primitives and CUDA layout controls

See `../../Docs/UNIFIED_HYDRO_FV_PRIMITIVES_2026-11-05.md` and the grouped
Astra plan `../../Docs/UNIFIED_HYDRO_FV_ALE_PLAN_2026-11-05.md`.

The authoritative final records are `v7/stdout_C.json`, `stdout_CXX.json`,
`stdout_CUDA.json`, `process_*.json`, and `memcheck/`. C/C++ cover 120 fan
fixtures plus analytic cone roots. The actual RTX 5000 Ada covers 96 mixed
runtime-dimensional inputs, with 28 HLL fallback cases. Final normalized
CPU/device discrepancy is 3.21e-15 and memory check reports zero errors.
The header uses a fixed five-slot payload with zero inactive momenta and
only d+2 independent physical variables. The new solver is not enabled.

Versions v1--v5 preserve failed packed-layout CUDA results. v3 contains an
optimization-disabled compiler control (passes) and memory check (zero
memory errors but scientific audit fails); v4 contains isolated constant/
runtime dimension and fan/cone call combinations. Forcing inlining alone
did not resolve the failed path. v6 is the first fixed-layout pass. v7 adds
explicit rejection of nonzero inactive-momentum rates and an unbounded-ray
guard, preserving the same passing physical results. The failed optimized
packed path is not used for a production trajectory.

All build commands, native diagnostics, source snapshots and CPU/device
trace arrays are retained. Sources/logs are losslessly gzipped. The manifest
records compressed and original SHA256s and the binaries kept only in GPFS.
Native run directories are
`/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/fv_primitives_v1` through `v7`.

For the current source, use a NEW GPFS directory and set `TMPDIR` there.
The exact recorded compiler commands use GCC `-O3 -std=gnu11`, G++
`-O3 -std=gnu++17`, and NVCC `-O3 -std=c++17 -arch=sm_89`. The CUDA audit
writes its trace with exclusive creation, so run a repeat from another new
directory. No fast-math/debug compiler option is needed for the final pass.

These are finite primitive/operator tests: zero new-FV PDE steps, no new
mesh/weight solve or restart conversion. They do not qualify actual-stage
positivity/GCL, full Euler accuracy, long Kepler, or GIZMO <=2x cost.
