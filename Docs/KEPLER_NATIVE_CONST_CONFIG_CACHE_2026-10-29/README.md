# Native constant-configuration cache candidate

The measured native RK4 timing showed repeated host configuration work outside the reported force/update intervals. This opt-in candidate caches the process-constant thermal-mode selection and the central-potential diagnostic parameters. It changes no force formula, integration step, face geometry, mass, thermal update or conservation-ledger placement.

Enable with `GFS_NATIVE_CONST_CONFIG=1` only when `GFS_THERMAL64`, `EUNHA_IC`, and the Kepler center/softening environment parameters remain fixed throughout a process. Default-off retains the original environment reads for the thermal mode and potential parameters. Atomic publication covers concurrent host initialization. The existing acceleration configuration cache is unchanged. Time-dependent Laguerre weights are not cached by this feature.

`BUILD.json` pins the two base files, zero-context reconstruction patch and frozen candidate. Reconstruct native wall-events v2 first, verify the base hashes, then apply `CONST_CONFIG_CACHE.patch` with `patch -p1` or `git apply --unidiff-zero`. Exact reconstruction against the tested source was verified. The native C and CUDA build succeeded on Ada; `build_ada.sh` records the modules and commands. Large binaries remain on GPFS.

`ABBA_AUDIT.json` records actual native GPU runs from the unchanged accepted Q1 t=.825 checkpoint to exact t=.82845. The order was original/cached/cached/original, followed by candidate-flag-off. All complete restart record bytes except the runtime next pointer, all 52 footer bytes and all printed numeric ledger lines matched. Source/input hashes, fixed masses, W0, finite/positive state and zero injections passed. No full-particle TSV was emitted.

After excluding the first two startup steps, ABBA median native RK4 wall time decreased from 0.758701 to 0.500953 seconds, a measured **1.5145x** speedup. Full MPI execution times, including initialization, were 7.185/7.301 seconds for the originals and 5.496/5.553 seconds for the cached runs. These are short same-geometry Ada measurements; they do not establish H200 speedup, MPI scaling, thermal64-enabled performance or cost at late distorted geometry.

`LATE_COLLISION_AUDIT.json` checks the gapped-ID accepted Q2 predecessor through exact t=48.4 with cache enabled. It passes the existing physical audit, one real collision and 13 rollbacks. Its complete stored state and footer match frozen v2 exactly. This guards the configuration cache across actual rejected attempts; it is not a new long-time disk result.

The frozen candidate executable is `/gpfs/kjhan/LagForce/Kepler_SlipWall_Integration_20261029/candidate_const_cache_v1/eunha2.exe`, SHA256 `7364a05c09eeb393d4e8528f01a82efb7f6d9c6fba9b394d2d1e62e6bd78345f`. It remains separate from H200 job 439333, which continues the originally qualified v2 checkpoint chain. The candidate is ready for a future accepted-checkpoint continuation; no active executable was replaced.
