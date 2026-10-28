# Same-checkpoint Ada GPU tessellation comparison at t=6.501

The [hash-guarded receipt](KEPLER_Q2_TESS_ADA_MATCHED_T6501_2026-10-27.json) compares three isolated CUDA builds on the same idle RTX 5000 Ada GPU, each restarted from the **identical accepted 5% site-jitter `t=6.5` snapshot** and evolved one RK4 step to exact `t=6.501`. Both W0 and projected-weight modes use the same physical parameters and 16 weight iterations. All accepted endpoints pass their independent one-step gates: W0 **11/11**, weighted **16/16**. The original weighted short auditor had checked an obsolete pressure-force banner; the corrected auditor rechecked the unchanged v3/v4 snapshots, all four logged weight stages, energy, power areas and physical fields before recording 16/16. No simulation was rerun for that audit correction.

| Mode | CUDA build | Kernel A median, four calls | Full tessellation median | Ten physical fields vs v3 |
|---|---|---:|---:|---|
| W0 | v3, 5×5 seed | 86.9895 ms | 105.976 ms | reference |
| W0 | v4, 3×3 seed | 112.211 ms | 131.7645 ms | bitwise equal |
| W0 | v7, original clipper + fixed limiter mode | 86.6170 ms | 105.048 ms | bitwise equal |
| weighted | v3, 5×5 seed | 87.6845 ms | 107.2245 ms | reference |
| weighted | v4, 3×3 seed | 109.334 ms | 128.647 ms | bitwise equal |
| weighted | v7, original clipper + fixed limiter mode | 86.1410 ms | 105.340 ms | bitwise equal |

The 3×3 seed is **29.0% slower in W0 Kernel A** and **24.7% slower in weighted Kernel A** despite fewer reported stencil expansions. A smaller seed leaves a larger intermediate polygon, so it may need more clipping or candidate work before the adaptive search; that mechanism is an inference, not directly counted. The v7 specialization lowers sm_89 kernel stack from **3,248 to 2,496 bytes/thread** and registers from 128 to 127, but improves Kernel A only **0.4% W0** and **1.8% weighted** in this four-call sample. Stack allocation is therefore not proven to be the dominant time cost.

The separate v6 branch combined the original 5×5 seed with an in-place exact clipper. Its very first native geometry pass reported **151 invalid cells** in both W0 and weighted runs, and neither produced an accepted endpoint. The exploratory CPU polygon prototype had matched 50,000 random clipping cases with nondegenerate output cycles, but the actual native result overrides that limited check. v6 is rejected pending a specific geometry fix; it is not an acceleration candidate for ongoing Kepler runs.

Because the same-card test establishes no useful speed gain from v4, its still-pending H200 duplicate **437439** was canceled before allocation (`Elapsed=0`). Existing long jobs and checkpoints were not changed. The current v3/v4/v7 measurements cover only one step on Ada; they do not establish relative H200 throughput or performance at late distorted geometries. Further speed work should profile actual candidate/clip counts or try a certified cached-neighbor path before replacing a validated native binary.

An attempted one-launch Nsight Compute profile of the accepted v3 W0 replay on Ada could not collect hardware counters: `ncu` returned `ERR_NVGPUCTRPERM` (GPU performance counters are not available to this user). The target still reached exact `t=6.501`. Its independently rerun short audit passed **11/11**, and all ten physical endpoint arrays were bitwise equal to the unprofiled v3 replay. The profiling command's nonzero exit is therefore a profiler-permission result, not a physics failure. Do not infer kernel occupancy or memory bottlenecks from this attempt.
