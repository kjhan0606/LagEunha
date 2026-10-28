# RK4 timing gap: isolated native measurements

The accepted A10 bounded run through `t=.825` spent 671.391 seconds in the
solver for 402 steps. Existing `CPU step` labels sum to 404.7 seconds, but
their `update` interval **contains** the GPU tessellation calls (264.034
seconds). Those two figures overlap. The existing `final` timer also begins
before the last `postStage` call, which is already timed as `post`; the
printed `CPU step` total double-counts that call. Therefore neither
subtracting GPU tessellation from wall time nor adding the old labels gives
a trustworthy exclusive time budget.

The [opt-in patch](KEPLER_RK4_GAP_PROFILE_2026-10-27.patch) applies cleanly
to the preserved [near32 source](KEPLER_FIXEDPOLY_NEAR32_SOURCE_2026-10-27/README.md)
(`patch --dry-run -p1` passed). It times the entire RK4 call and partitions
untracked intervals into prelude, after update, after force, and after post.
It also starts `final` immediately after the last timed `postStage`, removing
the existing overlap. `GFS_RK4_GAP_PROFILE=1` emits exact six-decimal-second
`[RK4 gap]` lines; closure is `wall - tracked - gaps` and should be near
timer overhead. No force, geometry, or state arithmetic is changed.

An independent syntax-local source clone at
`/tmp/Kepler_sparse_fixedpoly_profile_20261027` was built successfully with
the [reproducible build script](build_kepler_fixedpoly_gap_profile_20261027.sh).
Its patched `Exam/exam.c` SHA-256 is
`8b82b4a21a1fe3c7cf9df71c4d2a98bcb381ec9e4002088f8e7106930c0cc95d`.
The accepted GPU object is reused byte-for-byte (SHA-256
`e0d30024a9afc5abab3965ba39d01d04ea5171b01d527b7f9b905c5afd2ccd6c`),
and the separate profiled executable SHA-256 is
`a85b8f1349702db52d9372b5a79beddd1d57ee03ebb2943d244306e4bcf8f4dc`.
The live Slurm 437317 job has its own immutable staged baseline binary and
continues unmodified.

A [paired A10 runner](kepler_q2_bounded_rk_gap_pair_t085_a10_20261027.slurm)
ran as job **437330** from the same accepted `t=.825` checkpoint and
completed `0:0` in 49 seconds. The [pair receipt](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/rk_gap_pair_t085_a10_v1/pair_audit.json)
confirms both native physics gates, all ten parsed physical fields bitwise
equal at `t=.85`, and all 13 profile rows closing within `1e-5` seconds.
Excluding the first warmup step, median exclusive intervals per RK4 step are:

| Interval | Median seconds | Share of RK4 wall |
|---|---:|---:|
| Entire RK4 call | 1.6363 | 100% |
| Named update, force, post, final | 1.0491 | 64.1% |
| Unnamed after-force work | 0.5191 | 31.7% |
| Prelude before first update | 0.0657 | 4.0% |
| Other unnamed intervals | 0.0020 | 0.1% |

The after-force interval includes stage diagnostics, gravity, energy-audit
calls and RK4 state arithmetic. Source inspection finds `kepler_accel`
re-reading and parsing three unchanged environment values for every site
and every force/energy-audit stage. The isolated
[configuration-cache patch](KEPLER_ACCEL_CONFIG_CACHE_2026-10-27.patch)
changes only that repeated parameter lookup. A second profiled binary
(SHA-256 `58fa73d1b7f2ad36b478710b8ec4d0f2551300fe0e06ff530edb5c9ce664171b`)
was built from the same GPU object and submitted in the matched
[cache pair runner](kepler_q2_bounded_accel_cache_pair_t085_a10_20261027.slurm)
as A10 **437331**. It completed `0:0` in 48 seconds. Its
[pair receipt](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/rk_accel_cache_pair_t085_a10_v1/pair_audit.json)
reports all ten parsed physical fields bitwise equal at `t=.85`. Reading
both logs from **this same sequential pair**, the post-warmup median RK4
wall fell from `1.68204` to `1.41845` seconds (15.7%), and the after-force
interval from `0.550214` to `0.32048` seconds. The control ran first, so
some hardware warmup or node-load effect is possible; this is a short paired
result, not a sustained production speedup claim.

An independent later-state profile **437332** restarted the accepted
`t=1.651374332206546` bounded state and evolved both the unchanged and
timing-only binaries through 14 matched steps to `t=1.6785`. The Slurm
batch reported `FAILED 1:0` solely because its first auditor required the
entire serialized snapshots to have identical SHA-256 hashes. A read-only
[corrected audit](KEPLER_Q2_BOUNDED_RK4_GAP_LATE_RESTART_2026-10-27.json)
proved that the native headers, footers, and every byte after offset 8 in
all 55,351 particle records agree exactly; only the serialized linked-list
`next` pointer at record offset `0:8` differs across the two processes.
All 16 parsed snapshot fields and 14 step counts agree. The timing bins
close within `7.524e-6` seconds; median RK4 wall is `1.545312` seconds,
with `0.588602` seconds in the after-force interval. No simulation was
restarted to repair this packaging/auditor failure.

The gap accounting assumes `GFS_RK4_TRANSACTION` is off, as in the current
bounded run; rejected trial retries would need a separate interpretation.

## Why local stepping matters more than a larger global orbital parameter

The accepted bounded `t=0` native snapshot has 55,351 generators. With the
actual `EUNHA_KEPLER_DTETA=.06`, `eps=.01175` and the source formula
`dt_i=eta (r_i²+eps²)^(3/4)`, its smallest site radius `r=.10080469` gives
`dt_min=.00193984950`, matching the accepted native step `.00193984`.
This limiting generator is in the 1,476-site interior core (`r<.5`), which
contains only about 2.4% of the total initial mass. The 49,038 original
disk sites have minimum orbital limit `.0217480` (11.2 times the global
step); the 4,837 outer sites have minimum `.171091`.

The live GPU force kernel separately prints a face/CFL candidate near
`.004028` at `t≈3.5`, so increasing the **global** orbital coefficient
alone can at most roughly double this step before another measured limiter
acts. As an orbital-limit-only arithmetic upper bound, updating each site
at its own `dt_i` would reduce the nominal particle-update count by a
factor `N/dt_min / sum_i(1/dt_i)=38.3`. That number excludes hydro CFL,
neighbor synchronization, geometry rebuilds and conservative pair-flux
bookkeeping; it is **not** a predicted native speedup. It identifies the
small core as a measured reason to investigate conservative individual
time stepping after the bounded long-time physics gate, rather than relying
on one larger global `DTETA` or the 15.7% short cache optimization to reach
the paper's longest time.

The stronger [native local-step receipt](KEPLER_Q2_BOUNDED_LOCAL_DT_BUDGET_2026-10-27.json)
reads each cell's *actual serialized K4 hydro CFL* at record byte offset
128 and combines it with the orbital bound at both accepted `t=.001` and
`t=1.651374`. At both times only **eight cells** are within twice the
global step. For a hypothetical frozen `6.5`-time-unit span, power-of-two
block scheduling changes nominal active-cell updates from `185,481,201`
to `6,943,434`, an **arithmetic ratio 26.71**. The disk's minimum combined
local bound is `.021748`, and its minimum hydro bound is about `.0974`.
This still excludes asynchronous face fluxes, mesh prediction and GPU
work, so it cannot be used as a wall-time forecast. It corroborates the
earlier [full q1/q2 native timestep audit](KEPLER_TIME_STEP_SCALING_2026-10-26.md)
on the present sparse bounded IC and confirms that the local-step range
persists beyond initialization.

The [frozen native neighbor-graph audit](KEPLER_Q2_BOUNDED_FACE_BUDGET_2026-10-27.json)
reconstructs 166,014 undirected Delaunay edges at accepted `t=1.651374`;
twice that count is close to the native GPU's roughly 333,000 directed
face rows. Under the same power-of-two schedule, a face evaluated when
either endpoint is active has a nominal **26.24-fold** smaller update
count than evaluating every face at the global minimum step. Enlarging
each active set to its neighbors by one, two or three graph rings leaves
site-update arithmetic ratios of 25.32, 23.73 and 22.09 respectively.
These are topology-based counts at one instant. A correct local CUDA mesh
builder, asynchronous conservative flux register and wall interface are
still required before any wall-time or physical accuracy claim.

The remaining `after_force` cost is being resolved with a
[diagnostic-only stage timer patch](KEPLER_RK4_STAGE_DIAG_TIMER_2026-10-27.patch).
It was applied on top of the acceleration-cache and gap-timer branch and
built as SHA-256
`30d81c891dfaa394d67f37ea6069eefd46ff93eddd4f46286e2056d7c5cd3626`,
again with the exact accepted GPU object. The [paired short runner](kepler_q2_bounded_stage_diag_pair_t085_a10_20261027.slurm)
ran as A10 **437335** and completed `0:0`. Its
[pair receipt](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/rk_stage_diag_pair_t085_a10_v1/pair_audit.json)
passes both native gates with all parsed physical fields bitwise equal.
The ten full-state diagnostic calls take **0.274181 seconds median per
step**, within a `0.329734`-second after-force interval, or about 19% of
the `1.442126`-second RK4 wall time.

The diagnostic function calls `gfs_thermal64_requested()` through
`gfs_t64_U` and `gfs_t64_V` for each of 55,351 cells on each of ten stages.
The flag is fixed for a run. The
[local hoist patch](KEPLER_RK4_DIAG_THERMAL_HOIST_2026-10-27.patch) reads
it once per diagnostic call while retaining all ten full-state checks.
A new binary SHA-256
`d92b6cd3b8491678cdfafa7610608378db89ca9bd16585a52a51aaf26f3b3367`
was built with the same GPU object. A10 **437336** completed `0:0` in 43
seconds. Its [pair receipt](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/rk_diag_hoist_pair_t085_a10_v1/pair_audit.json)
passes both physical gates and all ten parsed physical-field comparisons
bitwise. Within the same A10 allocation, post-warmup median whole-RK4 wall
fell from `1.467731` to `1.1830635` seconds (19.4%). The ten full-state
diagnostic calls fell from `0.2787475` to `0.0075715` seconds. The optimized
[source overlay and exact patch chain](KEPLER_BOUNDED_FAST_RK4_2026-10-27/README.md)
have been preserved; neither active long-run binary was changed.
