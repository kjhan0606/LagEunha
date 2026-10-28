# Native polar Kepler time-step bottleneck (2026-10-26)

## Active source provenance (checked 2026-10-27)

The q1/q2 running binary is built from the isolated
`/tmp/Kepler_weighted_native_20261026` source, itself based on
`/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src`.
The [build receipt](KEPLER_POLAR_NATIVE_BUILD_2026-10-26.json) records the
modified `Exam/KH/kh.c` and `util.c` hashes and binary hash. The similarly
named files under this project's `src/Exam/KH/` have **different hashes**
from both the base and the active isolated source. Performance or physics
changes intended for these H200 comparisons must therefore be derived from
the build-receipt source and verified against the running binary's hash;
editing the project `src` alone would not affect the queued jobs. The
existing source patches in `Docs/` preserve the experimental changes.
Both current q1 and q2 continuation executables were independently hashed
and equal the receipt's SHA-256
`7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284`.

## Measured and computed inputs

The active 384² q1 continuation runs one MPI rank, two OpenMP threads and one
H200, using the native RK4/CUDA executable. Its log near t=0.064 shows four
GPU tessellations of about 26 ms each, a force-call `SoA` interval of about
35 ms, and a printed CPU-step breakdown of about 0.6–0.7 s. The `SoA` label
includes host periodic-padding refresh and GPU re-upload; it is not merely
array packing. At the same point the native global step is 7.64202e-5.
These are log intervals, not an end-to-end profiler or a sustained throughput
measurement.

`Exam/exam.c::kepler_dt_limit` takes the **minimum over every cell** of
`eta (r²+eps²)^(3/4)`. At eta=0.06 and eps=0.01175, the q1 and q2 polar
layouts each have one central site at r=0 that fixes this minimum. The q2c
layout starts with six sites at r=0.01104854, giving 1.2289879e-4. A full
step updates and tessellates all 147,456 cells, regardless of the cell that
set the minimum.

The [initial-layout calculation](KEPLER_ORBITAL_TIMESTEP_BUDGET_2026-10-26.json)
uses the exact native glass files and this C expression. For t=0→0.825:

| Layout | Global orbital dt | Sites within 2× minimum | First disk site orbital dt | Global steps |
|---|---:|---:|---:|---:|
| q1 | 7.6420e-5 | 1 | 0.02122 | 10,796 |
| q2 | 7.6420e-5 | 1 | 0.02175 | 10,796 |
| q2c | 1.22899e-4 | 6 | 0.02246 | 6,713 |

The q1 and q2 preflight logs report first-step **global hydro CFL candidates**
of 6.824e-4 and 2.026e-3 respectively. Therefore the explicit Kepler
orbital bound, rather than those hydro candidates, sets the actual initial
global dt. If every cell's local step were capped by the *observed global*
hydro candidate and the local orbital bound, a power-of-two active-cell count
would be about 8.0 times lower for q1 and 16.0 times lower for q2 than the
current all-cell work count. This deliberately conservative cap still assumes
local updates, so it is a **cell-update arithmetic ratio, not an achievable
wall-time speedup**. q2c has no native first-step hydro measurement yet.

The full t=0.001 native snapshots serialize each cell's final K4 hydro CFL
in the `dt` field at offset 128 of the 576-byte record. Its minimum agrees
with the final `GPU step 14` log: 6.7979e-4 for q1 and 2.2671e-3 for q2.
All 147,456 values are finite and positive. In the disk, the minimum local
bound `min(hydro CFL, orbital dt)` is about 0.02122 for q1 and 0.02175 for
q2, versus 7.6420e-5 globally. Only **one cell** lies within twice this
global bound in each layout. If these exact K4 per-cell bounds were frozen
for t=0→0.825, power-of-two scheduling would reduce the number of active
cell updates by factors of 832 and 809 respectively. The
[native snapshot calculation](KEPLER_NATIVE_LOCAL_DT_BUDGET_2026-10-26.json)
and [reproducer](kepler_native_local_dt_budget_20261026.py) record every block
level. These much larger ratios describe only the available local timestep
range at one instant. They exclude moving-face flux work on inactive cells,
re-tessellation, GPU launch/transfer overhead, positivity, and how dt changes
later; they must **not** be reported as forecast speedups.

The script [kepler_orbital_timestep_budget_20261026.py](kepler_orbital_timestep_budget_20261026.py)
reproduces those counts. Ratios without the hydro cap are much larger and
should not be used to forecast run time. The first-step hydro candidate need
not stay fixed later. Moving faces couple active and inactive cells, and a
full-mesh tessellation at every smallest block would erase most of the
arithmetic gain.

## Why the orbital cap cannot simply be removed

An independent [softened central-force RK4 control](KEPLER_SOFTENED_ORBIT_RK4_LIMIT_2026-10-26.json)
starts the innermost nonzero q1/q2 ring site on its exact circular orbit and
integrates the **same Plummer acceleration** as the native Kepler path, with
no mesh or hydro force. At t=0.825, position error divided by the initial
radius is:

| Inner ring | Current global orbital step | Whole domain at observed hydro CFL | Hydro macro step with gravity subcycled to current global orbital bound |
|---|---:|---:|---:|
| q1, r=0.015625 | 1.08e-5 | 0.550 | 1.02e-5 |
| q2, r=0.022097 | 1.42e-6 | 1.000 | 1.34e-6 |

The q2 large-step trajectory has lost essentially all angular momentum and
collapsed toward the softened center in this one-body control. At t=6 the
q1 large-step error is 0.736 initial radii, while the q2 trajectory remains
collapsed. Gravity-only subcycling to each ring's *local* orbital limit does
better than no subcycling but has 0.94% and 1.40% position errors at t=6 for
q1 and q2; the current global-bound subcycling recovers the baseline orbit
error (about 5.0e-4 and 5.9e-5). The disk inner edge at r=0.5 stays accurate
at the observed hydro CFL in this pure-gravity control. The
[reproducer](kepler_softened_orbit_rk4_limit_20261026.py) records all cases.

This establishes a numerical constraint: increasing the whole-domain RK4
step to the present hydro CFL without treating the inner orbits differently
is invalid. The subcycle column shows only that inexpensive local **gravity
arithmetic** can retain the old one-body accuracy while the outer macro step
is larger. It does not demonstrate a valid gas integrator. During each macro
step moving generators and faces must be predicted to the same time, and
their momentum/energy exchange must remain conservative and positive.

## Consequences for implementation

1. Preserve the q1/q2/q2c RK4 comparisons as physical controls. Changing the
   limiter alone would under-integrate the inner atmosphere and weaken their
   validity. The centerless q2c experiment is a controlled small improvement,
   not a route by itself to the paper's long endpoint.
2. Profile the native whole-step cost and validate the already implemented
   one-rank GPU local-periodic refresh on a polar checkpoint before using it
   in the long run. Earlier C010 A/B was bitwise equal through t=2.01, with
   the force-stage interval falling from roughly 35 ms to 16 ms, but it has
   not been checked on the new polar layout.
3. A block-step KDK/DKD implementation must drift all face generators to the
   same evaluation time, update each active/inactive face flux once with equal
   and opposite mass/momentum/energy increments, and retain pending fluxes
   until the neighbor's kick. The discrete face speed, gravity work, and
   snapshot times need the same time-centering. Check fixed cell mass and
   W=0 in the initial mode; introduce time-varying weights and any mass
   exchange only with their own conservation audit.
4. The GPU path must restrict tessellation/face work to active neighborhoods
   and measure the cost of their halos. A scheduler layered on top of four
   full-domain RK4 tessellations per smallest block cannot realize the
   active-cell work reduction. Preserve one-rank and multi-rank correctness
   separately; one-rank `migrateTreeVorork4Particles` already skips actual
   `pmigrate` when `NDDINFO=0`, so removing that call is not a useful first
   optimization.

These findings identify a large time-integration scaling problem independent
of whether polar alignment resolves the disk's shear-induced density error.
The native q1 `t=.825` density/energy audit has now passed all 15 current
checks; see the [control report](KEPLER_POLAR_NATIVE_CONTROL_2026-10-26.md).

## Independent polar local-sync A/B queued

Job 436902 was submitted on one H200 in
`/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1_local_sync_t001`.
It uses byte-verified copies of q1's native binary, glass, and parameters;
only `GFS_GPU_LOCAL_PERIODIC_SYNC` changes from 0 to 1. The original q1
preflight produced exact t=0 and t=0.001 snapshots and passed its physics
audit. The new job will run the same 14-step endpoint in a separate directory.
The [A/B audit](audit_kepler_polar_pairwise_20261026.py) will compare every
saved physical field by bits after sorting by cell ID, confirm the device
alias refresh actually activated, and compare median GPU-wrapper intervals.
At the last scheduler check job 436902 was pending GPU quota; the q1 long
run 436889 remained running near t=0.085 with no nonpositive internal energy
or floor injection. No A/B outcome is claimed yet.

**Completed A/B result, 2026-10-27:** Job **436902** completed the same
14-step q1 `t=.001` preflight using GPU-local periodic image refresh. The
[ID-wise pairwise audit](KEPLER_POLAR_Q1_LOCAL_SYNC_PAIRWISE_2026-10-27.json)
finds exact matching snapshot times and **bitwise equality** in every saved
`x,y,vx,vy,mass,volume,ie,den,pressure,w2` field at both `t=0` and
`t=.001`. The log confirms local alias refresh was used and host halo
refresh was not. Across 14 matched GPU-wrapper samples, the reported
`SoA` interval median fell from `36.25` to `0.1` ms and wrapper total
from `71.55` to `37.35` ms per RK stage. The `SoA` interval includes host
periodic padding refresh and particle packing, so the 362.5-fold ratio of
that interval is **not** a whole-step or wall-speedup. The wrapper reduction
is about 47.8%; a longer same-IC, same-time run is still needed before
promoting this option in the native Kepler chain. Existing queued physics
comparisons keep `GFS_GPU_LOCAL_PERIODIC_SYNC=0` for matched baselines.

An evolved-state confirmation is now staged as two independent low-priority
H200:1 jobs: **437180** (host refresh) and **437181** (device-local refresh).
Both use the same audited q1 accepted step-8119 checkpoint at
`t=0.6187644470884642`, the same halo-12 native binary and physical
settings, and stop at exact `t=.65`. The only runtime difference is
`GFS_GPU_LOCAL_PERIODIC_SYNC=0/1`; the
[restart selector](prepare_kepler_q1_local_sync_probe_20261027.py) checks the
real q1 halo-failure signature, accepted audit, binary SHA and input state
before copying. The [host](kepler_q1_halo12_local_sync_host_to065_20261027.slurm)
and [device](kepler_q1_halo12_local_sync_device_to065_20261027.slurm)
runners byte-check their staged scripts and run the same full physics and
radial audits. They will test alias-refresh fidelity after the outer seam
has already deformed, and allow whole-step timing comparison over hundreds
of steps. They have been submitted, not completed; do not extrapolate the
14-step result to `t=600` runtime.

## Whole-step wall timing gate

At 31m36s of job 436889's elapsed time, the native log had reached
approximately t=0.11647 after 1,511 continuation steps. This is about
1.25 wall seconds per early step, while the existing `[CPU step]` line
reports roughly 0.6–0.7 s in its named subregions. The discrepancy has not
been attributed; startup, untimed RK4 loops, logging, and filesystem costs
are all candidates. Extrapolating the present early pace to t=600 gives an
order-100-day single-H200 cost, which is **only a scale warning**, not a
measured long-time forecast. Block stepping and GPU active-neighborhood
work would need a validated implementation before that endpoint is practical.

An isolated copy of the user's native source now has optional
`GFS_KH_WALL_PROFILE=1` timing around the full RK4 call and the entire KH
loop iteration. The [minimal source patch](KEPLER_KH_WALL_PROFILE_2026-10-26.patch)
and [binary receipt](KEPLER_KH_WALL_PROFILE_BUILD_2026-10-26.json) show that
only `Exam/KH/kh.c` was changed; all active polar binaries remain intact.
H200 job 436905 was submitted with `afterany:436902` for the same q1 14-step
preflight in a separate directory. Its audit must establish the same
initial/endpoint physical fields before using its timing to decide whether
the missing interval is inside RK4 or in output/driver work. No profiler
result is yet available.

**Completed wall profile, 2026-10-27:** Job **436905** completed and the
[pairwise audit](KEPLER_POLAR_Q1_WALL_PROFILE_PAIRWISE_2026-10-27.json)
found every saved field bitwise identical to the original q1 preflight at
both snapshots. All 14 iterations have `[KH_WALL]` coverage. The median
whole KH iteration is `1.3202475` s, with only `1.45e-5` s outside its
evolution call; the median evolution call is `1.320233` s. The separately
printed `[CPU step]` named-subregion total has a coarsely rounded median
`0.7` s. Thus the previously missing roughly `0.6` s per step occurs
**inside** the RK4 evolution function, not in KH output scheduling or the
driver. The named timer excludes RK stage arithmetic, stage diagnostics,
some state loops and bookkeeping; this result does not yet allocate that
remainder among them. The GPU-wrapper median is `74.6` ms per stage on
the profiled run. A finer RK-internal timer is needed before optimizing a
specific loop. Disabling physical checks solely for apparent speed is not
justified by this measurement.

**Evolved-state wall sample, 2026-10-27:** The q1 halo-12 recovery completed
2,737 native steps from the accepted step-8119 `t=.618764447` checkpoint
to exact `t=.825`. Its measured simulation wall time was 4,502.544 s,
or **1.645 s per step** over this deformed outer-state interval. A purely
constant-step extrapolation at `dt=7.64202e-5` would cost about **149.5
days on one H200** to reach `t=600`. This is a cost scale, not a prediction:
future dt and per-step cost can change, and the sparse/core and active-cell
approaches may alter both. It confirms that merely extending this full-grid
baseline cannot meet the final benchmark economically. The q1 job's Slurm
exit was `FAILED 1:0` solely because the final radial auditor missed a
Python import in its staging directory; the native endpoint and subsequent
strict physics and rescued radial audits passed without replaying the run.

**What the queued local-sync A/B can buy:** In the live q2 halo-12 recovery
near `t=.755`, the current host-refresh `[GPU step] SoA` interval is about
`48 ms` per RK4 stage; four stages therefore expose about `0.19 s` of this
particular interval per accepted step. This interval includes host periodic
padding rebuild, checks, SoA fill, and upload, as documented in the
[local-fallback design](KEPLER_LOCAL_PERIODIC_FALLBACK_DESIGN_2026-10-27.md).
The q1 evolved halo-12 run took `1.645 s` per full step. Even eliminating
the q2-sized interval completely would address only order 10% of a step at
that measured wall scale, not the roughly 149-day full-grid `t=600` cost.
The evolved-state A/B jobs 437180/437181 remain useful for fidelity and
modest throughput, but their low scheduling priority should not delay the
corrected sparse IC's much larger predicted timestep gain. These are
cross-run scale estimates; the A/B's own matched whole-step timings are
required before claiming a speedup.

**RK4 internal profile queued, 2026-10-27:** An isolated source copy adds
wall timers for the ten RK4 stage-diagnostic calls, ten direct transaction
validity calls, and the full RK4 function. The existing named CPU-section
times are retained; their difference from the full wall after the two new
measurements is reported as a residual, not assigned speculatively to one
algorithm. The [source patch](KEPLER_RK_FINE_PROFILE_2026-10-27.patch),
[build receipt](KEPLER_RK_FINE_PROFILE_BUILD_2026-10-27.json),
[H200 runner](kepler_q1_rk_fine_profile_t001_20261027.slurm), and
[timing gate](audit_kepler_rk_fine_timing_20261027.py) record the opt-in
experiment. Job **437200** repeats the same 14-step q1 `t=.001` native IC
as the completed wall-profile control and has low scheduling priority. Its
runner requires native IC gates, bitwise equality of every saved physical
field with the earlier control, and 14 matched RK4/KH timing rows before
any performance inference. No timing result exists yet; the active q1/q2
physics binaries are unchanged.

**Halo-28 continuation capacity review, 2026-10-27:** The queued q1/q2
`.825→2.221` runners request 24 hours, and the `2.221→6.5` runners request
48 hours. The active source checks `n_real+n_pad` before each GPU use and
reinitializes the particle/face context when it exceeds capacity;
`gpu_build_cell_csr` likewise reallocates for the larger linked-cell grid and
entry count. Therefore the `max_particles=274722` printed by a halo-12 run is
not a fixed halo-28 ceiling. This is a source-level capacity check, not a
successful native halo-28 run. Halo-28 padding may substantially raise the
per-step cost, so the first actual continuation must supply a measured
throughput before the 24/48-hour limits can be judged safe. Preserve its
periodic checkpoints if the later 48-hour run approaches its limit; a timeout
cannot be interpreted as a physical failure or justify restarting from zero.

## Low-mass core oversampling and sparse-glass control

**Superseding IC gate, 2026-10-27:** The sparse preflight checks below missed
global periodic area conservation. Their `t=0` areas sum to about 30,470
instead of the box area 36, and their one-time masses are wrong even though
their disk annulus is initially unchanged. Their queued `t=.825`/`t=10`
continuations were canceled before execution. The
[mass-error diagnosis and corrected native preflights](KEPLER_SPARSE_IC_MASS_CORRECTION_2026-10-27.md)
supersede any earlier statement here that the sparse IC was physically
validated. No sparse long-time conclusion should use the old checkpoints.

The q=2 native `t=0` snapshot has 129 cells at `r<0.1`, although their
combined mass equals only **1.66 median disk-cell masses**. All 147,456 cells
in the 6×6 box carry the mass of about 53,959 such disk cells. Uniform
geometric sampling therefore spends most cells in low-density gas, and one
site at `r=0` fixes the global orbital step. These are measured native
masses, not an assumed equal-mass IC.

An isolated native-source [patch](KEPLER_SPARSE_GLASS_NATIVE_PATCH_2026-10-26.patch)
adds opt-in `EUNHA_ALLOW_SPARSE_GLASS=1`: a glass file may contain fewer
sites than `Nx*Ny` while `Nx=Ny=384` remains the nominal disk resolution.
The usual full-count behavior is unchanged, and an invalid sparse glass
fails rather than silently falling back to a Cartesian grid. The
[build receipt](KEPLER_SPARSE_GLASS_NATIVE_BUILD_2026-10-26.json) records
the isolated binary and hashes.

The first [generator control](make_kepler_sparse_core_glass_20261026.py)
keeps every q=2 disk and outer site, removes 123 sites inside `r<0.1`, and
retains six at `r=0.08838835`. It remains a gas-filled Voronoi domain;
this is neither a vacuum excision nor an equal-mass IC. The
[glass manifest](KEPLER_Q2_SPARSE_CORE_GLASS_2026-10-26.json) records 147,333
output sites. An independent initial tessellation finds all 49,038 disk
sites identical and maximum relative disk-cell-area difference `8.88e-11`
from original q=2; see the
[geometry receipt](KEPLER_Q2_SPARSE_CORE_INITIAL_GEOMETRY_2026-10-26.json).
An independent area-times-target-density estimate gives total mass inside
`r<0.5` only `0.0375%` higher than original q=2, although the tiny
`r<0.1` region itself is `9.66%` lower; see the
[core-mass geometry receipt](KEPLER_Q2_SPARSE_CORE_MASS_GEOMETRY_2026-10-26.json).
The minimum orbital bound would rise from `7.642e-5` to `1.598e-3`, a
**20.9-fold bound change**, but the actual global step may be smaller due
to hydro CFL or core geometry. H200 preflight job **436938** is queued to
measure the step and test fixed mass, W=0, positive state, and orbit
accuracy using the [sparse auditor](audit_kepler_sparse_core_preflight_20261026.py).
Do not promote this variant to a long run without that native gate. Even
a 20-fold step improvement does not by itself make `t=600` practical;
local time stepping and further GPU work reduction remain relevant.

A separate full-count A/B job **436939** uses the same patched binary with
the original q=2 glass and ordinary (unset) sparse flag. Its completed
`t=0` and `t=0.001` snapshots will be compared by ID and field bits to the
original q=2 preflight. A passing sparse-core physics audit alone cannot
prove the source patch preserved the default path; this independent control
closes that gap.

A static sparse-input path review found that `Exam/KH/util.c` allocates the
old `nx*ny` capacity but iterates only the glass header's actual site count,
sets contiguous IDs from those sites, and assigns the resulting `np` to
`VORO_NP` before migration, tessellation, RK4, and CUDA calls. Searches of
the active `Exam/exam.c` and `Exam/exam_gpu.cu` paths found no `nx*ny`
evolution loop; their cell work uses `VORO_NP`/`nbp`. The staged 54,312-site
glass has every coordinate strictly inside `[0,6)^2`. Standalone legacy
drawing utilities do contain `Nx*Ny` loops, but are separate executables
and are not used in these runs. This is a code-path review, not a substitute
for the pending native preflight and full-count A/B test.

## Profile-sampled 54k-site candidate

The [profile-sampled glass generator](make_kepler_profile_sampled_glass_20261026.py)
uses a target site mass equal to the measured q=2 median disk-cell mass.
It retains every one of the 49,038 q=2 disk generators, grades polar rings
in the low-density core and outside the disk, and uses a coarse periodic
Cartesian background. An initial version without transition buffers had a
**22.4% maximum disk-cell area change** at the inner edge and is rejected
as a matched disk IC; its [geometry receipt](KEPLER_Q2_PROFILE_SAMPLED_GEOMETRY_2026-10-26.json)
is retained to show that failure.

The improved [buffered candidate](KEPLER_Q2_PROFILE_BUFFERED_GLASS_2026-10-26.json)
keeps original q=2 sites in `0.48<=r<0.5` and `2<=r<2.02` as guards.
It has **54,312 sites** versus 147,456 in the reference and minimum radius
`0.10080469`. Independent [Voronoi geometry](KEPLER_Q2_PROFILE_BUFFERED_GEOMETRY_2026-10-26.json)
finds all disk generators identical, maximum relative initial disk-cell
area difference `1.18e-10`, and exact-orbit mass-weighted disk density
change `1.9699e-4` at `t=600`, matching the original q=2 geometry control.
An [independent q=1/q=2 long-geometry comparison](KEPLER_POLAR_Q1_Q2_LONG_GEOMETRY_2026-10-26.json)
finds q=2's `t=600` disk area-induced density error is `0.326` times q=1's
and remains below `0.374` times q=1 through a nominal 1000 inner-edge
periods. This favors q=2 as the geometric control, while native heating and
runtime remain unmeasured at those times.
The present orbital bound would increase **25.38-fold** to `0.001940` and
the site count falls **2.715-fold**. Their product is an idealized 68.9-fold
reduction in cell updates, **not a measured wall-speedup**; the hydro CFL,
CUDA stencil behavior, wall overhead, and native core dynamics remain
unknown. This is a profile-guided placement, not a proven exactly equal-mass
glass.

H200 native preflight **436940** uses the isolated sparse-glass binary,
exact `t=0` and `t=0.001` snapshots, and the
[buffered candidate auditor](audit_kepler_profile_buffered_preflight_20261026.py).
It checks site count, minimum radius, disk mass scatter, target initial
density, fixed mass, W=0, positive state, and disk orbital accuracy before
any long continuation is promoted.

The original narrow 8-cell halo continuations **436941** (to `t=.825`)
and **436942** (to `t=10`) were canceled before execution after the
2026-10-27 periodic-halo screen below. Their
[restart selector](prepare_kepler_profile_buffered_restart_20261026.py)
checks the audited checkpoint header, size, hash, and site count before use.
The [t=0.825 auditor](audit_kepler_profile_buffered_t0825_20261026.py)
requires positive state, exact fixed mass, W=0, disk density change below
0.005, mass-weighted orbit displacement below `1e-4`, and disk radial-speed
RMS below 0.01 of circular speed. The originally staged `t=10` runner
would have advanced only a passing `t=.825` checkpoint. Its
[restart selector](prepare_kepler_profile_buffered_to10_20261026.py) and
[physics auditor](audit_kepler_profile_buffered_t10_20261026.py) require
complete checkpoints and test GIZMO-annulus L1 velocity/density, disk
entropy, angular momentum, total energy, and zero floor injection.
These old jobs did not run; their 12→32-cell replacements are described
below. The actual t>6 and t=10 native outcomes remain unknown.

## Wider transition-buffer backup (geometry only)

The [full periodic Voronoi audit](KEPLER_Q2_PROFILE_BUFFERED_PERIODIC_GEOMETRY_2026-10-26.json)
of the staged 54,312-site input finds positive cells, total area
`35.9999999999984` in the box of area 36, and maximum cell circumradius
`0.1640`, well below the nominal eight linked-list-cell halo scale `0.5`.
It also reveals that the immediately adjacent **non-disk** guard rings
change area relative to full q=2: the inner `0.48<=r<0.5` ring's total
initial mass at prescribed density is `0.05009` versus `0.06271` in the
reference, and the outer `2<=r<2.02` ring's is `0.22388` versus `0.20785`.
The disk cells themselves remain unchanged, so this is a potential
transition-source concern for later gas evolution, not a demonstrated
disk failure.

A [wider-buffer candidate](KEPLER_Q2_PROFILE_WIDE_GLASS_2026-10-26.json)
keeps original q=2 sites over `0.44<=r<0.5` and `2<=r<2.06`, using
**55,931 sites** total (only 1,619 more). Its durable glass is
[here](KEPLER_Q2_PROFILE_WIDE_GLASS_2026-10-26.bin), SHA-256
`de50b97639a243f9384fc7b4bfc402f9eccb52386f30ecdd3450f7b0acaa9e8f`.
The [independent full-periodic audit](KEPLER_Q2_PROFILE_WIDE_PERIODIC_GEOMETRY_2026-10-26.json)
again partitions area 36 exactly to numerical precision and lowers the
maximum circumradius to `0.15525`. It restores the immediate inner and
outer guard-ring masses to `0.062709890657` and `0.207846299504`, equal
to full q=2 at the quoted precision. Its 49,038 disk generators are
bitwise identical to q=2, initial disk-cell areas differ by at most
`8.88e-11`, and the [exact-orbit t=600 density error](KEPLER_Q2_PROFILE_WIDE_DISK_GEOMETRY_2026-10-26.json)
is `1.9698956e-4`, the same as the staged candidate. The minimum orbital
bound remains approximately 25.38 times baseline; the count reduction is
2.636 times, for an idealized 66.9-fold cell-update reduction rather than
68.9-fold. Neither is a measured wall-speedup. This input is retained as a
backup and has **not** replaced or restarted pending native job 436940.
Native tests must decide whether the added transition rings matter.

## Focused Astra IC review (2026-10-26)

At the user's request, a single independent Astra review checked these
profile-sampled ICs against [Hopkins (2015), section 4.2.4](https://arxiv.org/html/1409.7395#S4.SS2.SSS4).
The paper starts with an effective uniform 256² particle placement and
removes particles inside `r<0.5` and outside `r>2` according to the target
density profile. Sparse sampling of low-density gas is therefore a
defensible benchmark design. The polar q=2 alignment, larger 6×6 periodic
box, and exact site layout remain different from the paper's IC; their
advantages must be labeled as IC effects rather than solver improvements.

The review distinguishes **fixed mass** from **equal mass**: the one-time
assignment `m_i=rho_target(r_i)*V_i(0)` followed by constant `m_i` is
mathematically consistent, though the core and atmosphere represent
different discrete mass and momentum distributions from full q=2. No site
inside `r≈0.1008` is a resolution choice, not a vacuum: Voronoi cells still
cover the origin. Their generator-centered orbit/gravity values need not
accurately integrate the steep softened potential over the whole cell.
The 54,312/55,931-site profile-sampled inputs actually have **eight** sites
on their first ring at `r≈0.10080469`; the separate sparse-core-only candidate
has six retained sites. The profile-sampled ring creates a symmetric
high-order Voronoi junction near the origin; the positive-area SciPy audit
does not prove native GPU topology stays quiet over many inner orbits.

The reviewer favors the **wide** 55,931-site buffer for a later native
comparison because it restores the immediate guard geometry at only about
3% more sites, while retaining the queued narrow run as a control. A
constant initial pressure means the narrow guard area differences are not
by themselves an initial pressure-force imbalance. The most consequential
long-time uncertainty is propagation of outer periodic-boundary artifacts
and core force/entropy errors into the disk. A [native seam audit](KEPLER_Q1_SEAM_T041254_2026-10-26.json)
at `t=.4125` places all 44 short pairs within `0.115` of a box edge, with
relative speed about unity at separations around `0.002`; near-edge
mass-weighted density and pressure changes are `0.0681` and `0.1298`.
This is a concrete outer-background signal to track, not evidence that it
has reached the disk or proof of a unique boundary cause.

The minimum discriminating native comparison is full q=2, narrow, and wide
with the same solver and physical time. Record ID-wise fixed mass and
`rho=m/V`, positivity and floor injection, radial bins of density,
pressure, entropy `P/rho^gamma`, radial velocity, angular momentum and
mass flux across `r=.5` and `r=2`, total energy, and wall/RHS cost.
Compare a common small timestep across ICs, then a half-step for a sparse
candidate to separate space and time errors. Evaluate the paper's annular
L1 metrics and density fields at native `t=120` and `t=600` before a
long-time stability claim. A solver-improvement claim additionally needs
old/new solver comparison on the **same IC**; IC acceleration alone does
not establish that claim.

## Initial sparse-core gravity quadrature (2026-10-26)

An [independent local-Voronoi quadrature](KEPLER_Q2_SPARSE_CORE_FORCE_2026-10-26.json)
measures the gravity source represented by the innermost **eight** cells.
Each narrow/wide generator at `r≈0.10080469` receives point acceleration
`≈96.44`; its cell's *area-averaged* inward Plummer acceleration is `≈190.72`,
a ratio of `1.9770–1.9779`. The corresponding angular sites in full q=2 lie
at `r≈0.11048543`; their ratio is only `1.0022–1.0044`. A sparse core cell
has median area `≈0.009057` versus `≈0.0002397` in the full q=2 comparison.
The prescribed target-density average over those sparse core polygons is
`≈1.1675` times the density evaluated at the generator. Two local Voronoi
neighborhood radii give polygon areas within `7.2e-15`; order-64 versus
order-128 quadrature differs in acceleration by at most `4.5e-12`.

This is a large **initial source-representation difference**, not a measured
native instability. The current code assigns every cell the circular speed
of its generator and applies gravity at that generator, so replacing only
the gravity source by a cell average would destroy this pointwise orbital
balance. The native preflight/continuations must measure core radial motion,
entropy and pressure, then whether a disturbance reaches the disk. If the
core limits stability, a consistent core velocity/source treatment or local
time stepping with restored core resolution is required; integrating gravity
over every cell at every stage is not implied by this audit.

An independent [local exact-orbit Voronoi control](KEPLER_Q2_SPARSE_CORE_BALLISTIC_2026-10-26.json)
separates this source issue from mesh geometry. It follows the innermost
eight sites at the exact softened circular frequency while holding their
initial masses fixed. Over 61 samples from `t=0` to `1.2`, plus `t=2.2214`,
`6.5`, `10`, `120`, and `600`, the sparse narrow/wide first-ring cells have
maximum individual density changes below `6.05e-4`; at `t=600` their
mass-weighted absolute change is about `8.4e-5`. Eight angularly analogous
full q=2 cells have maximum individual changes up to `5.29e-3` and
mass-weighted change `2.28e-3` at `t=600`. The local core polygons stay
within radius `0.157`, so sites beyond the `r<0.5` tessellation cannot
affect them. These are **geometry-only** results. They show that the sparse
core is geometrically quiet under exact orbits despite its large point-versus-
cell-mean gravity difference; only native evolution can establish which
effect controls heating or positivity.

The better guarded 55,931-site wide IC completed independent
CUDA+MPI+OpenMP H200 preflight **436967** to exact `t=.001`. Its
native audit checks the expected site count, fixed mass, positive state,
W=0, disk orbit, and both restored guard-ring mass sums. It uses the same
sparse-capable binary and physical parameters as the 54,312-site narrow
preflight; no wide long-run result exists yet.
Its originally staged 8-cell `t=.825` and `t=10` jobs **436969** and
**436970** were canceled before execution after the periodic-halo screen.
The replacements below preserve their physics, ID-wise mass, positive-state,
W=0, disk-accuracy, and radial audits. No wide native evolution result exists yet.

## Sparse q2 periodic-halo screen and replacement continuations (2026-10-27)

Both completed sparse preflights passed the then-current, incomplete native
`t=.001` checks; the global-area/mass failure is documented above. To
avoid investing H200 time in an insufficient periodic halo, the independent
[Voronoi auditor](audit_kepler_periodic_halo_geometry_20261027.py) was
extended to accept complete sparse native snapshots with exact ID coverage.
It wrapped each candidate's own `t=0` generators after exact softened
circular motion and reconstructed edge-cell polygons with a 2.5-unit image
margin at `t=.825` and `2.221`. Narrow (54,312 sites) and wide (55,931)
give the same extremal width demand and violation counts:

| Ballistic time | Minimum integer halo cells | Violations at 8 | At 12 | At 28 | At 32 |
|---:|---:|---:|---:|---:|---:|
| 0.825 | 12 | 8 | 0 | 0 | 0 |
| 2.221 | 29 | 32 | 24 | 4 | 0 |

For narrow at `t=6.5`, a wider 2.95-unit image margin gives minimum 26
cells; at `t=10`, minimum 16 cells. Both have zero 32-cell violations.
The margin exceeds twice the largest audited circumradius at each sampled
time. See the [narrow early](KEPLER_Q2_PROFILE_BUFFERED_BALLISTIC_HALO8_T0825_WIDTH_AUDIT_2026-10-27.json),
[narrow intermediate](KEPLER_Q2_PROFILE_BUFFERED_BALLISTIC_HALO8_T2221_WIDTH_AUDIT_2026-10-27.json),
[wide early](KEPLER_Q2_PROFILE_WIDE_BALLISTIC_HALO8_T0825_WIDTH_AUDIT_2026-10-27.json),
[wide intermediate](KEPLER_Q2_PROFILE_WIDE_BALLISTIC_HALO8_T2221_WIDTH_AUDIT_2026-10-27.json),
[narrow t=6.5](KEPLER_Q2_PROFILE_BUFFERED_BALLISTIC_HALO32_T65_WIDTH_AUDIT_2026-10-27.json),
and [narrow t=10](KEPLER_Q2_PROFILE_BUFFERED_BALLISTIC_HALO32_T10_WIDTH_AUDIT_2026-10-27.json)
receipts. These exact-orbit results are geometric risk screens; the outer
periodic continuation is not a physically valid gas solution, and sampled
times do not prove intermediate native safety.

Narrow **437099** and wide **437101** had been queued with 12-cell halos
to `t=.825`, followed by **437100** and **437102** with 32-cell halos to
`t=10`. Their scripts and dependency checks passed, but the original
preflight gate itself was incomplete: it accepted grossly wrong initial
cell areas and fixed masses. All four jobs were **canceled before execution**
on 2026-10-27. Their scripts
([narrow `.825`](kepler_q2_profile_buffered_to0825_h12_20261027.slurm),
[wide `.825`](kepler_q2_profile_wide_t0825_h12_20261027.slurm),
[narrow `10`](kepler_q2_profile_buffered_to10_h32_20261027.slurm),
[wide `10`](kepler_q2_profile_wide_t10_h32_20261027.slurm))
are historical records and must not be run from the invalid preflights.
The global halo widths remain geometry risk screens, not native successes.

The exact `t=0` site counts show the performance cost directly. For narrow
sparse q2 (54,312 real cells), 8/12/32-cell halos geometrically add
504/861/53,213 periodic images, giving 54,816/55,173/107,525 total sites.
For wide sparse q2 (55,931 real), the corresponding totals are
56,435/56,802/111,801. Thus 12-cell padding is nearly free on these sparse
outer layouts, but 32-cell padding approximately doubles the tessellation
site input. For comparison, full q2 has 147,456 real cells and 53,248
periodic images at 8 cells versus 221,330 at 28 cells, so its total input
grows from 200,704 to 368,786 sites. These are geometry counts from the
complete native initial snapshots, **not measured wall-speed ratios**.
Only a small set of outer cells violates the halo bound at the sampled
times; a future local fallback for those cells is preferable to a permanent
wide global halo if timing confirms a major cost.
The [local periodic-image fallback design](KEPLER_LOCAL_PERIODIC_FALLBACK_DESIGN_2026-10-27.md)
records an independently reconstructed reduced-image geometry at `t=2.221`,
the exact site counts, and the CUDA/host/MPI contracts needed before an
opt-in native implementation can be tested. It is not a measured speedup.

The existing GPU geometry-trial/RK4 transaction is not an equivalent
automatic halo repair for these accepted runs. The active q2 runners set
`GFS_VARPRESS_GEOM_CFL=.1` and leave `GFS_RK4_TRANSACTION` off, whereas
`gfs_rk4_transaction.inc` requires `C=.2` when enabled. On a geometry
rejection it halves the trial timestep; a missing periodic image is not
repaired by halving time. A suitable future implementation must distinguish
halo clipping from other invalid geometry, restore the accepted state, add
the needed periodic candidates, and retry at the same physical timestep
without changing the accepted CFL or EOS. This is a source-contract finding,
not a tested implementation.

## Current local-sync A/B queue state (2026-10-27)

Slurm still reports evolved-state probes **437180** (host periodic refresh)
and **437181** (device-local refresh) as `PENDING (QOSGrpGRES)`. They remain
paired H200:1 runs from the same accepted q1 step-8119 snapshot
(`t=0.6187644470884642`) to exact `t=.65`; the only requested physics-path
difference is `GFS_GPU_LOCAL_PERIODIC_SYNC=0/1`. The active q2 H200 baseline,
q1 H200 control, and refined A10 run are using the current GPU allocation.
The earlier 14-step `t=.001` pair's bitwise state match and 47.8% GPU-wrapper
reduction remain preliminary. Do not promote the local refresh to the long
q2 trajectory until the evolved checkpoint pair passes the exact endpoint,
ID-wise physical-field, and radial audits and supplies matched wall timing.

Before either pending job starts, an exact staging check found that the
radial auditor's imports `audit_kepler_seam_checkpoint_20261026.py` and
`audit_snapshot.py` were absent from both run directories. Both helpers have
now been copied into each directory and added to its `scripts.sha256` manifest.
The two runners pass `bash -n`, every staged checksum passes, and the exact
Slurm Python environment imports the radial auditor and both helpers. A fresh
Slurm check still showed 437180/437181 pending on QOS, so the repair happened
before execution.

At the latest queue/resource check, the H200 partition listed node `syn104`
with 8 H200 GPUs and `AllocTRES=gres/gpu=7`; 437076 and 437212 were using
H200s, while 437180/437181 still reported `QOSGrpGRES`. This is a point-in-time
resource snapshot, not a reservation guarantee. Keep the matched H200 pair
queued and avoid duplicating it on a different GPU model; neither probe has
started or produced an evolved-state result yet.

## Evolved q2 local-timestep budget at t=2.000 (2026-10-27)

The accepted native q2 checkpoint
`/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/halo12_after_failure_to0825/to2p221_h28/khout.026880.dat`
contains 147,456 cells at `t=2.0000159421`. A read-only analysis combines
the saved per-cell final-stage hydro CFL candidate with the existing softened
orbital limit, `min(dt_CFL,i, 0.06 (r_i^2 + 0.01175^2)^(3/4))`. Its receipt is
[`KEPLER_Q2_T200_LOCAL_DT_BUDGET_ALL_HORIZONS_2026-10-27.json`](KEPLER_Q2_T200_LOCAL_DT_BUDGET_ALL_HORIZONS_2026-10-27.json),
reproduced by [`kepler_local_dt_snapshot_budget_20261027.py`](kepler_local_dt_snapshot_budget_20261027.py).

At this evolved state the global orbital bound is `7.64202e-5`; the smallest
saved hydro CFL candidate is `1.22112e-4`. Only 5 cells have a combined bound
within twice the global step and 30 within four times. In the scientific disk
(`0.5 <= r < 2`), 49,038 cells have a combined local bound from `0.02175`
(minimum) to `0.10568` (median). This indicates substantial timestep
separation inside the disk, while the outer cells supply a second near-global
restriction through their hydro CFL values.

If these single-checkpoint values are frozen, power-of-two local scheduling
would give:

| Frozen interval | Global steps | All-cell updates at global dt | Block cell updates | Arithmetic ratio |
|---:|---:|---:|---:|---:|
| 0.221 (to t=2.221) | 2,892 | 426,442,752 | 8,469,338 | 50.35 |
| 4.5 (to t=6.5) | 58,885 | 8,682,946,560 | 171,228,127 | 50.71 |

The initial q2 `t=.001` frozen-state estimate was about 809-fold under the
same power-of-two cell-update metric. The drop to about 51-fold at `t=2`
shows that the outer low-CFL tail already erodes the ideal local-step
advantage; it is not evidence that the disk's timestep separation disappeared.
These counts assume each cell can independently use its saved bound. They do
not include neighbor synchronization, moving-face flux registers, remeshing,
GPU active-neighborhood costs, or changing CFL values, and are not predicted
wall speedups or evidence of a valid asynchronous hydro integrator. The
lineage audit below localizes the low-CFL tail to cross-shell pairs in the
outer background. The next useful checks are the exact limiting RK-stage face
for the `r=2.978` group and the mechanism that drives the outer tracers away
from circular motion; neither is established by the frozen-state arithmetic.

### Low-CFL cell lineage: outer-shell convergence, not yet a disk diagnosis

The low-CFL cells were traced by stable ID from q2 `t=.001` through the
existing `t=.825`, `1.000`, `1.500`, and `2.000` snapshots. The full time-series
receipt is
[`KEPLER_Q2_LOW_CFL_LINEAGE_TIMESERIES_2026-10-27.json`](KEPLER_Q2_LOW_CFL_LINEAGE_TIMESERIES_2026-10-27.json),
reproduced by [`kepler_q2_low_cfl_lineage_20261027.py`](kepler_q2_low_cfl_lineage_20261027.py).
At `t=.001`, all six tracked cells are on circular orbits (`v_r/v_c` near 0,
`v_phi/v_c` near 1), and each has the original `2.4414e-4` cell area. By
`t=.825`, IDs 143520 and 127494 already have large inward radial velocities
and strongly reduced or reversed azimuthal velocities. The values below are
`(r, v_r/v_c, v_phi/v_c)`:

| Cell ID | t=.001 | t=.825 | t=1.5 | t=2.0 |
|---:|---|---|---|---|
| 143520 | (3.2119, 0, 1) | (3.1689, -0.672, -0.750) | (2.9041, -0.694, -0.783) | (2.6947, -0.707, -0.813) |
| 127494 | (3.6360, 0, 1) | (3.3862, -0.964, -0.375) | (3.0146, -1.006, -0.397) | (2.7126, -1.038, -0.419) |
| 8512 | (3.3141, 0, 1) | (3.3141, 0, 1) | (3.2037, -0.782, -0.650) | (2.9780, -0.804, -0.674) |

By `t=2`, the first two tracers have formed close, converging pairs with
cells that remain close to circular:

| Pair at t=2 | Radius | Generator gap | Radial gap | Relative normal velocity | Saved local dt |
|---|---:|---:|---:|---:|---:|
| 10127 / 143520 | 2.6946 | 4.1497e-4 | 1.0417e-4 | -1.1771 | 2.74869e-4 |
| 76397 / 127494 | 2.7126 | 3.4502e-4 | 6.6408e-5 | -0.7242 | 3.65526e-4 |

For both pairs, `2*C*distance/v_sig` with the q2 `Courant=0.4` and
`v_sig=c_s,i+c_s,j-min(0,dv·e)` reproduces the saved local dt to the shown
precision. This is direct evidence that very close radial neighbors with
large converging velocities can create the outer CFL tail. The most restrictive
fourfold-symmetric group at `r=2.9780` also has a close pair, but its nearest
pair estimate is `3.963e-4`, larger than the saved minimum `1.221e-4`; the
exact limiting RK-stage face for that group remains unidentified.

For these six tracers the cell masses are unchanged, volumes fall to 0.50–0.60
of their `t=.001` values, and `P/rho^gamma` stays within about 0.3% of its
initial value. Their pressure and internal energy rise largely with the
compression, so these data do **not** show entropy-producing fake heating in
these particular cells. All cited pairs lie at `r>2.69`, outside the audited
scientific disk `0.5<=r<2`. They show that the outer background/guard region
develops near-coincident inter-shell pairs and loses circular motion; they do
not establish that the same mechanism causes inner-disk heating after `t=6`.

The evolved geometry confirms that close cross-shell pairs with strong
convergence can create a very small CFL, but it does **not** show that the IC
started with intrinsically too-close adjacent shells. At `t=.001`, IDs
10127/143520 began at radii 2.6947/3.2119 and IDs 76397/127494 at
2.7128/3.6360; both pairs had normal circular velocities and ordinary cell
areas. By `t=2`, they had converged to radial gaps of `1.04e-4` and `6.64e-5`
with large relative velocities. The initiating cause is therefore unresolved
outer-region radial drift/crossing, not a demonstrated initial spacing defect.
The next physics diagnostic is to repeat this lineage and pair audit at the
first `t>6` checkpoint and compare any gap collapse with the density/entropy
growth there. In parallel, the precise limiting stage-face for the `r=2.978`
group needs provenance before its saved dt is used as a block-step schedule
input.

### Same-shell versus adjacent-shell spacing in the scientific disk through t=2

The direct spacing comparison is now complete for q2 snapshots at
`t=.001`, `.825`, `1.000`, `1.500`, and `2.000`. The reproducible receipt is
[`KEPLER_Q2_RING_SPACING_AUDIT_2026-10-27.json`](KEPLER_Q2_RING_SPACING_AUDIT_2026-10-27.json),
from [`kepler_q2_ring_spacing_audit_20261027.py`](kepler_q2_ring_spacing_audit_20261027.py).
The IC has 68 initial-radius shells in `0.5<=r<2`, spaced uniformly by
`0.02209709`. In the disk, the median nearest distance within the same initial
shell is `0.01087755`; the median nearest distance to either adjacent initial
shell at `t=2` is `0.02218841`. Their ratio is **2.03985**, and **zero** of
the 49,038 current disk cells have an adjacent-shell candidate closer than
their same-shell neighbor. The ratio remains between `2.029` and `2.093` over
all cells and all five snapshots, with no adjacent-shell candidate closer.

So the proposed spacing imbalance is **not present in the q2 scientific disk
through t=2**: adjacent-shell neighbors are about twice as far as the nearest
same-shell neighbors, and this spacing pattern barely changes over the
measured interval. This is a generator-distance audit, not an exact face-flux
or t>6 causal test. It does not rule out later shell compression, but it gives
no basis for increasing inner resolution solely to widen the initial
inter-shell spacing. The next comparison should be made on the running q2
trajectory at and beyond `t=6`, before choosing that IC change.
