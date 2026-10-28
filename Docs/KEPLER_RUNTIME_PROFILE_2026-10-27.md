# Weighted Kepler runtime profile — 2026-10-27

## Measurement

Source: active Slurm job 437531, A10, one MPI rank, two OpenMP threads, v5
binary `feae0a48f737a7eab215f2e31be1d7938f45f894376dcdfda177f7e09f57d26b`.
The figures below come from the latest 100 complete `[CPU step]` records and
the corresponding 400 GPU tessellation records in the live log.

| Quantity | Recent measurement |
|---|---:|
| Printed `[CPU step] total` | median 1.00 s; mean 1.037 s |
| `update` phase | median 0.70 s |
| `force` phase | median 0.10 s |
| `final` phase | median 0.30 s |
| GPU tessellations per RK4 step | 4 |
| One GPU tessellation | median 169.8 ms; mean 170.9 ms |
| One `[GPU_VOLW2]` stage | median 3.38 ms |

The time marker reached `t=14.1415` at step 2136; the preceding energy ledger
at step 2135 reported relative total-energy drift `-3.339e-9`, zero
nonpositive internal-energy cells, and zero floor injection. These are
current-run diagnostics, not the final t=20 audit.

This measurement is from a different v5 executable than the archived source
bundle audited below. Its printed timer scope has not been verified against
that exact binary's source, so retain these as logged totals rather than
treating them as independently measured wall time.

## What the timing says

The four tessellation calls account for about 0.68 seconds per RK4 step, close
to the measured 0.70-second `update` phase. One weight stage takes about
fifty times less time than one tessellation. The next optimization should
target geometry preparation and tessellation data movement, rather than first
rewriting the weight solver.

The current source tree's RK4 path calls `updateDenW2Pressure2DBlend` at K1,
K2, K3, and K4. That function rebuilds particle cell lists and prepares GPU
inputs at each stage. The current tree also builds host cell-CSR data with
serial count/fill passes and packs the AoS particle records into host SoA
arrays. These paths are plausible targets for parallel packing/CSR generation
or persistent device-resident stage state.

The running v5 executable was built from an archived source bundle plus
several patches. The source tree is dirty and may have changed since that
bundle was built. Any driver code change must first compare the proposed edit
against the exact archived source used for the target binary. The runtime
measurements above remain valid regardless of that source comparison.

## Recommended driver investigation

1. Add separate timings for linked-list construction, CSR count/fill, AoS→SoA
   packing, host-device transfers, tessellation kernels, and downloads. The
   current aggregate timer cannot distinguish these costs.
2. Before adding parallel work, check the exact target source: this archived
   build already parallelizes SoA packing and per-cell CSR extraction. Profile
   any serial compaction or padding-refresh remainder first. Do not parallelize
   linked-list head insertion directly; that would race.
3. Consider keeping RK4 stage arrays resident on the GPU so each stage does not
   round-trip all particle fields through host AoS/SoA. Preserve the four
   geometry evaluations unless a mathematically equivalent integrator is
   separately implemented and validated.
4. Use the pending t=20→30 H200 OMP8 run and the matched one-rank/two-rank H200
   probe to quantify hardware/thread effects before changing the physical
   operator.

No source changes were made as part of this profile. The pending MPI probe is
documented in `KEPLER_MPI_SCALING_PROBE_2026-10-27.md`.

## Full-q2 H200 profile

The separate full-domain 384² W0 job **437076** ran on one H200, one MPI
rank, and two OpenMP threads. A sample of its last 100 completed
`[CPU step]` records near `t=1.93245` gives median/mean printed phase sum
`2.90/2.895 s`, with median update `1.6 s`, force `0.8 s`, and finalization
`0.5 s`. The corresponding 400 `[GPU tess]` records (four per RK4 step)
have median `302.5 ms` and mean `304.5 ms` each. The log separately reports
`[GPU step]` sub-timers; these counters have different scopes and should not
be summed as disjoint wall-time components without a finer profile.

At this sample, the step costs about 2.8 times the A10 55,351-site weighted
case above, while the full-q2 run has 2.66 times as many sites. These are
different GPU models, particle counts, and operators, so this ratio is not a
controlled hardware comparison. It does show that the four tessellation
calls and CPU update/force phases are both material costs. The next driver
profile should separate linked-list construction, force/moment work, host
packing and transfers, and CUDA synchronization in the exact archived source
used by job 437076. The present working source is dirty and was not changed.

## Live full-q2 timestep update

A later log read of job **437076** found it still running, with the latest
time marker at `t=1.95547`, global step 26,248. Over the last 100 time markers,
the accepted `dt` had median `7.642e-5`, mean `7.103e-5`, minimum `1.410e-5`,
and maximum `7.642e-5`. A few short low-dt excursions recovered to the prior
value within subsequent steps; this is not evidence of a sustained timestep
collapse. The job sets `GFS_VARPRESS_GEOM_CFL=0.1` but does not enable
`GFS_VARPRESS_GEOM_AUDIT`, so the log cannot identify which cell or limiter
caused the transient minima. Keep the accepted timestep and do not loosen the
geometric-CFL coefficient on this evidence. If low values become sustained,
use a matched checkpoint with the opt-in geometry audit to record the limiting
cell and swept-volume rate before changing the timestep rule.

## Exact timer scope and executable provenance correction

The archived source bundle used by the active polar runs is present at
`/tmp/Kepler_weighted_native_20261026`. Its executable SHA-256
`7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284`
matches the q1 437212, q2 437076, and refined A10 437771 run binaries.
Thus the timer interpretation below is tied to the actual executable source,
not to the unrelated project-root `src` tree.

In `Exam/exam_gpu_extract.c`, the printed `[GPU step]` columns map to these
intervals: `extract = t2_soa - t1_extract`; `SoA = t3_upload - t2_soa`;
`upload = t4_kernel - t3_upload`; `kernel = t5_download - t4_kernel`;
`download = t6_writeback - t5_download`; and `writeback = t7_done -
t6_writeback`. In the ordinary extraction path, `SoA` includes
`fillParticleSoA` and particle/face uploads. If GPU tessellation retained
faces on device, it can instead include host periodic-padding refresh,
packing, and particle upload. Thus `SoA` is a combined staging interval, not a
pure conversion timer. The `upload` label is misleading: that interval is
parameter preparation before the kernel call, not the host-to-device
transfer. `extract` is branch-dependent: it can cover CSR face extraction,
fused scatter, or only device-face bookkeeping. The `(full)` label reports
`fused == false`; it does not reveal `tess_faces_on_device`. The `writeback`
interval includes result copy-back, the MPI timestep reduction, and cleanup.

Accordingly, the live ~100 ms `SoA` value establishes a costly combined
packing/refresh/transfer interval, but it does not isolate which component
dominates or prove that the device-face reuse branch is inactive. A useful
next profile must split host padding refresh, `fillParticleSoA`, and each CUDA
upload, and print the branch flag. No active solver source or run binary was
changed for this audit.

The exact archived source already parallelizes `fillParticleSoA` with an
OpenMP static loop and `extractFaceCSR_2D` with an OpenMP dynamic face loop and
per-thread buffers. Therefore the earlier generic suggestion to first
parallelize packing/CSR does not describe missing parallel loops in this
binary. It remains possible that compaction, host padding refresh, or transfer
dominates, but the existing combined timer cannot tell. The OMP=2 production
setting and short A10 OMP probes also cannot establish the best H200 thread
count for this full q1/q2 branch.

## Full-q2 W=0 live profile update (2026-10-27 21:23 KST)

At `t=2.19632`, a fresh rolling sample from the same 147,456-cell H200 job
437076 gives 3.000 s median / 3.039 s mean for the printed phase sum across the last
100 records. Median phase times are `update=1.6 s`, `force=0.9 s`, and
`final=0.3 s`; 400 GPU tessellations have a 314.3 ms median and 315.4 ms
mean. The accepted `dt` is usually `7.64202e-5` but still has intermittent
short dips. The runner requests OMP=2, while the exact archived source contains
many OpenMP loops. This confirms an optimization target but does not prove
OMP=8 or multi-rank scaling; the pending 4-GPU `t=.003` startup test uses a
different 256² IC on a 24-by-24 domain (`Rin=2`, `Rout=8`, `epsilon=.047`)
and is not a representative scaling comparison. A
matched thread-count benchmark on this full-q2 state is still needed before
changing the long-run configuration.

## Matched full-q2 OMP pair staged

Jobs **437835** and **437836** are staged against the exact accepted q2
`t=2.221` restart, both to `t=2.231`. The first uses one MPI rank, one H200,
and OMP=2/2 CPUs; the second uses one MPI rank, one H200, and OMP=8/8 CPUs.
The second waits `afterok:437835`, and the first waits `afterok:437076`; both
are deliberately low priority. Their runner also requires the exact native
and radial endpoint audits. This is the first matched H200 thread-count
measurement for the full 384² W=0 branch; its result is not yet available.

## RK4 per-cell stage-loop audit in the exact archived source

For executable SHA-256
`7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284`, the K1–K4
per-particle RHS loops in `Exam/exam.c` have no OpenMP directives. Each stage
calls `kepler_accel` and `phase1_ie_stage` once per site. At 147,456 sites,
that is 589,824 calls to each helper per RK4 step. In the Kepler branch,
`kepler_accel` performs up to four `getenv` calls and three `atof` conversions
per call for run-invariant central-mass settings; `phase1_ie_stage` also
re-reads and parses `GFS_THERMAL64` per site and stage. These are source-level
operation counts, not a measured share of runtime.

The same source prints `[CPU step] total` as
`_t_update + _t_force + _t_post + _t_fin`. Those timers cover density/weight
updates, force calls, `postStage`, and finalization intervals, but exclude the
per-particle RK4 stage arithmetic between them. Therefore the q2 values above
are instrumented-phase sums, not full RK4 wall time. Full step wall time still
needs an external measurement over a known batch of steps.

Driver-side optimization candidates are to parse immutable environment
settings once before stepping and pass cached values to these helpers, then
profile OpenMP static scheduling for the independent per-site RK4 loops after
checking helper and diagnostic state for thread safety. No solver source was
changed, and no speedup is claimed from this audit.

## Full-domain Q1 H200 wall-clock sample (2026-10-27 22:36 KST)

The one-H200, OMP=2, 147,456-cell W=0 Q1 segment **437212** restarted at
global step 10,856 (`t=0.825`). At the 22:36 KST log read it had reached
global step 16,652 (`t=1.25808`), or 5,796 accepted steps. Slurm elapsed time
was 4:14:27, giving about **2.63 wall seconds per accepted step** averaged
over the whole segment. This includes startup, checkpointing and output; it
is not an isolated kernel benchmark. The latest printed `[CPU step] total`
was about 2.2 s, but the exact source shows that value is only the sum of
instrumented phases and excludes per-particle RK4 stage arithmetic.

At the observed average rate and usual `dt=7.642e-5`, the remaining
`t=1.258→6.5` interval would take roughly 51 wall hours if throughput and
timestep remain unchanged. This is only a capacity estimate. It reinforces
that profiling and speeding the native RK4 path matters for the long Kepler
benchmark; it is not evidence that the current full-domain branch has passed
the `t>6` physics gate.

## Fresh live samples for the t=2.221 checkpoints (2026-10-27 23:16 KST)

The H200 Q1 job **437212** is still RUNNING. Its latest log sample reached
`t=1.3168`, global step 17,452, from restart step 10,856 (`t=.825`). Slurm
elapsed time was 4:53:43, or about **2.67 wall seconds per accepted step**
averaged over this segment. The most recent RK4 ledger has relative total-
energy drift `-9.04e-11`, `dEhyd/Ehyd0=1.02e-4`, and zero nonpositive-
internal-energy and floor-injection counts. At the measured average rate,
the remaining `t=1.3168→2.221` segment is roughly 9 wall hours; this estimate
does not include the later continuation to `t>6`.

The refined W=0 A10 job **437771** is also RUNNING. Its latest log sample
reached `t=.470297`, global step 6,165, from restart step 14 (`t=.001`).
Elapsed time was 4:38:03, or about **2.71 wall seconds per accepted step**.
The latest ledger has relative total-energy drift `-1.37e-11`,
`dEhyd/Ehyd0=1.91e-5`, and zero nonpositive-internal-energy and floor-
injection counts. At this average rate, `t=.4703→2.221` is roughly 17.3 wall
hours. Its 24-hour wall limit leaves headroom if the rate remains similar.

The W=0 run's latest `[GPU step]` samples report 98–130 ms in the combined
SoA/staging interval and 1.7–7.8 ms in the force-kernel interval. The same
outer steps show four `[GPU tess]` records around 125–139 ms each. These
instrumentation scopes are not disjoint, so they must not be summed. They do
identify host staging and repeated tessellation as measured optimization
targets; the next useful profile must separate padding refresh, particle
packing, each transfer, tessellation, and force work for the exact run binary.
The job input manifests record binary hashes but not source hashes. The
separate archived bundle at `/tmp/Kepler_weighted_native_20261026` is mapped
to the live binary SHA-256 `7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284`,
and its `Exam/exam.c` hash matches the preserved overlay. The provenance comes
from that separate bundle mapping, not from the per-job manifests alone.

Neither live job has reached `t=2.221`, and neither is yet a new `t>6` result.

## Consequence for the GIZMO `t=600` endpoint

Both live segments advance about `0.10` code-time units per wall hour at
their measured one-GPU throughput. Naively extending that rate to the
Figure-7 physical endpoint `t=600` gives about **247–249 days on one GPU**,
before endpoint audits or any failures. This is a capacity estimate, not a
prediction of late-time timestep behavior. It shows that reducing the present
per-step overhead alone is unlikely to make the target practical: the plan
also needs a validated increase in stable timestep or local time integration,
with conservation and the historical `t>6` stability gates preserved. A
multi-GPU speedup must be measured, but cannot be presumed to close this gap.

## RK4 RHS loop safety review for the exact archived Q1 source

In the source overlay matched to the live Q1 executable, the four K1–K4
per-site RHS loops at `Exam/exam.c:10907, 10986, 11065, 11144` have no OpenMP
directive. Each iteration writes only the current site's RK4 stage fields and
uses local temporaries. The called `kepler_accel` helper only reads the
immutable environment and writes through that site's output pointers;
`phase1_ie_stage` reads that site and returns a scalar. The stage audit
accumulators run serially after each loop. This makes the four RHS loops
structurally parallel, provided particle storage and environment values remain
unchanged during the loop.

Two safeguards should precede an OpenMP implementation: parse and cache the
Kepler center/softening and thermal mode once per run instead of calling
`getenv`/`atof`/`atoi` per particle and stage; and call the lazy-cached
`sedov_phase1_on()` serially before entering a parallel region. Preserve the
implicit barrier between each stage loop and the following stage update.
This audit does not prove speedup or cover the other per-particle loops in the
RK4 integrator. The staged OMP=2/8 probe compares the unmodified baseline; it
cannot measure the benefit of adding parallel directives to these loops.

The two live Kepler runners **437212** and **437771** both explicitly unset
`GFS_LAGUERRE_ROTATION`. That switch changes a face-normal-speed correction
on weighted faces; it does not add intrinsic cell spin, a moment of inertia,
or a spin/energy ledger. The native spin acceptance gates remain outstanding,
and these baseline trajectories cannot validate a cell-spin cure for shear.

## Rolling live-log audit (2026-10-27 23:22 KST)

A read-only parse of the last 5,000 log lines yielded 294 accepted-time
markers for each live branch. Q1/H200 job **437212**, from `t=1.30924` to
`1.33109`, has median `dt=7.64202e-5`, mean `7.45918e-5`, and minimum
`8.3454e-6`; 11 of 294 steps were below `5e-5`, and only one was below
`1e-5`. The low sample is not a sustained timestep collapse. Its latest
ledger reports relative total-energy drift `-9.78e-11`,
`dEhyd/Ehyd0=1.055e-4`, and zero nonpositive-energy cells or floor injection.

Refined W0/A10 job **437771**, from `t=.462426` to `.484721`, has median
`dt=7.64202e-5`, mean `7.60947e-5`, and minimum `1.28395e-5`; only two of
294 steps were below `5e-5`. Its latest ledger reports total-energy drift
`-1.09e-11`, `dEhyd/Ehyd0=2.032e-5`, and zero nonpositive-energy cells or
floor injection.

Across this rolling sample the printed CPU phase-sum medians are 2.2 s for
Q1/H200 and 2.0 s for refined W0/A10. Their corresponding `[GPU step]`
medians are 236.8 ms and 196.3 ms; the combined SoA/staging medians are
127.2 ms and 98.4 ms, while the force-kernel medians are 1.7 ms and 7.7 ms.
These are distinct instrumentation scopes and different cell counts, so they
are not direct GPU comparisons. The small timestep dips do not explain the
current roughly 10 wall-hours-per-code-time throughput; repeated stage work
and host staging remain larger measured opportunities.

## Live job update (2026-10-27 23:30 KST)

An exact Slurm query still shows Q1/H200 **437212** and refined W0/A10
**437771** running. Their dependent continuations **437213**, **437772**, and
**437779** remain pending. The OMP2 and MPI2 probes (**437835**, **437892**)
are still queued for `QOSGrpGRES`, with their later probes dependency-held.
Neither live branch has reached its `t=2.221` gate.

The latest accepted Q1 marker is `t=1.34382` (step 17821). Around steps
6955–6958 its timestep briefly fell as low as `8.37e-8`, then recovered to
`6.50e-5` at step 6959 and the usual `7.642e-5` by step 6960. At the latest
marker, `dEtot/|Etot0|=-9.93e-11`, `dEhyd/Ehyd0=1.09e-4`, and the log reports
zero nonpositive internal-energy cells and zero floor injection. This was a
short timestep excursion, not sustained collapse; its cause is undetermined.

The latest refined W0/A10 marker is `t=0.497908` (step 6532), with
`dt=7.642e-5`, `dEtot/|Etot0|=-7.15e-12`, `dEhyd/Ehyd0=2.147e-5`, zero
nonpositive internal-energy cells, and zero floor injection. Continue the
existing runs and inspect the per-step timestep if the Q1 excursion recurs;
do not interpret this transient alone as a change in the physical stability
result.

## Source audit of the short Q1 timestep excursion (2026-10-27 23:35 KST)

The preserved source overlay matched to the live binary confirms that the
orbital limiter returns `eta (r^2+eps^2)^(3/4)` and reduces it over all cells.
For this run's `eta=0.06`, `eps=0.01175`, and central site at `r=0`, that
bound is `7.642e-5`. The accepted `dt` values around `1.3433` reached
`8.37e-8`, far below this bound, so the orbital limiter itself did not cause
the dip; a tighter force-step or integrator-side limit did.

The active CUDA force kernel forms a per-cell minimum from face-signal CFL,
viscous CFL, acceleration CFL, and entropy-change CFL (plus the optional
variable-pressure geometric CFL when compiled). It reduces these candidates
to one global timestep without reporting which candidate or cell won. The
CPU implementation has detailed `[HYDRO_CFL]`/`[VISC_CFL]` diagnostics, but the
live Q1 log shows the CUDA path and contains none of those records. Current
evidence therefore narrows the source to a non-orbital timestep restriction
but cannot identify its exact limiter or cell. A future code-side diagnostic
should record the winning limiter and cell index/state; no rerun is justified
from this brief, recovered excursion alone.

## Live progress update (2026-10-27 23:36 KST)

The next scheduler poll confirms **437212** (H200, elapsed 5:15:50) and
**437771** (A10, elapsed 5:00:10) are still running. The latest accepted-time
markers are `t=1.35114` and `t=0.505626`, respectively, both at
`dt=7.64202e-5`; neither has reached `t=2.221`. Their latest RK4 ledgers
report `dEtot/|Etot0|=-1.022e-10` and `-1.245e-11`, respectively, with zero
nonpositive internal-energy cells and zero floor injection. The held Q1/W0
continuations remain dependency-pending; OMP2/MPI2 remain `QOSGrpGRES`-pending.

## CUDA staging profile and pending local-sync gate (2026-10-27 23:39 KST)

A fresh eight-call tail sample from the active Q1/H200 log gives median
`SoA=126.6 ms`, wrapper total `232.2 ms`, and force kernel `1.7 ms` for
147,456 cells. The refined W0/A10 sample gives `SoA=93.4 ms`, wrapper total
`183.9 ms`, and kernel `7.6 ms` for 176,400 cells. These are per-wrapper
samples at different cell counts and on different GPUs; do not treat them as
a direct hardware comparison. In both, host-side data preparation occupies
about half of the wrapper, while the force kernel is a small share. The
separate full RK4 CPU-step timings are about 2.2 s (Q1) and 1.9 s (W0), so
wrapper time is not the whole-step wall cost.

Both live runners explicitly set `GFS_GPU_LOCAL_PERIODIC_SYNC=0`. Their
parameters satisfy the existing single-rank opt-in guard: `Kappa=0`,
`entropy_mode=0`, `av_mode=5`, and `cd_amax=cd_amin=0`. The exact source path
can replace the repeated host padding refresh with device-local periodic
alias refresh. The short `t=.001` pair already matched every saved physical
field bitwise and cut median wrapper time by 47.8%, but that is not enough to
promote it to the long trajectory. Matched evolved-state jobs **437180**
(host) and **437181** (device-local), from the same q1 checkpoint to `t=.65`,
remain queued on `QOSGrpGRES` as of 23:38 KST. Keep the active baselines
unchanged until that pair passes state and radial audits.

## Live horizon and wall-budget check (2026-10-27 23:40 KST)

At the latest exact poll, Q1/H200 **437212** had elapsed 5:20:03 and was at
`t=1.35778`; its restart began at `t=.825`. Refined W0/A10 **437771** had
elapsed 5:04:23 and was at `t=.51281`; its restart began at `t=.001`. Both
latest accepted-time markers use `dt=7.64202e-5`. Their latest RK4 ledgers
show `dEtot/|Etot0|=-9.67e-11` and `-1.29e-11`, respectively, with no
nonpositive internal energy or floor injection.

The measured averages are about `0.100` code-time/hour for Q1 and `0.101`
for W0. If those rates hold, Q1 needs roughly 8.6 more wall-hours to `t=2.221`
(about 14 hours total against its 24-hour limit); W0 needs roughly 16.9 more
(about 22 hours total against its 24-hour limit). These are straight-line
estimates only. W0 has limited wall-time margin, and its existing guarded
resume job remains the recovery path if the 24-hour segment ends first.

## Unreported CUDA-wrapper front-end cost (source audit, 2026-10-27)

The live wrapper prints its first subphase only after `t1_extract`; its total
clock starts earlier at `t0_total`. In the exact preserved
`exam_gpu_extract.c`, that pre-instrumentation interval includes the
non-fused linked-cell/padding setup, the per-cell `dt=1e10` initialization,
and context/capacity checks (source lines 692–755). Subtracting the displayed
subphase medians from the same eight-call samples leaves approximately 96 ms
per Q1 wrapper and 77 ms per W0 wrapper in this unreported front end. These
are approximate residuals, not isolated costs; the current log cannot assign
them among those operations.

This makes the optimization order clearer: keep the already paired local
periodic-refresh test for the measured SoA/padding cost, and separately split
the `t0_total`→`t1_extract` interval around linked-cell construction, `dt`
initialization, and context checks before changing those algorithms. This
would distinguish repeated host setup from actual CUDA kernel time.

## Latest scheduler sample (2026-10-27 23:45 KST)

The 23:45 KST scheduler poll still found **437212** and **437771** running;
**437180/437181** local-sync A/B remain queued for `QOSGrpGRES`. The latest
log markers read immediately afterward were Q1 `t=1.36548` and W0 `t=.520904`,
both at `dt=7.64202e-5`. Their latest ledgers remain free of nonpositive
internal energy and floor injection; relative total-energy drifts are
`-1.008e-10` and `-1.145e-11`. These runs have not yet crossed `t=2.221`.

## Native density image and separate mass/weight audit (2026-10-27 23:52 KST)

The stored `khmap.013200.ppm` corresponds to the Q1 checkpoint at
`t=1.0000125425`, but its grayscale output occupies only values 0–4 and is
not visually useful. I regenerated the figure from the validated native
checkpoint fields with the existing log-scaled density renderer:
[Q1 density at t=1.00001](KEPLER_Q1_H200_T1_DENSITY_2026-10-27.png).
The raster uses the cell `den` field and shows initial density, evolved
surface density, and relative disk-density change.

The complete-ID parser confirms all `w2` values are exactly zero in Q1 at
t=0 and t=1.00001. Cell masses are bitwise identical for all 147,456 cells
from t=0 to this checkpoint; independently, the refined W0 preflight has
bitwise-identical masses for all 176,400 cells from t=0 to t=.001. These are
configuration-specific observations: they confirm fixed mass in these two
runs, while `W=0` itself only specifies the geometry and does not imply a
mass-evolution rule.

For Q1 cells initially in `0.5 <= r < 2`, the t=1.00001 disk contains
48,680 cells. The mass-weighted signed density change is `3.73e-7`, the
mass-weighted absolute change is `4.97085e-4`, and the 95th percentile,
99th percentile, and maximum absolute cell changes are `1.677e-3`,
`2.936e-3`, and `5.993e-3`. This is an early-time checkpoint, not a t>6 or
GIZMO-endpoint result.

At this update, Q1/H200 **437212** is at `t=1.37659` and W0/A10 **437771** at
`t=.533054`; both remain at `dt=7.64202e-5`, with zero nonpositive internal
energy and zero floor injection. The t=2.221 gates and the evolved-state
local-sync pair remain pending.

## q2 t=2.221 endpoint and refreshed live state (2026-10-27 23:55 KST)

The native CUDA+MPI+OpenMP q2 H200 checkpoint `khout.030070.dat` at exact
`t=2.221` now has both the strict 15-check gate and the radial audit. In the
initial disk (`0.5<=r<2`), mass-weighted absolute density change is
`1.4081e-4`, absolute entropy-log change `3.5605e-4`, angular-momentum drift
`-6.39e-13`, and total-energy drift `-3.7761e-10`. In the published
`0.8<R<1.2` annulus, the velocity L1 error is `7.899e-9`; plateau-corrected
density L1 is `1.6166e-4`. All 15 listed checks pass, including exact fixed
mass and W=0, positivity, area tiling, and zero floor/dual-energy injection.
This establishes a good inner-disk checkpoint through `t=2.221`, not long-time
stability. The outer disk is already substantially different: mass-weighted
density change is `0.175` for `2.35<=r<2.75` and `0.334` for `2.75<=r<5`;
the latter bin's radial-speed change is `0.183`. The outer-boundary region
must remain visible in subsequent audits and GIZMO comparisons.

The checkpoint density figure is
[Q2 H200 density at t=2.221](KEPLER_Q2_H200_T2.221_DENSITY_2026-10-27.png).
It is rendered from the validated native `den` field with the existing
complete-ID, fixed-mass, and W=0 checks. The image shows the initial and evolved
surface density plus the disk's cellwise density ratio.

At the 23:54 KST scheduler sample, q2 continuation **437250** (`2.221→4.5`)
was pending for `QOSGrpGRES`, and **437251** remained dependent on it. Thus no
`t>2.221` continuation has started and no `t>6` result exists. Separate
controls are still running: q1/H200 **437212** had reached `t=1.38112`, step
18,330, after 5:35 elapsed; refined W0/A10 **437771** had reached `t=.537904`,
step 7,056, after 5:19 elapsed. Both retain `dt=7.64202e-5` and their latest
RK4 ledgers show zero nonpositive internal-energy cells and zero floor
injection. The local-sync evolved-state pair and OMP/MPI scaling probes remain
queued, so there is not yet a long-trajectory performance result for those
CUDA changes.

## Midnight live check and long-horizon throughput requirement (2026-10-28 00:00 KST)

The next Slurm sample shows q2 segment **437250** (`2.221→4.5`) still pending
with `QOSGrpGRES`; **437251** remains dependency-gated. The q1/H200 control
**437212** is RUNNING at `t=1.3879`, global step 18,422, 5:40:17 into its
restart from `t=.825`. Its latest ledger has `dEtot/|Etot0|=-1.084e-10`,
`n_ie_le0=0`, and zero floor injection. Refined W0/A10 **437771** is RUNNING
at `t=.546234`, step 7,165, with the same `dt=7.64202e-5`, no nonpositive
internal energy, and zero floor injection. The matched local-sync jobs
**437180/437181**, OMP probes **437835/437836**, and MPI probes **437892/437893**
are still queued behind H200 group GRES or their declared dependencies.

The measured q1 restart rate is `0.09925` code-time/hour, based on advancing
from `t=.825` to `1.3879` in 5:40:17 (7,566 accepted steps, 2.70 wall seconds
per accepted step). If that rate held, reaching `t=6.5` from this marker would
take about 51.5 hours, while reaching the paper's `t=600` would take about
6,031 hours (251 days). These are constant-rate projections, not runtime
promises; timestep changes can make them worse. They show that the present
baseline is not yet a practical path to the full paper endpoint. The queued
OpenMP/MPI and host/device-local-sync comparisons are therefore performance
gates to prioritize before expanding the long chain beyond its existing
staging. No throughput improvement is claimed yet.

## H200 node and group-limit audit (2026-10-28 00:03 KST)

Slurm reports the only node in the `h200` partition, `syn104`, configured with
8 H200 GPUs; 7 GPUs, 54/64 CPUs, and 784 GiB are allocated. The q2 continuation
**437250** still requests one H200 and remains pending with reason
`QOSGrpGRES`, even though the node report shows one GPU unallocated. This is a
Slurm group-resource-limit constraint rather than evidence that all physical
H200 cards are occupied. The exact association/QOS cap is not exposed by this
job-state output. Preserve the existing dependency chain and avoid duplicate
GPU submissions while that limit blocks the benchmark and its scaling probes.

## Source-level hot-path candidate (2026-10-28)

Read-only inspection on `syntax` reconfirmed that archived executable
`eunha2.exe` has SHA-256
`7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284`.
In its matching `Exam/exam.c`, `kepler_accel` at line 4528 reads
`EUNHA_IC`, `EUNHA_KEPLER_CX`, `EUNHA_KEPLER_CY`, and `EUNHA_KEPLER_EPS`
and converts the three numeric strings on every call. The four serial K1–K4
RHS loops call it once per real cell per stage. At 147,456 cells, that is
589,824 calls, 2,359,296 `getenv` calls, and 1,769,472 `atof` calls per RK4
step in the ordinary Kepler path. These are exact source-level operation
counts, not a measured fraction of wall time.

The same four loops call `phase1_ie_stage`, which calls
`gfs_thermal64_requested()` for every cell. In this binary's
`Exam/gfs_thermal64.h`, that inline helper performs an uncached `getenv` and
`atoi`. This adds 589,824 environment lookups and 589,824 integer parses per
RK4 step. The combined source count is therefore about 2.95 million
`getenv` calls, 1.77 million `atof` calls, and 0.59 million `atoi` calls per
step, before other helpers. Again, this is an operation count, not a timing
measurement.

The first driver-side optimization candidate is to parse the immutable Kepler
settings and thermal-mode flag once during run initialization, then pass
cached parameters/flags into the gravity and internal-energy helpers.
Preserve the current defaults and `EUNHA_IC` gating exactly. Benchmark the
same accepted checkpoint and compiler settings before and after; do not infer
a speedup from operation counts alone.

The K1–K4 loops are also serial and write per-cell, disjoint RK4 fields, so
OpenMP static scheduling is a plausible second candidate. Before adding it,
ensure optional `phase1_ie_stage` mode selection is initialized before the
parallel region; in particular `sedov_phase1_on()` uses a lazy static cache
that must be initialized serially first. Keep diagnostics/counters outside or
behind safe reductions. `fillParticleSoA` already has an OpenMP static loop and
`extractFaceCSR_2D` already uses OpenMP; those should not be treated as missing
parallelization. No solver source or active binary was changed in this review.

At 00:21 KST, SHA-256 of the executable in the live q1 job directory
`/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/halo28_from0825_to2221`
was checked directly and matched the hash above. Thus the source-count analysis
applies to running job **437212**, rather than only to an archived build. Its
captured CPU step reports 2.7 s total (`update=1.3 s`, `force=0.9 s`,
`final=0.4 s`); this identifies the host update/force path as a concrete
optimization target, while the source operation counts still do not quantify
how much of that measured time the repeated environment parsing consumes.

## H200 QoS cap confirmation (2026-10-28 00:05 KST)

`syntax` accounting query `sacctmgr show qos part_h200` returns
`GrpTRES=gres/gpu=6`. At the same check the queue had six RUNNING H200 jobs,
while `scontrol show node syn104` reported 8 installed GPUs and
`AllocTRES=gres/gpu=7`. The q2 root and scaling probes still report
`QOSGrpGRES`. The spare physical card therefore does not bypass the
partition-level six-GPU group cap. This replaces the earlier uncertainty
about the precise QoS limit; no unrelated allocation was changed.

## Live checkpoint (2026-10-28 00:07 KST)

q1/H200 **437212** remains RUNNING at `t=1.39653`, global step 18,535. Its
latest RK4 ledger at step 7,679 has `dEtot/|Etot0|=-1.110e-10`,
`dEhyd/Ehyd0=1.236e-4`, no nonpositive internal energy, and zero floor,
dual-energy, or source-floor injection. The refined 420² A10 job **437771**
is RUNNING at `t=0.55678`, step 7,303. q2 root **437250**, OMP2 probe
**437835**, and MPI2 probe **437892** remain pending with `QOSGrpGRES`; all
later intervals are still dependency-gated. The six-GPU H200 QoS cap remains
fully occupied, and no full-domain periodic result beyond `t=2.221` exists.

## Live checkpoint (2026-10-28 00:08 KST)

q1/H200 **437212** reached `t=1.39867`, global step 18,563; refined A10
**437771** reached `t=0.559295`, step 7,338. Both remain RUNNING and below
their exact `t=2.221` audit. The H200 queue still has six jobs RUNNING; q2
root **437250**, matched OMP2 **437835**, and MPI2 **437892** remain pending
with `QOSGrpGRES`. No interval beyond the already audited periodic q2 parent
has started.

## Live checkpoint (2026-10-28 00:10 KST)

q1/H200 **437212** reached `t=1.40050`, step 18,587; refined A10 **437771**
reached `t=0.561588`, step 7,368. Both remain RUNNING. The q2 root, OMP2,
and MPI2 jobs remain `QOSGrpGRES`; all later segments remain gated. At the
q1 restart's observed average rate (about `0.099` code-time/hour), its
remaining `t=1.4005→2.221` interval projects to about 8.3 wall hours. This is
only a rate-based estimate; the next H200 opening is determined by Slurm and
the six-GPU QoS cap, not by this projection.

## Live checkpoint (2026-10-28 00:11 KST)

q1/H200 **437212** reached `t=1.40211`, step 18,609; refined A10 **437771**
reached `t=0.563699`, step 7,396. Both remain RUNNING. The H200 partition
still has six active jobs; q2 root **437250**, OMP2 **437835**, and MPI2
**437892** remain `QOSGrpGRES`. The periodic q2 branch is still at its
audited `t=2.221` parent, with no `t>6` result yet.

## Live checkpoint (2026-10-28 00:12 KST)

q1/H200 **437212** reached `t=1.40402`, step 18,634; refined A10 **437771**
reached `t=0.565915`, step 7,425. Both remain RUNNING and below their
`t=2.221` gates. Six H200 jobs still occupy the partition QoS allowance;
q2 root **437250**, OMP2 **437835**, and MPI2 **437892** remain
`QOSGrpGRES`, with downstream segments dependency-held. The periodic q2
branch has not advanced beyond its audited `t=2.221` parent.

## Live checkpoint (2026-10-28 00:13 KST)

q1/H200 **437212** reached `t=1.40554`, step 18,654; refined A10 **437771**
reached `t=0.567749`, step 7,449. Both remain RUNNING. q2 root **437250**,
OMP2 **437835**, and MPI2 **437892** remain `QOSGrpGRES`, with downstream
segments dependency-held. The accepted periodic q2 endpoint remains `t=2.221`;
no later periodic checkpoint has been produced.

## Live checkpoint (2026-10-28 00:15 KST)

q1/H200 **437212** reached `t=1.40709`, step 18,675; refined A10 **437771**
reached `t=0.569736`, step 7,475. Both remain RUNNING. q2 root **437250**
and the exact-q2 OMP2/MPI2 probes **437835/437892** remain `QOSGrpGRES`;
the H200 partition still shows six running jobs. All post-root intervals
remain dependency-gated, and the periodic q2 branch has no checkpoint beyond
the audited `t=2.221` parent.

## Live checkpoint (2026-10-28 00:17 KST)

q1/H200 **437212** reached `t=1.41089`, step 18,726; refined W=0/A10
**437771** reached `t=0.574321`, step 7,535. Both remain RUNNING and have not
reached their `t=2.221` audit gates. The periodic q2 root **437250** and the
OMP2/MPI2 probes **437835/437892** remain pending with `QOSGrpGRES`; their
dependent segments are still held. The periodic branch therefore still has
no accepted checkpoint after `t=2.221`.

## H200 QoS turnover estimate (2026-10-28 00:18 KST)

The partition still has six of six QoS GPU slots occupied. In the live H200
queue sample, the two oldest allocations had each used about 1 day 15 hours
of a 2-day limit, leaving roughly 9 hours to their nominal wall limits. The
q1/H200 job **437212** is projected to reach its `t=2.221` stop in about 8.3
hours; if it completes normally, it could free a QoS slot around 08:40 KST.
Otherwise, either older allocation may free a slot near 09:20 KST if it runs
to its limit. These are turnover estimates only; early completion and Slurm's
next start are not guaranteed. The q2 root **437250** remains queued and has
not been duplicated or changed.

## Live checkpoint (2026-10-28 00:19 KST)

q1/H200 **437212** reached `t=1.41324`, step 18,758; refined W=0/A10
**437771** reached `t=0.577352`, step 7,575. Both remain RUNNING. q2 root
**437250** and OMP2/MPI2 probes **437835/437892** remain `QOSGrpGRES`; their
dependent segments remain held. No accepted periodic checkpoint beyond
`t=2.221` or `t>6` result is available yet.

The same log sample records a q1 CPU step at 2.7 s total (`update=1.3 s`,
`force=0.9 s`, `final=0.4 s`) and a refined A10 GPU step at 193.6 ms total
(`SoA=95.6 ms`, `kernel=7.6 ms`, `tessellation=115.2 ms`). These are
per-run observations from different resolutions/configurations, so they are
not a controlled device comparison; they do show that the q1 CPU update/force
path is a substantial optimization target for the requested end-to-end speed.

## Live checkpoint (2026-10-28 00:21 KST)

q1/H200 **437212** reached `t=1.41485`, step 18,779; refined W=0/A10
**437771** reached `t=0.579339`, step 7,601. Both remain RUNNING. The q2 root
**437250** and OMP2/MPI2 probes **437835/437892** remain `QOSGrpGRES`; the
dependent periodic segments remain held. The accepted periodic q2 endpoint is
still `t=2.221`, with no `t>6` result yet.

## Follow-up live log read (2026-10-28 00:21 KST)

The next accepted markers are q1/H200 **437212** at `t=1.41562`, step 18,790,
and refined W=0/A10 **437771** at `t=.580333`, step 7,614. Their latest RK4
ledgers report relative total-energy changes `-1.056e-10` and `-1.787e-11`,
respectively, with no nonpositive internal-energy cells or floor injection.
q2 **437250** and the OMP2/MPI2 probes **437835/437892** remain QoS-pending;
the main periodic branch has not advanced past its `t=2.221` checkpoint.

## Slurm start estimate and live checkpoint (2026-10-28 00:22 KST)

`squeue --start` estimated MPI2 probe **437892** could start at 09:20 KST,
still with reason `QOSGrpGRES`; it gave no start estimate for q2 root
**437250** or OMP2 probe **437835**. Treat 09:20 as a scheduler estimate,
not a reservation. In the same check, q1/H200 **437212** reached `t=1.41662`,
step 18,804, and refined W=0/A10 **437771** reached `t=0.581785`, step
7,633. Both remain RUNNING; no periodic checkpoint beyond `t=2.221` exists.
The latest RK4 ledgers in those logs report `dEtot/|Etot0|=-1.035e-10` and
`-1.988e-11`; both report zero nonpositive internal-energy cells and zero
floor injection. This confirms short-interval conservation only; it does not
establish orbital/density stability at `t>6`.

## Live checkpoint (2026-10-28 00:24 KST)

q1/H200 **437212** reached `t=1.41938`, step 18,841; refined W=0/A10
**437771** reached `t=0.585224`, step 7,678. Both remain RUNNING. Their latest
RK4 records show `dEtot/|Etot0|=-1.042e-10` and `-2.225e-11`, respectively;
both still report `n_ie_le0=0` and `floor_inj=0`. The q2 root and OMP2/MPI2
probes remain `QOSGrpGRES`, with no accepted periodic result past `t=2.221`.

## Scheduler start estimate correction (2026-10-28 00:24 KST)

A fresh `squeue --start` read at 00:24 returned `START_TIME=N/A` for q2 root
**437250** and the OMP2/MPI2 probes **437835/437892**; all three still show
`QOSGrpGRES`. The earlier 09:20 MPI2 estimate is no longer present and must
not be used as a current forecast. Slurm presently provides no reliable start
time for the H200 work; the q2 chain remains held at its audited `t=2.221`
parent.

## Live checkpoint (2026-10-28 00:25 KST)

The direct log read advances q1/H200 **437212** to `t=1.42058`, step 18,858;
its latest RK4 record has `dEtot/|Etot0|=-1.043e-10`, `n_ie_le0=0`, and
`floor_inj=0`. Refined W=0/A10 **437771** reached `t=.586905`, step 7,700.
Both remain RUNNING. Slurm still shows q2 **437250** and OMP2/MPI2
**437835/437892** pending with `QOSGrpGRES`; no periodic checkpoint after
`t=2.221` has been produced.

## Live checkpoint (2026-10-28 00:26 KST)

q1/H200 **437212** reached `t=1.42154`, step 18,872; refined W=0/A10
**437771** reached `t=0.588197`, step 7,717. Both remain RUNNING. A current
`squeue --start` query still returns no start time for q2 root **437250** or
the OMP2/MPI2 probes **437835/437892**, all pending with `QOSGrpGRES`.

## Live checkpoint (2026-10-28 00:28 KST)

q1/H200 **437212** reached `t=1.42344`, step 18,897. Its latest RK4 ledger
reports `dEtot/|Etot0|=-1.029e-10`, `n_ie_le0=0`, and `floor_inj=0`. Refined
W=0/A10 **437771** reached `t=.5907`, step 7,750. Both remain RUNNING. q2
root **437250** and OMP2/MPI2 **437835/437892** remain `QOSGrpGRES`; no start
time or periodic result beyond `t=2.221` is available.

## Live checkpoint (2026-10-28 00:27 KST)

q1/H200 **437212** reached `t=1.42351`, step 18,898; refined W=0/A10
**437771** reached `t=0.590776`, step 7,751. Both remain RUNNING. q2 root
**437250** and OMP2/MPI2 probes **437835/437892** remain `QOSGrpGRES`; the
dependent periodic segments remain held, with no result beyond the accepted
`t=2.221` endpoint.

## Live checkpoint (2026-10-28 00:29 KST)

Direct log-tail read at 00:28:50 KST shows q1/H200 **437212** at
`t=1.42498`, step 18,918, and refined W=0/A10 **437771** at `t=0.592763`,
step 7,777. Both remain RUNNING. q2 root **437250** and OMP2/MPI2 probes
**437835/437892** remain pending with `QOSGrpGRES`; no accepted result past
`t=2.221` is available.

## Live checkpoint (2026-10-28 00:30 KST)

q1/H200 **437212** reached `t=1.42666`, step 18,940; refined W=0/A10
**437771** reached `t=0.59498`, step 7,806. Both remain RUNNING. Their latest
RK4 records report total-energy relative changes `-1.045e-10` and
`-2.493e-11`, respectively, with zero nonpositive internal-energy cells and
zero floor injection. q2 root **437250** and OMP2/MPI2 probes **437835/437892**
remain QoS-pending; no result past `t=2.221` exists.

## Live checkpoint (2026-10-28 00:31 KST)

q1/H200 **437212** reached `t=1.42834`, step 18,962. Its latest RK4 ledger
reports `dEtot/|Etot0|=-1.047e-10`, `n_ie_le0=0`, and `floor_inj=0`.
Refined W=0/A10 **437771** reached `t=.597043`, step 7,833. Both remain
RUNNING. q2 root **437250** and OMP2/MPI2 probes **437835/437892** remain
`QOSGrpGRES`; no accepted periodic result after `t=2.221` exists.

## Live checkpoint (2026-10-28 00:32 KST)

q1/H200 **437212** reached `t=1.42865`, step 18,966; refined W=0/A10
**437771** reached `t=0.597501`, step 7,839. Both remain RUNNING. Latest
RK4 ledgers show total-energy relative changes `-1.052e-10` and `-2.824e-11`,
with zero nonpositive internal-energy cells and zero floor injection. q2 root
**437250** and OMP2/MPI2 probes **437835/437892** remain `QOSGrpGRES` pending.

## Live checkpoint and gate-time projection (2026-10-28 00:33 KST)

The direct queue/log check at 00:32:53 KST shows q1/H200 **437212** at
`t=1.43018`, step 18,988, and refined W=0/A10 **437771** at `t=0.599637`,
step 7,867; both remain RUNNING. Latest RK4 ledgers report
`dEtot/|Etot0|=-1.055e-10` and `-2.821e-11`, with `n_ie_le0=0` and
`floor_inj=0` in both. q2 root **437250** and OMP2/MPI2 probes
**437835/437892** are still `QOSGrpGRES` pending.

Using only the 00:24-to-00:32:53 progress as a linear estimate, q1 would need
about 10.8 more wall hours and A10 about 16.7 hours to reach `t=2.221`. Against
their 24-hour limits, A10's estimated margin is about 1.4 hours. This is a
rough rate projection, not a completion guarantee; the next checkpoints will
show whether its pace holds.

## Live checkpoint (2026-10-28 00:34 KST)

q1/H200 **437212** reached `t=1.43186`, step 19,010; refined W=0/A10
**437771** reached `t=0.601853`, step 7,896. Both remain RUNNING. Their
latest RK4 ledgers show `dEtot/|Etot0|=-1.046e-10` and `-2.862e-11`, with
zero nonpositive internal-energy cells and zero floor injection. `squeue
--start` still gives no estimate for q2 **437250** or OMP2/MPI2 probes
**437835/437892**, all pending with `QOSGrpGRES`.

The 00:24-to-00:34 linear pace still projects roughly 10.7 hours for q1 and
16.5 hours for A10 to reach `t=2.221`; A10 retains only about 1.6 hours of
estimated margin under its 24-hour limit. These estimates are approximate and
will be revised from later checkpoints.

## Live checkpoint (2026-10-28 00:35 KST)

At 00:35:18 KST, q1/H200 **437212** was at `t=1.43339`, step 19,030, and
refined W=0/A10 **437771** at `t=0.603840`, step 7,922; both were RUNNING.
Their latest RK4 records report `dEtot/|Etot0|=-1.017e-10` and
`-2.987e-11`; both still have `n_ie_le0=0` and `floor_inj=0`. q2 jobs remain
QoS-pending. Using the 00:24-to-00:35:18 average, estimated time to `t=2.221`
is about 10.6 h for q1 and 16.4 h for A10, leaving roughly 1.6 h before A10's
24-hour limit if that average continues.

## Live checkpoint (2026-10-28 00:37 KST)

At 00:36:57 KST, q1/H200 **437212** was at `t=1.43556`, step 19,059, and
refined W=0/A10 **437771** at `t=0.606591`, step 7,958; both remain RUNNING.
Latest RK4 records report `dEtot/|Etot0|=-1.026e-10` and `-2.968e-11`, with
zero nonpositive internal-energy cells and zero floor injection. The 00:24-to-
00:36:57 pace implies about 10.5 h for q1 and 16.3 h for A10 to reach
`t=2.221`; this leaves A10 roughly 1.7 h of its 24-hour limit, if the average
continues. q2 root and probes remain QoS-pending.

## Live checkpoint (2026-10-28 00:38 KST)

At 00:38:10 KST, q1/H200 **437212** reached `t=1.43716`, step 19,080, and
refined W=0/A10 **437771** reached `t=0.608655`, step 7,985. Both remain
RUNNING. Latest RK4 records show `dEtot/|Etot0|=-1.017e-10` and
`-2.918e-11`; both still report zero nonpositive internal-energy cells and
zero floor injection. q2 root **437250** and probes **437835/437892** remain
QoS-pending, and no result past `t=2.221` is available.

## Live checkpoint (2026-10-28 00:39 KST)

At 00:39:18 KST, q1/H200 **437212** reached `t=1.43869`, step 19,100, and
refined W=0/A10 **437771** reached `t=0.610642`, step 8,011. Both remain
RUNNING. Latest RK4 records show `dEtot/|Etot0|=-9.998e-11` and
`-2.918e-11`; both report `n_ie_le0=0` and `floor_inj=0`. The start-time
query returned `N/A` for q2 **437250** and probes **437835/437892**. At the
00:24-to-00:39 average pace, reaching `t=2.221` projects about 10.3 h for q1
and 16.2 h for A10, with about 1.8 h of A10's 24-hour limit left as margin.

## Live checkpoint (2026-10-28 00:40 KST)

At 00:40:27 KST, q1/H200 **437212** reached `t=1.44022`, step 19,120, and
refined W=0/A10 **437771** reached `t=0.612552`, step 8,036. Both remain
RUNNING. Latest RK4 ledgers show `dEtot/|Etot0|=-9.794e-11` and
`-2.857e-11`, with `n_ie_le0=0` and `floor_inj=0`. The recent overall pace
still suggests roughly 10.1 h for q1 and 15.9 h for A10 to reach `t=2.221`;
the q2 root and MPI/OpenMP probes remain queued at `QOSGrpGRES`.

## Live checkpoint (2026-10-28 00:42 KST)

At 00:42:08 KST, q1/H200 **437212** reached `t=1.44243`, step 19,150, and
refined W=0/A10 **437771** reached `t=0.615380`, step 8,073. Both remain
RUNNING. Latest RK4 ledgers show `dEtot/|Etot0|=-1.002e-10` and
`-3.070e-11`; each has `n_ie_le0=0` and `floor_inj=0`. `squeue --start`
continues to report no start estimate for q2 root **437250** or probes
**437835/437892**; all remain `QOSGrpGRES` pending. No result beyond
`t=2.221` is available.

## Matched OMP/MPI probe resource audit (2026-10-28 00:44 KST)

The staged short probes use the same accepted q2 `t=2.221` parent and stop at
`t=2.231`. The OMP submission script selects only `OMP_NUM_THREADS` (2 or 8);
although its `#SBATCH` default is two CPUs, the effective Slurm request for
job **437836** is eight CPUs. The `squeue` resource fields confirm **437835**
is 2 CPUs/1 H200 and **437836** is 8 CPUs/1 H200, so the OMP=8 case is not
oversubscribed. The MPI manifests confirm **437892** uses 2 ranks, 2 CPUs per
rank, and 2 H200 GPUs; **437893** uses 4 ranks, 2 CPUs per rank, and 4 H200
GPUs. Both point to the same parent receipt, initial condition, and `2.221 →
2.231` interval. Each pair has exact endpoint and radial-audit gates.

At the resource check, OMP2 (**437835**) and MPI2 (**437892**) were pending on
`QOSGrpGRES`; their OMP8/MPI4 successors (**437836/437893**) were dependency-
held. No scaling result has been measured yet. The configurations are
resource-matched as intended, and their low priority keeps them behind the
long q2 continuation.

## Live checkpoint (2026-10-28 00:43 KST)

At 00:43:16 KST, q1/H200 **437212** reached `t=1.44396`, step 19,170, and
refined W=0/A10 **437771** reached `t=0.617290`, step 8,098. Both remain
RUNNING. Latest RK4 records show `dEtot/|Etot0|=-1.020e-10` and
`-3.359e-11`; both have `n_ie_le0=0` and `floor_inj=0`. The long-window
average still projects about 10.2 h for q1 and 16.1 h for A10 to `t=2.221`;
q2 root and OMP/MPI probes remain `QOSGrpGRES` pending.

## Live checkpoint (2026-10-28 00:44 KST)

At 00:44:24 KST, q1/H200 **437212** reached `t=1.44549`, step 19,190, and
refined W=0/A10 **437771** reached `t=0.618994`, step 8,122. Both remain
RUNNING. Latest RK4 records show `dEtot/|Etot0|=-1.036e-10` and
`-3.362e-11`; both report `n_ie_le0=0` and `floor_inj=0`. q2 jobs remain
`QOSGrpGRES` pending, and no result beyond `t=2.221` is available.

## Live checkpoint (2026-10-28 00:46 KST)

At 00:45:52 KST, q1/H200 **437212** reached `t=1.44739`, step 19,215, and
refined W=0/A10 **437771** reached `t=0.621439`, step 8,154. Both remain
RUNNING. Their latest RK4 ledgers show `dEtot/|Etot0|=-1.022e-10` and
`-3.623e-11`, with zero nonpositive internal-energy cells and zero floor
injection. q2 root **437250** and OMP2/MPI2 probes **437835/437892** remain
`QOSGrpGRES` pending. The long-window average still projects roughly 9.9 h
for q1 and 15.8 h for A10 to `t=2.221`; these remain rough estimates.

## Live checkpoint (2026-10-28 00:47 KST)

At 00:46:58 KST, q1/H200 **437212** reached `t=1.44875`, step 19,234, and
refined W=0/A10 **437771** reached `t=0.623197`, step 8,177. Both remain
RUNNING. Latest RK4 ledgers show `dEtot/|Etot0|=-1.026e-10` and
`-3.619e-11`, with `n_ie_le0=0` and `floor_inj=0` for both. q2 root and
OMP2/MPI2 probes are still `QOSGrpGRES` pending; no result beyond `t=2.221`
is available.

## Live checkpoint (2026-10-28 00:48 KST)

At 00:48:06 KST, q1/H200 **437212** reached `t=1.45028`, step 19,254, and
refined W=0/A10 **437771** reached `t=0.625088`, step 8,202. Both remain
RUNNING. Latest RK4 ledgers show `dEtot/|Etot0|=-1.031e-10` and
`-3.529e-11`; both report `n_ie_le0=0` and `floor_inj=0`. q2 root and probes
remain QoS-pending. The 00:24-to-00:48 linear projection is about 9.9 h for
q1 and 15.8 h for A10 to `t=2.221`, leaving roughly 2 h of A10's 24-hour
limit if the rate stays similar.

## Live checkpoint (2026-10-28 00:49 KST)

At 00:49:14 KST, q1/H200 **437212** reached `t=1.45189`, step 19,275, and
refined W=0/A10 **437771** reached `t=0.626922`, step 8,226. Both remain
RUNNING. Latest RK4 ledgers show `dEtot/|Etot0|=-1.043e-10` and
`-3.670e-11`, with zero nonpositive internal-energy cells and zero floor
injection. q2 root and probes remain `QOSGrpGRES` pending. The 00:24-to-00:49
average projects about 9.8 h for q1 and 15.8 h for A10 to `t=2.221`; this is
a linear rate estimate only.

## Q1 short-timestep tail audit (2026-10-28 00:49 KST)

The last 100 `[RK4E]` records in q1 job **437212** have median `dt=7.642e-5`,
mean `7.54549e-5`, and minimum `1.474e-5` at step 8,377; only one of 100
records is below half the maximum. The local sequence is `7.642e-5` at step
8,375, `5.294e-5` at 8,376, `1.474e-5` at 8,377, then back to `7.642e-5`
at 8,378 and thereafter. This is a brief two-step reduction, not a sustained
timestep collapse. The runner sets `GFS_VARPRESS_GEOM_CFL=0.1` but does not
enable `GFS_VARPRESS_GEOM_AUDIT`; the log therefore does not identify the
limiting cell or cause. Keep the current CFL setting until a matched diagnostic
can identify a persistent limiter problem.

## Live checkpoint (2026-10-28 00:50 KST)

At 00:50:25 KST, q1/H200 **437212** reached `t=1.45340`, step 19,295, and
refined W=0/A10 **437771** reached `t=0.628831`, step 8,251. Both remain
RUNNING. Latest RK4 ledgers show `dEtot/|Etot0|=-1.037e-10` and
`-3.838e-11`; both report `n_ie_le0=0` and `floor_inj=0`. q2 root and probes
remain `QOSGrpGRES` pending. The 00:24-to-00:50 average still projects
roughly 9.8 h for q1 and 15.9 h for A10 to reach `t=2.221`.

## Live checkpoint (2026-10-28 00:52 KST)

At 00:51:33 KST, q1/H200 **437212** reached `t=1.45491`, step 19,315, and
refined W=0/A10 **437771** reached `t=0.630665`, step 8,275. Both remain
RUNNING. Latest RK4 ledgers show `dEtot/|Etot0|=-1.046e-10` and
`-3.848e-11`; both report `n_ie_le0=0` and `floor_inj=0`. q2 root and probes
remain `QOSGrpGRES` pending. Long-window estimates to `t=2.221` remain close
to 9.8 h for q1 and 15.9 h for A10, with A10 still within the 24-hour limit
if the observed average pace continues.

## Live checkpoint (2026-10-28 00:53 KST)

At 00:52:39 KST, q1/H200 **437212** reached `t=1.45625`, step 19,334, and
refined W=0/A10 **437771** reached `t=0.632346`, step 8,297. Both remain
RUNNING. Latest RK4 ledgers show `dEtot/|Etot0|=-1.048e-10` and
`-3.860e-11`; both report zero nonpositive internal-energy cells and zero
floor injection. q2 root and probes remain `QOSGrpGRES` pending. The long
window still projects roughly 9.8 h for q1 and 15.9 h for A10 to `t=2.221`.

## Q1 short-window throughput sample (2026-10-28 00:53 KST)

Between the 15:50:59 UTC and 15:52:53 UTC log reads, q1 **437212** advanced
from `t=1.45409`, global step 19,304, to `t=1.45647`, step 19,337: 33
accepted steps and `0.00238` code-time units in 114 wall seconds. This is
`0.0752` code-time units per wall hour, or about 3.45 wall seconds per
accepted step. Linear extrapolation of only this short window to `t=600` is
about 332 days. The restart-wide average, using `t=0.825` and the current
6:32:38 elapsed allocation, is faster at about `0.0965` code-time units per
hour, or about 258 days to `t=600`. These divergent estimates show material
rate variability; neither is a production ETA. Both confirm that the present
throughput is far below a practical `t=600` campaign without measured
optimization.

## Live checkpoint (2026-10-28 00:54 KST)

At 00:54:10 KST, q1/H200 **437212** reached `t=1.45812`, step 19,360, and
refined W=0/A10 **437771** reached `t=0.634571`, step 8,327. Both remain
RUNNING. Latest RK4 records show `dEtot/|Etot0|=-1.044e-10` and
`-3.905e-11`, with `n_ie_le0=0` and `floor_inj=0` for both. q2 root and
OMP2/MPI2 probes remain `QOSGrpGRES` pending, with dependent segments held;
no result beyond `t=2.221` exists.

## Q1 H200 recent phase profile (2026-10-28 00:54 KST)

From the most recent 100 `[CPU step]` records in job **437212**, the printed
phase medians are `update=1.4 s`, `force=0.9 s`, `final=0.5 s`, and
`total=2.8 s` (means `1.389`, `0.948`, `0.498`, and `2.838 s`). The
corresponding 400 GPU tessellation records have median `249.161 ms` and mean
`248.788 ms`; the 100 `[GPU step]` totals have median `237.2 ms` and mean
`234.756 ms`. The instrumentation scopes can be nested, so these values must
not be summed as independent costs. They show that the current H200 path still
spends substantial time in both update/geometry preparation and force work;
the matched OMP and MPI probes remain necessary before choosing a scaling
change.

## Live checkpoint (2026-10-28 00:55 KST)

At 00:55:19 KST, q1/H200 **437212** reached `t=1.45965`, step 19,380, and
refined W=0/A10 **437771** reached `t=0.636329`, step 8,350. Both remain
RUNNING. Latest RK4 ledgers show `dEtot/|Etot0|=-1.017e-10` and
`-4.056e-11`; each has `n_ie_le0=0` and `floor_inj=0`. q2 root and probes
remain `QOSGrpGRES` pending. The t=2.221 gate remains roughly 9.8 h away for
q1 and 16 h for A10 by the 00:24-to-00:55 average; still only a projection.

## Live checkpoint (2026-10-28 00:58 KST)

At 00:57:47 KST, q1/H200 **437212** reached `t=1.46293`, step 19,423, and
refined W=0/A10 **437771** reached `t=0.639920`, step 8,397. Both remain
RUNNING. Latest RK4 records show `dEtot/|Etot0|=-1.017e-10` and
`-4.368e-11`, with zero nonpositive internal-energy cells and zero floor
injection. q2 root and probes remain `QOSGrpGRES` pending. These early-time
conservation records do not establish long-term Kepler disk stability.

## Live checkpoint (2026-10-28 00:59 KST)

At 00:59:19 KST, q1/H200 **437212** reached `t=1.46492`, step 19,449, and
refined W=0/A10 **437771** reached `t=0.642169`, step 8,427. Both remain
RUNNING. Latest RK4 ledgers show `dEtot/|Etot0|=-1.018e-10` and
`-4.451e-11`; both report `n_ie_le0=0` and `floor_inj=0`. q2 jobs remain
`QOSGrpGRES` pending.

## Performance implication for the full Kepler target

The measured q1 throughput still puts a direct `t=600` campaign at roughly
260-330 wall days, depending on whether the restart-wide or short-window rate
is used. This is too slow for the requested campaign without code-level
acceleration. The measured q1 phase profile places substantial time in host
update and force work, while source inspection identifies repeated per-cell
environment parsing and serial K1-K4 RHS loops as concrete candidates. The
next implementation work should remain physics-preserving: cache immutable
run parameters first, then benchmark safe OpenMP static scheduling of the
independent RHS loops with initialized lazy state and thread-safe diagnostics.
The matched OMP/MPI probes are still queued, so scaling changes remain
unmeasured.

## Live checkpoint (2026-10-28 01:00 KST)

At 01:00:25 KST, q1/H200 **437212** reached `t=1.46629`, step 19,467, and
refined W=0/A10 **437771** reached `t=0.643851`, step 8,449. Both remain
RUNNING. Latest RK4 ledgers show `dEtot/|Etot0|=-1.017e-10` and
`-4.701e-11`, with `n_ie_le0=0` and `floor_inj=0`. q2 root and OMP2/MPI2
probes still have no Slurm start estimate and remain `QOSGrpGRES` pending.
The long-window projection is about 9.7 h for q1 and 16.2 h for A10 to
`t=2.221`; A10's 24-hour time-limit margin is roughly 1.4 h if the pace holds.

## Live checkpoint (2026-10-28 01:02 KST)

At 01:01:52 KST, q1/H200 **437212** reached `t=1.46799`, step 19,491, and
refined W=0/A10 **437771** reached `t=0.646143`, step 8,479. Both remain
RUNNING. Latest RK4 ledgers show `dEtot/|Etot0|=-1.038e-10` and
`-4.719e-11`; both report `n_ie_le0=0` and `floor_inj=0`. q2 root and
OMP2/MPI2 probes remain pending with `QOSGrpGRES`; dependent segments remain
held.

## Live checkpoint (2026-10-28 01:05 KST)

Slurm reported q1/H200 **437212** and refined W=0/A10 **437771** RUNNING at
01:04:41 KST. Their live logs read at 01:05:07 reached respectively
`t=1.47216` (step 19,548; RK4 ledger step 8,692) and `t=0.651398` (step
8,549; RK4 ledger step 8,535). The relative total-energy ledgers were
`-1.034e-10` and `-4.214e-11`; both reported `n_ie_le0=0` and
`floor_inj=0`. The A10 run is the refined zero-weight reference: setting its
Laguerre weights to zero recovers ordinary Voronoi geometry, while these
energy checks do not establish fixed cell masses. The `t=2.221` gate remains
incomplete; neither run has reached it.

## Live checkpoint (2026-10-28 01:06 KST)

At 01:06:23 KST, Slurm still showed q1/H200 **437212** and refined W=0/A10
**437771** RUNNING (elapsed 6:46:08 and 6:30:28). Their logs reached
`t=1.47354` (step 19,566; RK4 ledger step 8,710) and `t=0.653232` (step
8,573; RK4 ledger step 8,559). Relative total-energy drift remained
`-1.039e-10` and `-4.434e-11`; both still had `n_ie_le0=0` and
`floor_inj=0`. q2 root **437250** and OMP2/MPI2 probes **437835/437892** remain
pending for `QOSGrpGRES`; q2 continuation **437251** remains dependency-held.

## Live checkpoint and accepted q2 gate (2026-10-28 01:08 KST)

At 01:08:16 KST, Slurm reported q1/H200 **437212** and refined W=0/A10
**437771** RUNNING. Their latest logs reached `t=1.47587` (step 19,598;
RK4 ledger step 8,742) and `t=0.656213` (step 8,612; RK4 ledger step
8,598). Relative total-energy changes were `-1.051e-10` and `-4.753e-11`;
both had `n_ie_le0=0` and `floor_inj=0`. q2 continuation **437250** is
PENDING for `QOSGrpGRES`; **437251** remains dependency-held.

The full-domain periodic 384² W=0 q2 parent at exact `t=2.221` has both its
15-check global receipt and radial receipt. In the scientific disk
`0.5<=r0<2`, mass-weighted density change is `1.4081e-4`, absolute
entropy-log change `3.5605e-4`, orbit-position error `2.4628e-8`, angular
momentum change `-6.394e-13`, and total-energy change `-3.776e-10`; mass,
W=0, positivity, tiling, and no-injection checks pass. The same receipt shows
substantial change outside the disk: for `r0>=2.75`, density change is
`0.3343` and absolute radial speed is `0.1826`; the current outer periodic
edge has angular-momentum change `-0.603`. This is not yet proof that the
physical disk is contaminated, but it makes the matched radial bins and edge
metrics at `t=4.5` and `t=6.5` essential. The exact `t=600` GIZMO endpoint
remains untested.

## Live checkpoint (2026-10-28 01:10 KST)

At 01:09:53 KST, `sacct` confirmed full-domain q2 parent **437076**
COMPLETED after 16:30:41; its exact `t=2.221` global and radial receipts are
present and pass. The next segment **437250** is PENDING specifically for
`QOSGrpGRES`, while **437251** and later segments remain held on dependencies.
The independent OMP2/MPI2 probes **437835/437892** are also QOS-pending.

At 01:10:01 KST, q1/H200 **437212** and refined W=0/A10 **437771** logs
advanced to `t=1.47813` (step 19,628; RK4 ledger 8,772) and `t=0.65904`
(step 8,649; RK4 ledger 8,635). Their relative total-energy changes were
`-1.037e-10` and `-4.841e-11`, with `n_ie_le0=0` and `floor_inj=0` for both.
Neither has reached its own `t=2.221` endpoint. Thus the accepted q2 parent
supports disk stability through `t=2.221` only; the historical `t>6` gate and
the GIZMO `t=600` endpoint remain outstanding.

## Live checkpoint (2026-10-28 01:11 KST)

At 01:11:14 KST, q1/H200 **437212** and refined W=0/A10 **437771** remained
RUNNING at `t=1.47961` (step 19,648; RK4 ledger 8,792) and `t=0.660951`
(step 8,674; RK4 ledger 8,660). Relative total-energy changes were
`-1.038e-10` and `-4.876e-11`; both still report `n_ie_le0=0` and
`floor_inj=0`. The A10 recovery segment **437772** remains dependency-held
behind its parent, and q2's `t=4.5` segment **437250** remains pending on
`QOSGrpGRES`. No state beyond the accepted periodic q2 `t=2.221` endpoint has
been produced by this chain.

## Recovery guard and live checkpoint (2026-10-28 01:12 KST)

At 01:12:38 KST, q1/H200 **437212** and refined W=0/A10 **437771** remained
RUNNING at `t=1.48136` (step 19,671; RK4 ledger 8,815) and `t=0.663153`
(step 8,703; RK4 ledger 8,689). Their relative total-energy changes were
`-1.029e-10` and `-4.802e-11`; both had zero nonpositive internal-energy
cells and zero floor injection. q2's `t=4.5` job **437250** and matched
OMP2/MPI2 probes remain QOS-pending.

Inspected the actual A10 recovery runner for **437772**: it is `afterany` the
current **437771**, exits without a second evolution if the parent's exact
target audit passes, and resumes only after a wall-time stop from the newest
complete, hash-verified checkpoint. It refuses other parent failure states
and refuses to repeat from the startup checkpoint. The two-day recovery path
is therefore staged to preserve progress if the current 24-hour allocation
ends before `t=2.221`.

## Live checkpoint (2026-10-28 01:14 KST)

At 01:14:19 KST, q1/H200 **437212** and refined W=0/A10 **437771** were
still RUNNING at `t=1.48346` (step 19,699; RK4 ledger 8,843) and
`t=0.665828` (step 8,738; RK4 ledger 8,724). Their total-energy relative
changes were `-1.029e-10` and `-4.840e-11`, with zero nonpositive internal
energy and zero floor injection. q2 segment **437250** and OMP2/MPI2 probes
remain `QOSGrpGRES` pending; **437251** and the A10 recovery **437772** remain
dependency-held. Both active logs continue advancing, but neither has reached
their own `t=2.221` endpoint.

## Live checkpoint (2026-10-28 01:15 KST)

At 01:15:35 KST, q1/H200 **437212** and refined W=0/A10 **437771** were
RUNNING at `t=1.48514` (step 19,721; RK4 ledger 8,865) and `t=0.667891`
(step 8,765; RK4 ledger 8,751). Relative total-energy changes were
`-1.044e-10` and `-4.812e-11`; both still report positive internal energy
and zero floor injection. The q2 `t=4.5` root **437250** and OMP2/MPI2 probes
remain `QOSGrpGRES` pending; no downstream q2 segment has started.

## Exact binary/source provenance and performance candidate (2026-10-28)

The archived `syntax:/tmp/Kepler_weighted_native_20261026` bundle was found.
SHA-256 of its `eunha2.exe` and the executable copies in both active q1
(437212) and refined A10 (437771) directories is exactly
`7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284`.
Thus the co-located archived source snapshot is the correct executable-era
reference for these runs. Its Git HEAD is `472e7fcc3419b060d8fa0c74f0961f6c1f9f188b`,
but tracked files including `Exam/exam.c` and `Exam/exam_gpu.cu` have local
modifications; the executable is associated with that modified tree, not the
clean HEAD alone. The current workspace `src/Exam/exam.c` has a different
SHA-256 (`58629863fe38f2f6477863d9e3239c001c6b721f60b0b246976f05e1963cdc0f`)
from the archived exact-build-era `Exam/exam.c`
(`f38f780d7fbb1f6325417a2dec6fc4110685699192c6087e54cbfcb98a05134a`), so
source changes must start from the archived tree or first be reconciled
explicitly; no source or binary was changed here.

In the exact-build-era `Exam/exam.c`, `kepler_accel()` reads and parses the
Kepler environment variables on every call. It is called once per cell in
each of the four RK4 RHS stages; the four corresponding `for(i=0; i<np; i++)`
loops have no OpenMP pragma. At q1's 147,456 cells this is 589,824 helper
calls per step, at least 2.36 million `getenv` lookups and 1.77 million
`atof` conversions per step for the active Kepler profile. This is a
physics-preserving optimization candidate: parse immutable center, softening,
and mode once and pass cached values into the force helper. A second candidate
is OpenMP-static scheduling of independent RHS cell loops, after ensuring the
lazy `phase1_ie_stage` mode is initialized serially and diagnostics/counters
remain thread-safe. The exact source is now identified, but measured speedup
and unchanged physics still require a controlled rebuild and matched benchmark;
no such code change is made in this run.


## Controlled Kepler timestep-cap pair (2026-10-28 01:34 KST)

The short `dteta` comparison is now staged and submitted from the verified
384² q1 checkpoint at step 13,200 (`t=1.0000125425182276`) to exact `t=1.1`.
Both runs use the same archived executable (`7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284`),
parameters (`2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43`),
checkpoint (`24ed4e62d0160f26659a71f490360c42e196e6592c1152db2ea06dec9f01bea8`),
initial reference, 2-CPU A10 request, and audit path. Only
`EUNHA_KEPLER_DTETA` differs: baseline `0.06` versus `0` (no Kepler-specific
cap). The exact build's source documents that zero leaves the hydro timestep
alone and guards the cap with `eta > 0`.

Before submission, both run manifests passed `sha256sum -c`; restart header,
step, time, target, initial-file digest, and per-run parameter digest were
checked. The stale eta=0 parameter digest was corrected to the actual parameter
file hash, and both staged runners were synchronized with the local runner's
hash preflight. No solver source or binary was modified.

Slurm jobs **438003** (baseline) and **438004** (eta=0) were submitted to A10,
excluding `syn02` where the existing W=0 run is active. At the 01:34:16 KST
poll, 438003 was RUNNING on `syn03` (elapsed 17 s); 438004 was dependency-held
with `afterok:438003`. This is only a short timestep/performance diagnostic,
not a long-time physics result. The baseline landed on a mixed-use node, so
absolute throughput may be noisy; the sequential pair avoids overlap between
the two probe runs. Do not infer a `t>6` or GIZMO `t=600` result from this gate.


### Baseline probe checkpoint (2026-10-28 01:35 KST)

At 01:35:27 KST, job **438003** remained RUNNING on `syn03` at elapsed
1:28. Its restarted trajectory had reached step 13,227 / `t=1.00202` (27
steps beyond the checkpoint). Each accepted step was still `7.64202e-5`, the
central-cell Kepler cap; the latest RK4 ledger had
`dEtot/|Etot0|=1.619e-12`, `n_ie_le0=0`, and `floor_inj=0`. This confirms the
cap is active in the baseline, but says nothing yet about the eta=0 run or
long-time stability.


### Central-cell limiter audit (2026-10-28 01:39 KST)

Read the exact t=0 and t=1.0000125425 snapshots with the existing native
snapshot reader. Cell ID 73920 is at `r=0` initially and `r=8.26e-11` at the
restart, with density at the configured `0.01` floor and mass
`2.11432e-6`, only `1.61e-7` of total disk mass. Because the Kepler limiter
includes all cells in its global minimum, this floor cell sets
`dt=7.64202e-5`. The nearest cell with `rho>=0.02` is at `r=0.10937` and its
Kepler cap would be `2.189e-3` (28.6 times larger). However, there is real
mass inside `r=0.5` (2.47% at the restart), so dropping the entire inner region
from the limiter is not justified. The submitted `eta=0` branch will show
whether the native hydro CFL alone can advance the state while retaining the
physical audit gates; no density mask is introduced in this probe.

At 01:38:46 KST, baseline job **438003** was RUNNING at elapsed 4:47 and had
reached `t=1.00684`, step 13,292. Its accepted `dt` remained `7.64202e-5`.
The eta=0 comparison **438004** remains held behind `afterok:438003`.


## Isolated RK4/OpenMP performance candidate (2026-10-28)

To address the measured CPU-side cost without changing any active executable,
created `/gpfs/kjhan/LagForce/Kepler_weighted_native_20261028_perf_candidate`
from the verified build-era snapshot. The saved prebuild binary hash is the
same as the q1 baseline (`7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284`);
the candidate source began from exact `Exam/exam.c` hash
`f38f780d7fbb1f6325417a2dec6fc4110685699192c6087e54cbfcb98a05134a`.

The isolated source change caches immutable Kepler center/softening/mode once
per RK4 step, passes the existing thermal64 and Sedov-phase flags into the
cell helper and energy audit, and adds static OpenMP scheduling to the four
independent K1–K4 per-cell RHS loops. Per-cell equations and reductions are
unchanged. The candidate compiled and linked successfully with the active
Intel/CUDA toolchain; new binary SHA-256 is
`2fc4e12c98d89f498bb038682a10b60ce759942a894409448bd10ac6ce94d9f5`. `ldd`
resolves CUDA 13, Intel MPI, and Intel OpenMP runtime libraries. The compiler
reported one existing format warning at the unchanged `VORO_NPAD` diagnostic;
no candidate-specific compile warning appeared. Runtime equivalence and speed
are not yet established.

A matched run directory using this binary, the same `t=1.0000125425` restart,
`eta=0.06`, and exact `t=1.1` target passed all input hash checks. Slurm job
**438007** was submitted with `afterany:438004` and `--exclude=syn02`, so it
runs only after the timestep pair has ended and avoids the existing W=0 A10
job. This short run is the performance/code-path gate; it must pass the same
global and radial physics audits before using the candidate for a longer
continuation.


### Baseline transient-minimum checkpoint (2026-10-28 01:51 KST)

At 01:51:47 KST, job **438003** was RUNNING on `syn03` at elapsed 17:48,
`step=13550`, `t=1.02590`. The last four unusually small accepted steps were
`6.17e-5`, `1.39e-5`, `3.50e-6`, and `5.03e-6`; the next step recovered to
`7.642e-5`. Thus the low-dt episode is transient in this window, while most
steps still land exactly on the Kepler cap. The latest energy ledger remained
near-roundoff (`dEtot/|Etot0|=3.28e-13`, no nonpositive internal-energy cells,
no floor injection). Job **438004** remains afterok-held; the isolated RK4
candidate job **438007** remains afterany-held. The q2 t=4.5 segment **437250**
still waits on H200 `QOSGrpGRES`.


### Baseline timestep-cap occupancy

A streaming log analysis of the latest 200 accepted baseline steps (steps
13,394–13,593, `t=1.01456–1.02914`) found median `dt=7.64202e-5`, minimum
`3.02758e-6`, and 186/200 steps within 1 ppm of the `eta=0.06` cap
`7.64201789e-5`. The other 14 were lower transient CFL steps. This confirms
the Kepler cap limits the regular-step cadence in this window; the eta=0 pair
is needed to measure how much the native hydro candidate can safely take over.
It does not yet establish a speedup or accuracy at later times.


### Probe status (2026-10-28 01:55 KST)

At 01:55:33 KST, baseline **438003** remained RUNNING on `syn03` at elapsed
21:34, `t=1.03101`, step 13,620, with ordinary accepted `dt=7.64202e-5`.
The eta=0 run **438004** and the isolated OpenMP/cache candidate **438007**
remain dependency-held. Existing q1/H200 **437212** and refined W=0/A10
**437771** remain RUNNING. Full-periodic q2 continuation **437250** remains
pending on H200 `QOSGrpGRES`; **437251** is still dependency-held.


### Baseline probe and candidate review (2026-10-28 01:58 KST)

At 01:58:06 KST, baseline **438003** was still RUNNING on `syn03` at elapsed
24:18. Its log had reached `t=1.03498`, step 13,673, with the usual
`dt=7.64202e-5`. The latest RK4 energy ledger reported
`dEtot/|Etot0|=-4.58e-12`, `n_ie_le0=0`, and `floor_inj=0`. The four GPU
tessellations were each about 0.20–0.22 s; the printed CPU-step total was
2.4 s. These are local timing samples, not a sustained speedup estimate.
Jobs **438004** and **438007** remain dependency-held; the H200 q2 continuation
**437250** remains pending on `QOSGrpGRES`.

Reviewed the isolated candidate's exact diff against the active source on
`syntax`: its four new OpenMP loops fill only the corresponding per-cell RK4
stage fields, read that cell's state plus a read-only cached Kepler config,
and perform no shared reductions. The intended loop body therefore has no
obvious write race. Runtime equivalence, full-step timing, and any speedup
remain unverified until job **438007** runs and its audits are read.


### Pairwise audit staged; latest allocation check (2026-10-28 02:01 KST)

Added `audit_kepler_dt_eta_omp_candidate_pair_20261028.py` to compare the
matched baseline and OMP/cache candidate endpoints at exact `t=1.1`. It
requires both global and radial physics receipts to pass, equal step-sample
counts, bitwise equality in `x,y,vx,vy,mass,volume,ie,den,pressure,w2`, and a
lower candidate median CPU-step time before it reports the candidate gate as
passed. Python syntax compilation succeeded; the script has not run because
the candidate output does not exist yet.

At 02:01:37 KST, baseline **438003** remained RUNNING on `syn03` at elapsed
27:38, `t=1.03990`, step 13,738, with `dt=7.64202e-5`. Its latest energy
ledger reported `dEtot/|Etot0|=-3.56e-12`, `n_ie_le0=0`, and `floor_inj=0`.
The H200 partition had six running GPU allocations at this check (four
`gpu:1` and two `gpu:H200:1`), matching the six-GPU QoS cap; q2 continuation
**437250** therefore remains `QOSGrpGRES` despite spare physical cards on the
node. Jobs **438004** and **438007** remain dependency-held.


### Live run checkpoint (2026-10-28 02:03 KST)

At 02:03:50 KST, baseline **438003** remained RUNNING on `syn03` at elapsed
29:51, `t=1.04318`, step 13,781, with `dt=7.64202e-5`. Its latest RK4 ledger
reported `dEtot/|Etot0|=-4.67e-12`, `n_ie_le0=0`, and no floor injection; the
printed CPU-step total remained about 2.4 s. q1/H200 control **437212** was
RUNNING at `t=1.5436`, step 20,521. Refined W=0/A10 **437771** also remained
RUNNING. The eta=0 comparison **438004**, OpenMP/cache candidate **438007**,
and full-periodic q2 continuation **437250** had not started; 437250 remained
held by the H200 QoS cap. This checkpoint is not a long-time stability result.


### Geometry-moment reference gate (2026-10-28 02:07 KST)

Added a host `long double` reference for polygon area, centroid and central
second moments, with exact rectangle/asymmetric-triangle results, translation
invariance, and reversed/degenerate polygon rejection. It compiled warning-
free under C++17 and all controls passed. This is the independent oracle for a
later isolated CUDA-moment patch; the native solver still has no spin state or
GPU-returned Q tensor.

At 02:07:31 KST, baseline **438003** remained RUNNING on `syn03` at elapsed
33:32. Its latest log sample reached `t=1.04874`, step 13,854, with
`dt=7.64202e-5`, `dEtot/|Etot0|=-6.48e-12`, `n_ie_le0=0`, and no floor
injection. Jobs **438004/438007** remain dependency-held; q2 segment **437250**
remains QoS-pending on H200.


### Isolated GPU polygon-moment audit staged (2026-10-28 02:32 KST)

The isolated source worktree on `syntax`, based on commit
`472e7fcc3419b060d8fa0c74f0961f6c1f9f188b`, now contains an opt-in CUDA
polygon-moment diagnostic. It accumulates the signed area, centroid and
mass-central Q tensor in the tessellator's existing polygon-corner pass, then
copies the first GPU polygon array once and compares every local cell against
the independent host `long double` reference. It requires
`GFS_GPU_RAW_FACES=1`; it does not add spin state or alter the force update.

The CUDA object compiled, `Exam/libexam.a` rebuilt, and the full
`eunha2.exe` linked successfully with `make -j2 USE_CUDA=1
GFS_VARPRESS_W0=1`. The candidate executable SHA-256 is
`abd32739bc77ed0ede429ec1f14f795b5c00f50b0898d08c999503ddb1de3eb1`.
`ldd` resolved its runtime libraries. Build output contained existing format
warnings and the Intel Fortran warning for `-fno-strict-aliasing`; no new
moment-code warning appeared.

The one-step audit run **438014** was staged at
`/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/polygon_moment_audit_20261028_a10`.
It reuses the accepted 384² q1 checkpoint at step 13,200,
`t=1.0000125425182276` (checkpoint hash
`24ed4e62d0160f26659a71f490360c42e196e6592c1152db2ea06dec9f01bea8`) and
targets `t=1.0000625425182277`, so the test should execute only one short
step. Its inputs and binary passed the staged SHA-256 preflight, and all
dynamic libraries resolved. It requests one A10, one MPI rank and two OpenMP
threads, and is dependency-held with `afterany:438007`; `syn02` is excluded to
avoid the existing W=0 A10 run. The audit result is pending. This is a
geometry implementation gate only, not evidence of long-time Kepler stability.

The staging and batch files are
[`stage_kepler_polygon_moment_audit_20261028.py`](stage_kepler_polygon_moment_audit_20261028.py)
and [`kepler_polygon_moment_audit_20261028.slurm`](kepler_polygon_moment_audit_20261028.slurm).

### Dependency-chain checkpoint (2026-10-28 02:37 KST)

At 02:37:16 KST, baseline **438003** was still RUNNING on `syn03`, elapsed
1:03:17, at `t=1.09169`, step 14,423, with `dt=7.64202e-5`. Its latest
printed CPU phase sum was 2.4 s and GPU tessellations were about 216–230 ms.
The exact `t=1.1` endpoint remained ahead. **438004** and **438007** were still
dependency-held, so the new one-step audit **438014** had not started; it
remains afterany-held behind 438007. The existing q1 H200 **437212** and W=0
A10 **437771** were also RUNNING. No polygon-moment result is available yet.

### Baseline endpoint reached (2026-10-28 02:42 KST)

At 02:42:59 KST, job **438003** had reached its exact target `t=1.1`, step
14,534, with the final RK4 energy record reporting
`dEtot/|Etot0|=-4.395e-12`, `dEhyd/Ehyd0=2.013e-5`, `n_ie_le0=0`, and
`floor_inj=0`. Slurm still showed 438003 RUNNING while its endpoint audit
scripts ran; 438004, 438007, and the dependent polygon audit 438014 remained
held. This is a short baseline endpoint, not a t>6 physics result.

### Baseline endpoint audits passed (2026-10-28 02:44 KST)

The completed **438003** endpoint and radial audits both passed. At `t=1.1`,
the independently recomputed total-energy relative change is
`-6.63e-11`; the W0 and exact fixed-mass checks pass, the disk mass-weighted
density change is `4.78e-4`, entropy-log change `1.96e-4`, position error
`2.08e-8`, and angular-momentum change `-8.53e-13`. The GIZMO annulus velocity
L1 is `1.26e-8`; its literal density L1 is `1.00e-2` because of the floor
plateau, while the plateau-corrected value is `6.29e-4`. All listed endpoint
gates passed, but these cover only `t=1.1`. Slurm had started **438004** at
02:44:25 KST; **438007** and the polygon audit **438014** still waited on the
dependency chain.

### Eta=0 comparison progress (2026-10-28 02:47 KST)

At 02:47:29 KST, **438004** was RUNNING on `syn03`, elapsed 4:27, at
`t=1.02238`, step 13,285. With the Kepler timestep limiter disabled it was
accepting `dt=2.91689e-4`; the latest energy record had
`dEtot/|Etot0|=-9.535e-10`, `n_ie_le0=0`, and zero floor injection. The target
remains the matched `t=1.1` checkpoint; the endpoint audit and the following
OMP candidate **438007** have not completed. Polygon audit **438014** remains
afterany-held behind 438007.

### Eta=0 live update (2026-10-28 02:51 KST)

At 02:51:22 KST, **438004** remained RUNNING on `syn03`, elapsed 8:20, at
`t=1.04231`, step 13,360, with accepted `dt=2.68728e-4`. The latest energy
record reported `dEtot/|Etot0|=-1.620e-9`, `dEhyd/Ehyd0=8.229e-6`,
`n_ie_le0=0`, and zero floor injection. The exact `t=1.1` endpoint remains
ahead; **438007** and **438014** remain dependency-pending.

### Eta=0 progress and corrected performance gate (2026-10-28 02:56 KST)

At 02:56 KST, **438004** was still RUNNING on `syn03`, elapsed 13:06, at
`t=1.06948`, step 13,458. The latest energy record reported
`dEtot/|Etot0|=-2.444e-9`, `n_ie_le0=0`, and zero floor injection. Jobs
**438007** and **438014** remained dependency-pending.

The matched eta/dt/OMP candidate comparison now gates performance on the full
`[RUN_TIME] wall=...` interval around `mpirun`, parsed from the unique Slurm
stdout file for each run. The instrumented `[CPU step]` phase sum remains a
diagnostic only: it omits per-particle RK4 arithmetic and therefore cannot
stand in for complete run time. The candidate must reduce full `mpirun` wall
time to pass the speed gate. This corrects the earlier partial-timer gate;
the candidate pair comparison is still pending completion of job **438007**.

### Eta=0 endpoint audit (2026-10-28 03:03 KST)

Job **438004** completed with exit code 0 after 18:50 wall-clock (the timed
`mpirun` interval was 1126.546 s). It reached the exact `t=1.1` endpoint at
step 13,571. Its audit passed fixed mass, exact W0, positive volume/internal
energy, area tiling, density consistency, and the configured disk/GIZMO
endpoint gates. Total-energy change was `-3.7473e-9`; within the initial disk
`0.5 <= r <= 2`, mass-weighted density change was `4.7820e-4`, entropy-log
change `1.8082e-4`, orbit-position error `2.0847e-8`, and angular-momentum
change `-8.593e-13`. The GIZMO annulus velocity L1 was `1.256e-8`; plateau-
corrected density L1 was `6.293e-4`. No floor or dual-energy injection
occurred. These remain short-time checks only.

The radial audit also exposes a serious outer low-density/edge issue hidden by
the disk-averaged gates: cells initially at `2.75 <= r <= 5` (50,128 cells,
initial mass 0.1303) show mass-weighted density change `0.1488`, radial-speed
change `0.0945`, and angular-momentum change `-0.3159`; the current edge subset
has angular-momentum change `-0.5819`. This is outside the main disk annulus,
so it does not negate the quoted disk-interior metrics, but it must be tracked
and explained before interpreting long-time global behavior. The eta=0 run is
not a matched performance optimization: its 1126.546 s `mpirun` time is only
about 3.67x shorter than the eta=0.06 baseline's 4139.427 s, with a different
timestep policy. Job **438007** is now running for the matched eta=0.06 OMP
candidate comparison; polygon audit **438014** remains dependency-held behind
it.

### Outer-edge pair interpretation correction (2026-10-28)

Directly comparing the completed eta=0 endpoint receipt with the eta=0.06
baseline **438003** shows that the large outer/edge changes are effectively
identical: for initial `2.75 <= r <= 5`, density change is
`0.1487975848` versus `0.1487976109`, radial-speed change is
`0.0945012273` versus `0.0945012299`, and angular-momentum change is
`-0.3159031495` versus `-0.3159031449`. The current-edge subset angular-
momentum changes are `-0.5819034148` versus `-0.5819034014`. Thus this
feature is not caused by setting eta to zero; it is common to the matched
initial setup and is consistent with the already documented periodic seam /
low-density outer-background sensitivity. Causality still requires a
boundary-matched control. The disk density, position, angular momentum and
GIZMO annulus norms agree to the shown precision, but the absolute
entropy-log change is `1.9572e-4` at eta=0.06 versus `1.8082e-4` at eta=0,
and total-energy relative change is `-6.63e-11` versus `-3.75e-9`. Both runs
remain within the current short-time gates and have zero floor injection;
these small differences show the endpoint is not bitwise invariant to the
timestep policy. Do not treat a simple rotation-to-rest potential taper as
the remedy: the existing [seam-control analysis](KEPLER_PERIODIC_SEAM_CONTROL_2026-10-26.md)
shows such a transition is Rayleigh-unstable for this constant-pressure
inviscid disk. A paper-matched isolated/free-surface branch needs its own
bounded-cell force, work, mass, energy, and angular-momentum boundary closure;
see the [boundary-gap audit](KEPLER_ISOLATED_BOUNDARY_GAP_2026-10-27.md).

### Existing long-horizon artifacts checked

The archived run
`/gpfs/kjhan/LagForce/Kepler_res256_t120_tcommon_ringIC_rho0p7_1d0ba7e`
contains an actual trajectory through `t=65.3418` (step 7,887); its
`khout.005872.dat` endpoint is at `t=62.000547`. It used the older
`schemeH_1d0ba7e` executable, the 48,254-cell `glass_kepring_384` ring IC,
`EUNHA_VORO_RHO=0.7`, eight OpenMP threads, and unset `EUNHA_KEPLER_DTETA`.
The run was Slurm-cancelled while still progressing, not at the known t=6–10
wall. At its last logged step the ledger had `dEtot/|Etot0|=3.191e-6` and
`dEhyd/Ehyd0=7.111e-3`; this is historical evidence, not a current-build
validation. It is not an interchangeable restart for the current 147,456-cell
polar q1 setup.

The archived `KepQdump*_t10` files are at `t=10.02605` but contain 12,064
cells (192² glass) and were produced by separate Q2/Q3 reconstruction
diagnostic builds. They are useful for past geometry/torque inspection, not as
matched t=10 endpoints for the current q1 trajectory. Keep these prior
artifacts in view while avoiding a false direct comparison across ICs,
resolutions, or solver builds.

### Current long-run queue checkpoint (2026-10-28 03:10 KST)

The existing full-periodic q1 H200 continuation **437212** is still RUNNING at
`t=1.62388`, step 21,620, after 8:50:12; its latest RK4 ledger has
`dEtot/|Etot0|=-1.101e-10`, `dEhyd/Ehyd0=1.950e-4`, `n_ie_le0=0`, and no
floor injection. The q1 `t=2.221→6.5` segment **437213** remains correctly
dependency-held. The full-periodic q2 continuation **437250** is still PENDING
on `QOSGrpGRES`, with **437251** held behind it. The fresh OMP/cache comparison
**438007** is RUNNING on A10 from the exact q1 t=1.0000125 checkpoint and has
reached `t=1.01296`; its matched endpoint comparison is not yet available.
This profile pair and the dependent one-step moment audit **438014** are
supporting runtime/geometry checks; the already staged q1/q2 long continuations
remain the path for crossing the historical `t≈6` interval.

At 03:12 KST, **438007** was still RUNNING on `syn01` at elapsed 11:43,
`t=1.01702`, step 13,429. The latest ledger reported
`dEtot/|Etot0|=6.380e-12`, `n_ie_le0=0`, and no floor injection. The q1 H200
continuation **437212** remained RUNNING; q2 **437250** remained `QOSGrpGRES`
pending, with no change to its protected dependency chain.

At 03:17 KST, **438007** remained RUNNING at elapsed 15:38, `t=1.02276`,
step 13,506. Its latest energy record was `dEtot/|Etot0|=7.124e-13`, with
`n_ie_le0=0` and zero floor injection. The matched endpoint audit and
dependent one-step GPU moment audit **438014** are still pending.

The q1 two-day recovery continuation **438008** is also staged with
`afterany:437213` and targets `t=6.5`. Its `prepare_resume.py` hashes the
accepted executable/parameters/reference, resumes only after a scheduler
timeout/preemption/node failure, and selects the latest complete positive
fixed-mass W0 checkpoint in `[2.221,6.5]`. It becomes a no-op if the parent
already passed its exact `t=6.5` audit and refuses physics-failed parent
states. It is a guarded scheduler recovery, not a duplicate t=6.5 evolution.

At 03:22 KST, **438007** remained RUNNING at elapsed 20:16, `t=1.02953`,
step 13,598, with `dEtot/|Etot0|=-4.378e-13`, no nonpositive internal energy,
and zero floor injection. q1 H200 **437212** remained RUNNING; its latest
available log reached `t=1.63778`, step 21,807. q2 **437250** remained pending
on `QOSGrpGRES`.

At 03:25 KST, **438007** remained RUNNING at elapsed 24:11, `t=1.03468`,
step 13,669. The latest ledger reported `dEtot/|Etot0|=-4.545e-12`,
`n_ie_le0=0`, and zero floor injection. q1 H200 **437212** had advanced to
`t=1.64249`, step 21,872, with `dEtot/|Etot0|=-1.139e-10`,
`dEhyd/Ehyd0=2.015e-4`, and no nonpositive internal energy or floor
injection. q2 **437250** and the matched OMP=2 H200 probe **437835** remained
QoS-pending; their dependency-protected follow-ons remain held.

## Verified syntax queue and progress (2026-10-28 03:28 KST)

This Codex shell is on `grammar`, whose local Slurm accounting is a different
cluster view. The Kepler run IDs in this document were therefore rechecked
through the configured `syntax` SSH alias (`syntax.kias.re.kr`); that host's
Slurm records confirm they belong to `kjhan` and remain active. Use `syntax`
for subsequent polling of these runs rather than the local `grammar` Slurm
view.

At 03:28 KST, A10 candidate **438007** was RUNNING on `syn01`, elapsed 26:21.
Its latest log marker reached `t=1.03852`, step 13,720; the RK4 ledger gave
`dEtot/|Etot0|=-4.326e-12`, `dEhyd/Ehyd0=7.478e-6`, `n_ie_le0=0`, and zero
floor injection. Its `t=1.1` matched-endpoint audit has not run yet. Baseline
**438003** and eta=0 run **438004** both completed successfully; the GPU
polygon-moment audit **438014** remains dependency-held after 438007.

Full-periodic q1 H200 **437212** was RUNNING on `syn104`, elapsed 9:08:05, at
`t=1.64566`, step 21,916. Its latest ledger reported
`dEtot/|Etot0|=-1.150e-10`, `dEhyd/Ehyd0=2.026e-4`, `n_ie_le0=0`, and zero
floor injection. The `t=2.221` continuation **437213** and guarded `t=6.5`
recovery **438008** are still dependency-held.

The 420² refined W=0 A10 continuation **437771** was RUNNING on `syn02`,
elapsed 8:52:18, at `t=0.863622`, step 11,364. Its latest ledger reported
`dEtot/|Etot0|=-4.748e-11`, `dEhyd/Ehyd0=6.682e-5`, `n_ie_le0=0`, and zero
floor injection. It has not reached its `t=2.221` endpoint. q2 full-periodic
H200 segment **437250** and matched OMP=2 probe **437835** remain pending for
`QOSGrpGRES`; downstream jobs **437251** and **437836** remain dependency-held.

At 03:31 KST, candidate **438007** remained RUNNING at `t=1.04227`, step
13,769, elapsed 29:32. Its latest energy ledger had
`dEtot/|Etot0|=-4.163e-12`, `dEhyd/Ehyd0=8.223e-6`, and no nonpositive
internal energy or floor injection; the matched `t=1.1` endpoint and full
runtime are still pending, so no speedup conclusion is available yet.

The q1 H200 continuation **437212** remained RUNNING at `t=1.64882`, step
21,960, elapsed 9:11:09. Its latest ledger was
`dEtot/|Etot0|=-1.194e-10`, `dEhyd/Ehyd0=2.037e-4`, with `n_ie_le0=0` and
zero floor injection. The q2 H200 segment **437250** became RUNNING on
`syn104` at 03:30:29 after the QoS wait; by this check it had advanced from
its audited `t=2.22138` restart to `t=2.22199`, step 30,084. Its latest
ledger was `dEtot/|Etot0|=-1.970e-13`, `dEhyd/Ehyd0=5.175e-7`, with no
nonpositive internal energy or floor injection. This is only the segment's
startup interval, not a long-horizon result. The independent 420² refined W=0
A10 branch **437771** was at `t=0.867519`, step 11,415, with
`dEtot/|Etot0|=-4.300e-11`, `dEhyd/Ehyd0=6.745e-5`, `n_ie_le0=0`, and zero
floor injection. The downstream q1/q2 continuations and the one-step audit
remain dependency-held; matched OMP=2 remains QoS-held.

At 03:33 KST, candidate **438007** was RUNNING at `t=1.04553`, step 13,812,
elapsed 31:02; its latest ledger had `dEtot/|Etot0|=-4.727e-12`,
`dEhyd/Ehyd0=8.875e-6`, `n_ie_le0=0`, and zero floor injection. The q1 H200
run **437212** reached `t=1.65180`, step 22,001, with
`dEtot/|Etot0|=-1.228e-10`, `dEhyd/Ehyd0=2.048e-4`, and no positivity/floor
events. q2 H200 **437250** advanced to `t=2.22465`, step 30,123, with
`dEtot/|Etot0|=-7.512e-13`, `dEhyd/Ehyd0=1.964e-6`, and no such events. The
420² W=0 A10 branch **437771** reached `t=0.870947`, step 11,460; its ledger
was `dEtot/|Etot0|=-4.727e-11`, `dEhyd/Ehyd0=6.801e-5`, with no such events.
All four remain short of their respective endpoint gates; this checkpoint
does not establish long-time stability or a performance win.

The face-speed implementation was checked against the production source:
W=0 automatically uses the Voronoi face-centroid correction on CPU and GPU,
even when `GFS_LAGUERRE_ROTATION` is unset. That environment switch extends
the geometry correction to nonzero-weight Laguerre faces. Neither path
evolves intrinsic cell spin; see the [spin integration gates](KEPLER_SPIN_NATIVE_INTEGRATION_GATES_2026-10-27.md).

At 03:37 KST, the matched performance candidate **438007** was still RUNNING
at `t=1.05019`, step 13,873, elapsed 35:12. Its latest RK4 ledger had
`dEtot/|Etot0|=-6.028e-12`, `dEhyd/Ehyd0=9.812e-6`, `n_ie_le0=0`, and zero
floor injection. The exact `t=1.1` endpoint and full wall time remain pending.
q1 H200 **437212** was at `t=1.65607`, step 22,059, with
`dEtot/|Etot0|=-1.244e-10`, `dEhyd/Ehyd0=2.063e-4`, and no positivity/floor
events. q2 H200 **437250** was at `t=2.22810`, step 30,177, with
`dEtot/|Etot0|=4.713e-13`, `dEhyd/Ehyd0=3.845e-6`, and no such events. The
420² W=0 A10 run **437771** was at `t=0.875761`, step 11,523, with
`dEtot/|Etot0|=-4.262e-11`, `dEhyd/Ehyd0=6.879e-5`, and no such events.
All four remained below their endpoint gates. The benchmark-duration audit
now distinguishes the plotted GIZMO `t=600` image from the paper's separate
approximate 1,000-orbit statement; see the [boundary and benchmark audit](KEPLER_ISOLATED_BOUNDARY_GAP_2026-10-27.md).

At 03:41 KST, q2 H200 **437250** provided a clear follow-up to the small
timestep seen at step 143 (`dt=5.230e-6`): by local step 178, `dt` had
recovered to `7.273e-5`, and at step 179 it was back to `7.642e-5`. The run
reached `t=2.23300`, step 30,249; its latest ledger showed
`dEtot/|Etot0|=1.345e-12`, `dEhyd/Ehyd0=6.532e-6`, `n_ie_le0=0`, and zero
floor injection. The episode is transient in the observed interval, not a
sustained stall. `GFS_VARPRESS_GEOM_AUDIT` is disabled in this runner, so the
limiting cell remains unidentified; keep the existing coefficient unchanged.

At the same 03:41 poll, performance candidate **438007** was RUNNING at
`t=1.05630`, step 13,953, with `dEtot/|Etot0|=-7.254e-12`,
`dEhyd/Ehyd0=1.105e-5`, and no positivity/floor events. q1 H200 **437212** was
at `t=1.66106`, step 22,134, with `dEtot/|Etot0|=-1.287e-10` and
`dEhyd/Ehyd0=2.081e-4`; the 420² W=0 A10 branch **437771** was at `t=0.881916`,
step 11,604, with `dEtot/|Etot0|=-4.778e-11` and
`dEhyd/Ehyd0=6.980e-5`. All three had zero nonpositive internal-energy and
floor-injection counts. No endpoint or matched performance comparison is
available yet.

### Orbital timestep-cap controlled pair (2026-10-28 03:44 KST)

Completed A10 job **438035** compared W=0 q1 restarts from the same `t=0.825`
checkpoint and executable, changing only `EUNHA_KEPLER_DTETA` from `0.06` to
`0`. Both reached the exact `t=0.845` endpoint with 147,456 cells and
bitwise-identical per-cell masses. With the limiter enabled, the segment used
263 steps, ended with `dt=4.5008e-5`, and took 1,101.990 s. With it disabled,
the segment used 64 steps, ended with `dt=2.0524e-4`, and took 271.563 s:
4.56x larger endpoint timestep and 4.06x lower measured wall time on the same
A10 allocation.

The disk-only errors over cells active at the restart are effectively the
same at this short endpoint: RMS relative radius error `6.3776416e-10` vs
`6.3776367e-10`, RMS `|vr|/vc` `4.0267945e-8` vs `4.0267806e-8`, mean relative
azimuthal-velocity error `5.7743141e-9` vs `5.7743065e-9`, and active-region
relative `Lz` changes `8.42437e-13` vs `8.41993e-13` (enabled vs disabled).
Relative total-energy change from the restart was `1.58e-12` enabled and
`9.40e-10` disabled. Both are small over this interval, but the uncapped run's
energy error is about 595x larger. These results support a short local
accuracy/performance comparison only; they do not establish safe uncapped
evolution through the historical `t=6-10` failure window. Keep the cap in the
long branches until a matched longer-horizon physics gate passes.

### Live continuation checkpoint (2026-10-28 03:47 KST)

The scheduler still showed **438007**, **437212**, **437250**, and **437771**
RUNNING; their endpoint audit jobs remain dependency-held. A subsequent read of
their shared logs found:

- A10 performance candidate **438007** reached `t=1.06738`, step 14,100;
  latest RK4 ledger: `dEtot/|Etot0|=-5.273e-12`,
  `dEhyd/Ehyd0=1.331e-5`, no nonpositive internal energy or floor injection.
  It has not reached its exact `t=1.1` comparison endpoint.
- q1 H200 **437212** reached `t=1.67068`, step 22,271, with latest ledger
  `dEtot/|Etot0|=-1.282e-10`, `dEhyd/Ehyd0=2.115e-4`, and no positivity/floor
  events. In the latest 177 time markers, 11 fell below half the nominal
  `7.642e-5` step; the minimum `9.36e-7` occurred in a short burst, and the
  step had recovered to `7.642e-5` by step 22,275. The runner has no geometric
  limiter-cause audit, so the transient's cause is unknown; this sample does
  not justify changing the CFL coefficient.
- q2 H200 **437250** reached `t=2.24194`, step 30,377; latest ledger
  `dEtot/|Etot0|=-2.540e-12`, `dEhyd/Ehyd0=1.145e-5`, with no positivity/floor
  events.
- Refined 420² W=0 A10 **437771** reached `t=0.892957`, step 11,752; latest
  ledger `dEtot/|Etot0|=-5.012e-11`, `dEhyd/Ehyd0=7.164e-5`, with no
  positivity/floor events.

The short `DTETA=0` pair is encouraging for local disk accuracy and runtime,
but the live long branches retain `DTETA=0.06`; none of these checkpoints yet
reaches the historical `t=6-10` window.

### Follow-up log read (2026-10-28 03:51 KST)

The 03:50:26 scheduler poll still found **438007**, **437212**, **437250**, and
**437771** RUNNING; **438014**, **438008**, **437213**, **437251**, and the
matched OMP probes remain held by dependencies or the H200 QoS limit. The
next log read showed:

- **438007** at `t=1.07098`, step 14,147; latest ledger
  `dEtot/|Etot0|=-8.179e-12`, `dEhyd/Ehyd0=1.405e-5`, no positivity/floor
  events. The `t=1.1` matched audit is still pending.
- q1 **437212** at `t=1.67389`, step 22,313, with `dt=7.642e-5` again after
  the brief small-step burst documented above; latest ledger
  `dEtot/|Etot0|=-1.284e-10`, `dEhyd/Ehyd0=2.127e-4`, no positivity/floor
  events.
- q2 **437250** at `t=2.24476`, step 30,418; latest
  `dEtot/|Etot0|=-1.767e-12`, `dEhyd/Ehyd0=1.301e-5`, no positivity/floor
  events.
- Refined 420² W=0 **437771** at `t=0.896354`, step 11,798; latest
  `dEtot/|Etot0|=-4.975e-11`, `dEhyd/Ehyd0=7.221e-5`, no positivity/floor
  events.

No exact endpoint or long-horizon conclusion was available at this checkpoint.

### q1 disk-versus-ambient internal-energy audit (snapshots through t=1.5)

Compared the fixed-ID cohorts in q1 snapshots at `t=0.825`, `1.0000125`, and
`1.5000267`. The effective `Gamma_hydro` in this run is `5/3`; cell masses
remain bitwise identical by ID. Cohorts are selected by their radius at
`t=0.825`, so the disk and exterior are followed Lagrangianly rather than
re-binned by their later positions.

For the disk cohort initially at `0.5 <= r < 2` (48,819 cells), summed
internal energy relative to `t=0.825` was `0.99999618` at `t=1.0` and
`0.99998346` at `t=1.5`. At `t=1.5`, its RMS fractional radius drift was
`2.35e-8` and its mass-weighted mean entropy-proxy change
`d ln(P/rho^(5/3))` was `-1.41e-5`. These snapshots show no measurable net
inner-disk heating through `t=1.5`.

The global sum of cell internal energies rose by `0.519%` at `t=1.0` and
`2.465%` at `t=1.5`, but the increase is almost entirely in the cohort
initially outside `r=2.5`. That exterior cohort's internal-energy sum rose
`1.126%` and `5.347%`, while at `t=1.5` its mass-weighted changes were
`d ln rho=0.0605701`, `d ln P=0.1009504`, and `d ln u=0.0403803`. These match
the adiabatic compression relations for `Gamma=5/3`;
`d ln(P/rho^(5/3))=1.93e-7`. The global internal-energy growth therefore
does not indicate shear heating of the Kepler disk in these snapshots; it is
dominated by nearly isentropic compression in the exterior gas. The cause of
that exterior compression remains to be examined, and this early interval
does not establish stability through `t=6-10`.

The repeatable audit is `Docs/audit_kepler_lagrangian_entropy.py`; its
snapshot receipts for `t=1.0` and `t=1.5` are in
`Docs/KEPLER_Q1_LAGRANGIAN_ENTROPY_AUDIT_T15_2026-10-28.json`. It verifies
fixed mass by cell ID before computing these cohort metrics.

### Live continuation checkpoint (2026-10-28 04:01 KST)

At 04:01:07 KST, the `syntax` scheduler still showed **438007**, **437212**,
**437250**, and **437771** RUNNING. A subsequent shared-log read found:

- **438007** at `t=1.08647`, step 14,352; latest ledger
  `dEtot/|Etot0|=-7.504e-12`, `dEhyd/Ehyd0=1.727e-5`, no positivity/floor
  events. Its exact `t=1.1` endpoint comparison is still pending.
- q1 **437212** at `t=1.68710`, step 22,495; latest ledger
  `dEtot/|Etot0|=-1.356e-10`, `dEhyd/Ehyd0=2.174e-4`, no positivity/floor
  events. The latest accepted step was `5.55e-5`, with previous short dips
  having recovered to nominal values.
- q2 **437250** at `t=2.25713`, step 30,595; latest
  `dEtot/|Etot0|=1.947e-13`, `dEhyd/Ehyd0=1.988e-5`, no positivity/floor
  events.
- Refined 420² W=0 **437771** at `t=0.911745`, step 12,001; latest
  `dEtot/|Etot0|=-5.422e-11`, `dEhyd/Ehyd0=7.482e-5`, no positivity/floor
  events.

All remain short of their respective endpoint gates; these values do not
establish `t=6-10` stability.

### q1 short-step burst recovery (2026-10-28 04:04 KST)

In the latest 200 q1 **437212** time markers, short-step excursions recurred
around global steps 22,391-22,393, 22,496-22,497, 22,524-22,525, and
22,528-22,533. The minimum was `1.90e-7` at `t=1.68929`; by step 22,534 the
accepted step had returned to `7.642e-5`. The energy ledger remained
conservative with no nonpositive internal energy or floor injection. The
runner still has `GFS_VARPRESS_GEOM_AUDIT` disabled, so the limiting cell and
cause are unknown. This was a brief recovery, not a sustained timestep stall;
leave the accepted geometric-CFL coefficient unchanged.

In the same log read, candidate **438007** had reached `t=1.08937`, step
14,391, with `dEtot/|Etot0|=-8.658e-12`, `dEhyd/Ehyd0=1.788e-5`, and no
positivity/floor events. It had not yet reached the `t=1.1` audit endpoint.

The archived native overlay confirms that a targeted diagnostic run with
`GFS_VARPRESS_GEOM_AUDIT=1` reports each RK4 stage's swept-volume `rate`,
`hSoverV`, and limiting cell ID in `[RK4_GEOM_CFL]`. The active runner has
this disabled. Reserve it for a matched diagnostic checkpoint if small steps
become persistent; the present bursts recovered without changing the physics
setting.

### Matched OMP/cache candidate audit (438007, t=1.1)

The matched baseline **438003** and candidate **438007** both completed with
exit code 0 at exact `t=1.1`. Global and radial physics receipts passed, step
sample counts matched, and all ten audited endpoint fields
(`x,y,vx,vy,mass,volume,ie,den,pressure,w2`) were bitwise identical. The full
`mpirun` wall times were `4139.427 s` baseline and `4190.382 s` candidate; the
candidate was 50.955 s (1.23%) slower in this pair. Its median GPU-stage
timer was lower (`172.8 ms` vs `197.1 ms`), but the median CPU-step total was
unchanged at `2.4 s`, and total runtime did not improve. Do not promote this
OMP/cache candidate as an acceleration based on this evidence. The result is
an exact short-horizon equivalence gate, not a long-time stability claim.
The printed CPU subphases (rounded to `0.1 s`) shifted from update/force/final
medians `1.2/0.8/0.4 s` to `1.1/0.7/0.6 s`; total median remained `2.4 s`.
The runs used different A10 nodes, so the 1.23% wall-time difference does not
prove a robust regression either. It does show no measured end-to-end gain at
OMP=2; use the already queued OMP-scaling probes for the next performance
decision instead of promoting this patch.
Full receipt: `Docs/KEPLER_OMP_CACHE_MATCHED_AUDIT_20261028.json`.

### GPU polygon area/centroid/Q audit (438014)

The one-step A10 CUDA audit completed with exit code 0 from the accepted q1
restart at `t=1.0000125425` to `t=1.0000625425`. It examined 147,456 cells;
all were valid and the maximum normalized GPU-versus-host-reference moment
error was `8.397862e-16`, against tolerance `1e-8`. This validates the
polygon area, centroid, and central-Q calculation for this W=0 geometry
checkpoint, providing a geometric prerequisite for spin work. It does not add
or validate a native spin evolution state.

### Live long-branch checkpoint after short gates (2026-10-28 04:14 KST)

At the 04:14:17 KST scheduler poll, q1 **437212**, q2 **437250**, and refined
420² W=0 **437771** remained RUNNING. Their `t=6.5` continuations and
matched OMP probes were still dependency/QoS-held. Shared logs then showed:

- q1 at `t=1.70164`, step 22,705; `dEtot/|Etot0|=-1.367e-10`,
  `dEhyd/Ehyd0=2.227e-4`, with no nonpositive internal energy or floor
  injection.
- q2 at `t=2.27035`, step 30,793; latest ledger
  `dEtot/|Etot0|=-9.730e-13`, `dEhyd/Ehyd0=2.729e-5`, with no positivity/floor
  events. Its accepted `dt` was below nominal (`3.8876e-5`) but not yet shown
  to be a sustained restriction.
- refined 420² W=0 at `t=0.929111`, step 12,234; latest ledger
  `dEtot/|Etot0|=-6.025e-11`, `dEhyd/Ehyd0=7.783e-5`, no positivity/floor
  events.

The short matched OMP/cache and one-step polygon-moment gates are complete;
none of these long branches has reached its planned endpoint or the historical
`t=6-10` window.

### Archived `rho=0.7` Lagrangian cohort heating history (2026-10-28 04:18 KST)

The archived 48,254-cell `schemeH_1d0ba7e` trajectory was audited by fixed
cell ID, grouping cells by their initial radius at `t=2.02621`. The cell masses
are bitwise unchanged by ID. This is a radial localization diagnostic for that
older run, not a matched validation of the current 147,456-cell q1 build.

| Target t | Whole-disk internal-energy sum / initial | Initial r=0.5–1 cohort | Initial r=1–2 cohort | Inner RMS fractional radial drift |
|---:|---:|---:|---:|---:|
| 4.016 | 1.00195 | 1.01089 | 0.99969 | 7.29e-6 |
| 6.004 | 1.00404 | 1.02253 | 0.99936 | 1.43e-5 |
| 8.029 | 1.00633 | 1.03477 | 0.99915 | 2.06e-5 |
| 10.017 | 1.00879 | 1.04704 | 0.99913 | 2.71e-5 |
| 14.026 | 1.01463 | 1.07444 | 0.99952 | 5.76e-5 |
| 26.006 | 1.09638 | 1.48298 | 0.99871 | 6.00e-4 |
| 40.002 | 1.89918 | 5.46334 | 0.99865 | 2.53e-3 |
| 62.001 | 490.332 | 2426.91 | 1.03317 | 7.21e-2 |

The mass-weighted mean entropy-proxy change `ln(P/rho^gamma)` in the inner
cohort grows from `0.0134` at `t=6.004`, to `0.0232` at `t=10.017`, `0.0986`
at `t=26.006`, `0.2746` at `t=40.002`, and `3.412` at `t=62.001`. The outer
cohort stays near its initial internal energy through `t=40`; by `t=62` its
energy sum is only 3.3% higher. This shows that the late catastrophic energy
growth is strongly concentrated in the inner ring of this old trajectory.
It does not identify the cause (e.g. artificial-viscosity heating, local
resolution/noise, or another inner-disk instability), and its smooth but small
inner heating from `t=6` to `t=10` differs from the separate previously
reported abrupt `t=6-10` failure. Use current q1/q2 radial cohort and shock/
viscous-work diagnostics at their planned endpoints to establish whether the
same localization and onset occur in the present setup.

Audit: `KEPLER_OLD_RHO07_T2_T62_COHORT_HISTORY_2026-10-28.json`. The single
t2-to-t62 aggregate receipt is `KEPLER_OLD_RHO07_T2_TO_T62_LAGRANGIAN_AUDIT_2026-10-28.json`.

### Live q1/q2 and no-cap checkpoint (2026-10-28 04:20 KST)

The `syntax` scheduler showed q1 H200 **437212**, q2 H200 **437250**, refined
420² W=0 A10 **437771**, and matched no-cap A10 **438040** RUNNING. The q1
and q2 continuations toward their longer gates remain held behind their
respective parent jobs; no job was restarted or duplicated.

- q1 **437212** reached `t=1.70912`, step 22,808; the latest ledger had
  `dEtot/|Etot0|=-1.362e-10`, `dEhyd/Ehyd0=2.255e-4`, and no nonpositive
  internal energy or floor injection.
- q2 **437250** reached `t=2.27746`, step 30,898; the latest ledger had
  `dEtot/|Etot0|=-3.049e-12`, `dEhyd/Ehyd0=3.131e-5`, and no positivity/floor
  events. The last accepted `dt=1.0443e-5` followed two short steps; this
  sample alone does not establish a sustained timestep restriction.
- Refined **437771** reached `t=0.937859`, step 12,358; its latest ledger had
  `dEtot/|Etot0|=-5.565e-11`, `dEhyd/Ehyd0=7.936e-5`, and no positivity/floor
  events.
- The new matched comparison **438040** restarted from the eta=0 branch's
  common `t=0.845` snapshot and reached `t=0.935530`, global step 11,236. Its
  latest ledger had `dEtot/|Etot0|=-3.919e-9`, `dEhyd/Ehyd0=1.503e-5`, and no
  positivity/floor events. Its run script changes only
  `EUNHA_KEPLER_DTETA=0`; in the source, positive eta imposes a central-orbit
  timestep cap and eta=0 disables it. The W=0, fixed-mass 384² IC is otherwise
  retained. The exact t=1.5 gate compares active-annulus orbit errors, angular
  momentum, energy, and radial density to the capped H200 reference. It has
  not reached the endpoint, so the timestep-cap effect remains unknown.

The existing paired probe at `t=0.845` found exact per-ID mass equality and
nearly identical active-annulus state metrics for eta=0.06 and eta=0: RMS
fractional radius error was `6.378e-10` in both, and RMS `vr/vc` was
`4.027e-8` in both. Across the subsequent common-time log interval
`t=0.845-0.9355`, the capped reference accepted 1,207 steps (median
`dt=7.642e-5`), while the no-cap run accepted 316 (median `dt=3.049e-4`),
about 3.8 times fewer steps over nearly the same elapsed simulation time.
This is a step-count reduction, not yet a wall-clock speedup or an accuracy
result. The runs use different GPUs (H200 reference, A10 candidate); the
planned exact t=1.5 physical audit is still required before adopting eta=0.

These checkpoints are all before the historical `t=6-10` interval and do not
validate long-time stability or native spin evolution.

### No-cap energy tradeoff and native-spin code audit (2026-10-28 04:22 KST)

At nearly matched times around `t=0.935`, the capped H200 log at `t=0.934990`
reported `dEtot/|Etot0|=-2.012e-11`, `dEhyd/Ehyd0=1.808e-5`, and
`E_int=5.479941e-5`. The no-cap A10 log at `t=0.935094` reported
`-3.899e-9`, `1.495e-5`, and `5.479963e-5`, respectively. Both had positive
internal energy and zero floor injection. Together with the 3.8x smaller
accepted-step count over `t=0.845-0.9355`, this is a real efficiency-versus-
energy-ledger tradeoff: the no-cap total-energy drift is still small in
absolute terms, but about 190x larger than the capped branch at this early
checkpoint. Do not accept `eta=0` based on step count alone; retain the
planned exact `t=1.5` orbit, density, angular-momentum, and energy audit.

The two saved `t=0.845` branch states used by the short probe are not bitwise
identical: IDs, fixed masses, and `w2=0` match exactly, while maximum absolute
position differences are `9.40e-7`, velocity-component differences `3.44e-4`,
volume differences `6.30e-9`, and internal-energy differences `1.61e-11`.
These small differences arise in the preceding eta-dependent evolution from
the shared `t=0.825` restart. The planned `t=1.5` audit uses that common
`t=0.825` initial state and therefore measures the full timestep-policy
divergence; it is not a counterfactual that switches eta at an identical
`t=0.845` state. A matched continuation from one exact `t=0.845` checkpoint
would only be needed if the full-branch comparison leaves the cap's isolated
effect ambiguous.

At the 04:22:44 KST scheduler poll, q1 **437212**, q2 **437250**, refined
420² W=0 **437771**, and no-cap **438040** were all still RUNNING. Their
latest log times were `1.71225`, `2.27978`, `0.940611`, and `0.943666`,
respectively. Energy ledgers remained free of nonpositive internal energy and
floor injection. The no-cap run had not reached its `t=1.5` audit.

The current worktree's `voronoi_face_rotation` is a face-kinematics
correction: it adjusts normal face speed using neighboring generator
velocities and the centroid's offset from the Voronoi midpoint. The W=0 path
is applied regardless of `GFS_LAGUERRE_ROTATION`; that flag only extends the
correction to nonzero Laguerre weights. This is not a stored cell-spin state.
The particle/RK4 structures have no intrinsic-spin or spin-energy fields;
`w2` is the Laguerre geometry weight. Current q1/q2 scripts are W=0 and unset
the weighted-face extension flag. Existing isolated affine-spin and
two-cell-ledger prototypes are documented in
`KEPLER_AFFINE_SPIN_RECONSTRUCTION_2026-10-27.md`,
`KEPLER_SPIN_CENTROID_CLOSURE_2026-10-27.md`, and
`KEPLER_VARIATIONAL_SPIN_PAIR_TORQUE_2026-10-27/README.md`. They are not
connected to the native RK4 solver; centroid evolution, a physically specified
face torque, and non-rigid kinetic-energy closure remain the native-integration
gaps.

At the 04:26:07 KST poll, **437212**, **437250**, **437771**, and **438040**
were still RUNNING. Latest times were q1 `t=1.71631`, q2 `t=2.28388`, refined
420² `t=0.94557`, and no-cap `t=0.954500`. The corresponding relative total-
energy ledgers were `-1.388e-10`, `-5.380e-12`, `-6.031e-11`, and
`-4.600e-9`; all still reported zero nonpositive internal-energy cells and
zero floor injection. The no-cap candidate had not yet produced its exact
`t=1.0` checkpoint or reached the planned `t=1.5` comparison gate.

At the 04:28:18 KST poll the same four jobs remained RUNNING. Their latest
times were q1 `t=1.71893` (`dEtot/|Etot0|=-1.364e-10`), q2 `t=2.28649`
(`-6.381e-12`), refined 420² `t=0.948286` (`-6.014e-11`), and no-cap
`t=0.963370` (`-4.948e-9`). Each latest record still had `n_ie_le0=0` and
`floor_inj=0`; the no-cap `t=1.0` snapshot was not yet present.

At the 04:30:49 KST poll all four jobs remained RUNNING. Latest times were q1
`t=1.72222` (`dEtot/|Etot0|=-1.341e-10`), q2 `t=2.28932`
(`-7.937e-12`), refined 420² `t=0.951805` (`-6.042e-11`), and no-cap
`t=0.970943` (`-5.223e-9`). The latest ledgers still had no nonpositive
internal-energy cells or floor injection. The no-cap branch had not yet
produced the `t=1.0` snapshot.

### Paired timestep-cap audit at t≈1 (2026-10-28 04:41 KST)

The capped reference snapshot is at `t=1.000012543`; the no-cap snapshot is
at `t=1.000074674` (offset `6.21e-5`). Both contain the same 147,456 cell
IDs, exact per-ID masses, and W0 geometry. Using common initial state
`t=0.825`, the active `0.5≤r≤2` metrics are nearly unchanged between runs:

| Metric | eta=0.06 H200 | eta=0 A10 |
|---|---:|---:|
| RMS fractional radius error | `5.7981e-9` | `5.8000e-9` |
| RMS `vr/vc` | `4.5701e-8` | `4.5687e-8` |
| Mean absolute `vphi` error / `vc` | `7.4172e-9` | `7.4172e-9` |
| Disk relative angular-momentum change | `5.6943e-12` | `5.6626e-12` |
| Mean absolute annular-density change | `5.8090e-6` | `5.6277e-6` |
| Relative total-energy change | `-3.0748e-11` | `-7.1858e-9` |

The fixed-mass Lagrangian cohort audit likewise finds disk internal-energy
sum ratios `0.99999618` (capped) and `0.99999775` (no-cap); the disk
mass-weighted mean entropy-proxy changes are `-2.75e-6` and `-1.17e-6`.
At this early checkpoint the no-cap branch reproduces the disk orbit and
thermodynamic profile to the measured precision, with a still-small but
larger total-energy drift. Its earlier common-time sample used about 3.8x
fewer steps. These snapshots are about one no-cap step apart in time, and the
candidate is on A10 while the reference is on H200; do not interpret this as
a wall-time benchmark or a t>6 stability result. The no-cap job is still
running toward the exact `t=1.5` endpoint audit.

The snapshot parser's metadata key `n` is the particle count, not a global
time-step index. The paired audit scripts and t≈1 JSON receipt now report the
cell count explicitly; this label correction does not change any physics
metrics above.

Reproducible receipts: `KEPLER_DTETA_T1_KINEMATIC_AUDIT_2026-10-28.json`,
`KEPLER_DTETA_CAPPED_VS_NOCAP_T1_COHORT_2026-10-28.json`, and
`audit_kepler_dteta_pair_snapshots.py`.

At the 04:41:38 KST poll the four jobs remained RUNNING: q1 **437212** at
`t=1.73579`, q2 **437250** at `t=2.30177`, refined **437771** at `t=0.966963`,
and no-cap **438040** at `t=1.01261`. The latest no-cap ledger was
`dEtot/|Etot0|=-6.780e-9`, `dEhyd/Ehyd0=2.926e-5`, with no nonpositive
internal energy or floor injection. It continues toward the exact `t=1.5`
endpoint.

At the 04:43:04 KST poll the same four jobs were still RUNNING. Latest times
were q1 `t=1.73759`, q2 `t=2.30356`, refined 420² `t=0.968928`, and no-cap
`t=1.01810`. The no-cap ledger remained free of positivity/floor events, with
`dEtot/|Etot0|=-7.039e-9` and `dEhyd/Ehyd0=3.032e-5`.

At the 04:44:05 KST poll all four remained RUNNING. The latest log markers
were q1 `t=1.73891` (`dEtot/|Etot0|=-1.297e-10`), q2 `t=2.30445`
(`-9.872e-12`), refined 420² `t=0.970539`, and no-cap `t=1.02172`
(`-7.193e-9`). All energy logs still had zero nonpositive internal-energy
cells and zero floor injection. In the refined 420² run, the last 100 accepted
steps contained four short-`dt` steps (minimum `2.22e-6`); the next step had
already recovered to nominal `7.642e-5`. This is a transient, not a sustained
timestep collapse; keep the existing geometric-CFL setting.

At the 04:45:46 KST Slurm poll, jobs **437212**, **437250**, **437771**, and
**438040** were still RUNNING. Fresh log reads around 04:46 KST gave q1
`t=1.74232` (`dEtot/|Etot0|=-1.325e-10`, `n_ie_le0=0`, `floor_inj=0`), q2
`t=2.30763` (`-9.699e-12`, no positivity/floor events), refined 420²
`t=0.97354` (`-6.808e-11`, no positivity/floor events), and the no-cap branch
`t=1.03142` (`-7.578e-9`, no positivity/floor events). The no-cap branch has
not reached its planned `t=1.5` comparison endpoint; keep it running and do
not draw a long-time conclusion from this update.

At the 04:48:31 KST Slurm poll, the same four jobs remained RUNNING. Fresh log
reads gave q1 `t=1.74385` (`dEtot/|Etot0|=-1.337e-10`), q2 `t=2.30974`
(`-1.290e-11`), refined 420² `t=0.976196` (`-6.882e-11`), and no-cap
`t=1.03925` (`-7.846e-9`); all four still reported `n_ie_le0=0` and
`floor_inj=0`. In their last 60 `Time=` records, the minimum `dt` and count
below `1e-5` were respectively q1 `5.45e-8` / 9, q2 `4.99e-6` / 2, refined
420² `5.81e-5` / 0, and no-cap `2.53e-6` / 3. The q1 and q2 latest `dt`
values had recovered to `7.642e-5`; these are recurring brief dips, not a
demonstrated sustained timestep collapse.

The GPU force kernel computes several candidate limits (face signal/CFL,
acceleration, and entropy-change `dK`), then retains only the minimum in
`dt_out`; it does not report which candidate or cell set the global minimum.
The earlier candidate list omitted a fifth limiter: with
`GFS_VARPRESS_GEOM_CFL=.1`, the W=0 kernel also bounds the sum of absolute face
swept-volume rates by `var_geom_cfl / (sum(|Vdot_face|)/V_cell)`. The existing
logs still cannot attribute q1's sub-microsecond steps to a specific limiter.
Preserve all stability limits until reason/argmin telemetry identifies the
active one. The earlier comparison of the central-orbit cap only at the active
annulus edge (`R=0.5`, `dt_eta≈2.12e-2`) was incomplete: `kepler_dt_limit`
takes the minimum over every cell, including low-density cells inside the
annulus.

At the 05:08 KST log read, q1's original run had advanced to `t=1.77123`.
Parsing all 12,825 accepted-step records from its log found 122 records with
`dt<1e-5`; the minimum was `5.45063e-8` at `t=1.74353`. This confirms
recurrent short steps, but not yet their cause.

An opt-in diagnostic was added to the isolated reconstruction of the
authenticated 2026-10-26 native source. `EUNHA_GPU_DT_DIAG=1` records the
winning per-cell limiter (face CFL, viscous CFL, acceleration CFL, entropy
`dK` CFL, or W=0 swept-volume geometry CFL), cell ID/position, and the face
neighbor when applicable. It is disabled by default. Patch script:
`patch_kepler_gpu_dt_reason_diag_20261028.py`; source hash receipt is in the
candidate tree as `dt_reason_diag_patch_receipt.json`. The instrumented binary
SHA-256 is
`b6103c30effc0dcbef927c48a49fac9eebf4b59c0cb5555fdbecb6c17551e29b`.

Diagnostic job **438070** was submitted in a separate directory from the
complete q1 checkpoint `khout.019926.dat` (`step=19926`, `t=1.50003`) through
`t=1.77`. It requests one H200 and keeps the original q1 run and its files
untouched. Its current Slurm state is **PENDING**, reason `QOSGrpGRES`; no
diagnostic evolution has started yet. The immediate gate is to capture the
limiter at the known event near `t=1.74353`, then compare the diagnostic run's
energy/timestep history against the original continuation from the same
checkpoint. No CFL change is authorized by this evidence alone.

At the 05:15 KST read, original q1 had advanced to `t=1.77659` with
`dEtot/|Etot0|=-1.45e-10`, `n_ie_le0=0`, and `floor_inj=0`. Across the last
200 logged steps, median elapsed time was `2.9 s/step`: `update=1.5 s`,
`force=0.9 s`, and `final=0.5 s` (the logged post phase rounded to zero).
This gives a concrete profiling target alongside the short-`dt` diagnosis;
the phase labels alone do not yet establish which work inside `update` or
`force` dominates or can be parallelized further.

### Integrator identity audit (2026-10-28)

The effective q1 `params.dat` says `Hydro time-stepping way = 1 # KDK`, but
the numeric value is authoritative: `Params/params.h` defines `0` as leapfrog
and `1` as RK4, and `Exam/Kepler/kepler.c` dispatches `GAS_EVOLMETHOD==1` to
`exam2d_vph_rk4_int`. The active run log also prints `[RK4E]`. Thus q1 and the
matched diagnostic restart are RK4; the `# KDK` suffix is a stale/misleading
comment and does not select KDK. Do not interpret these runs as evidence about
the KDK integrator or individual time stepping.

### GIZMO long-horizon reference gate

The requested long-time target is much later than `t=6`. In Hopkins (2015),
the Kepler disk is shown at `t=120` in Fig. 6 and `t=600` in Fig. 7; the text
also reports an MFM integration to roughly 1000 inner-disk orbits. The paper's
Fig. 8 convergence test measures the velocity error at 10 orbits at `R≈1`
(about `t=62.83` for `GM=1`). Its Fig. 7 caption says the inner disk has made
`>250` orbits by `t=600`, while the body text describes about 160; use
absolute simulation time for the principal visual gate and keep orbit counts
as a separately audited quantity. Reference: [Hopkins 2015, §4.2.4 and Figs.
6–8](https://arxiv.org/abs/1409.7395).

The paper's comparison setup is 2D, `gamma=5/3`, no self-gravity, constant
`P=1e-6`, circular velocity in the softened Kepler potential, and the Eq. 48
surface-density profile: floor `0.01`, constant annulus `0.5≤r≤2`, cubic
inner taper, and outer taper width `0.1`. The active q1 initialization
(`EUNHA_KEPLER_PROFILE=hopkins`, `RIN=.5`, `ROUT=2`, floor `.01`, width
`.1`, pressure `1e-6`) matches that profile. Its softening is explicitly
`epsilon=.01175`; the paper does not prescribe one unique small epsilon.
However, q1 has `147456=384^2` generators, whereas the paper calls its particle
setup an effective `256^2` resolution with particles removed beyond the ring.
Treat q1 as a higher-resolution evolution of the same profile, not as an
exact-resolution reproduction; compare initialization/sampling before
claiming direct image-level agreement.

At the measured q1 cadence (`dt≈7.642e-5`, median `2.9 s/step`), a constant
cadence estimate from `t=1.7766` to `t=600` is about 7.83 million steps or
263 days on one H200. This is only a planning estimate, not a runtime
prediction, but it makes timestep stability and GPU/CPU phase optimization
necessary before releasing a continuous long-horizon run. The staged gates
are therefore: resolve the recurrent short-`dt` cause; verify quantitative
errors at the paper's 10-orbit point (`t≈62.83`); validate density and orbital
metrics at `t=120`; then advance to the figure endpoint `t=600`. Treat the
reported ~1000-orbit continuation as an additional endurance target, not a
substitute for the `t=600` comparison.

The later q1 log read recorded another short-step episode at `t=1.78317`
(`dt=1.0564e-5`, then `6.8378e-6`), followed by recovery to
`dt=4.09597e-5` at `t=1.78499`. The latest ledger remained
`dEtot/|Etot0|=-1.434e-10`, `n_ie_le0=0`, and `floor_inj=0`. The recurring
small steps persist without an energy or positivity failure. Job **438070**
remains pending under `QOSGrpGRES` and still targets the earlier confirmed
minimum around `t=1.74353`.

At the 05:24 KST check, the original q1 job remained RUNNING at `t=1.78721`
with nominal `dt=7.64202e-5`, `dEtot/|Etot0|=-1.444e-10`, and no nonpositive
internal energy or floor injection. Diagnostic job **438070** remained
PENDING. `sacctmgr` reports the `part_h200` QoS group cap as
`GrpTRES=gres/gpu=6`; the H200 partition has one 8-GPU node, so the queue
reason reflects the configured QoS allocation cap, not a failure in the
diagnostic binary or its input checkpoint.

### GPU-path source audit for the long-run cost

The RK4 timer's `update` phase includes four calls per accepted step to
`updateDenW2Pressure2DBlend`. In the CUDA branch, each call still rebuilds the
host linked list, builds the GPU cell CSR, serially packs all real and padded
particles from AoS to SoA, performs a weight MPI reduction, uploads particle
and cell arrays, runs GPU tessellation/gradient work, downloads tessellation
results, and writes them back to AoS. This explains why tessellation and data
movement are repeated four times even though force evaluation is on the GPU.
The code already has `GFS_GPU_TESS_PROFILE=1` instrumentation that separates
CSR build, pack, MPI weight reduction, upload, GPU, download, and writeback;
enable it in the next profiling run before changing this path.

A custom two-force-evaluation KDK blend integrator exists as
`exam2d_vph_kdk_int_blend` (`EVOLMETHOD=3`) for other 2D drivers, but
`Exam/Kepler/kepler.c` dispatches only method 1 to RK4 and all other values to
`exam2d_vph_int`. Therefore changing Kepler's parameter to `3` would not select
the custom KDK routine. A future KDK/individual-step phase needs an explicit
Kepler dispatcher and boundary/post-stage contract; it must be tested as a
separate continuation from a verified checkpoint.

At the 05:27 KST log read, q1 had advanced to `t=1.79083` with
`dt=5.86909e-5`, `dEtot/|Etot0|=-1.433e-10`, and zero nonpositive-energy or
floor-injection events. The short-step episodes remain transient and the
baseline continues to advance. Job **438070** is still pending on the
`part_h200` group-GRES cap.

### Exact-zero-weight GPU fast path benchmark (2026-10-28)

Three matched A10 runs restarted from the same q1 checkpoint at `t=0.825` and
advanced 263 steps to `t=0.845`: baseline job **438011**, W=0 fast-path job
**438015**, and an independent fast-path rebuild/repeat **438016**. All three
completed successfully. After sorting by particle ID, the final position,
velocity, mass, volume, internal energy, density, pressure, and weight arrays
were bitwise identical between baseline and both fast-path outputs. Thus the
fast path preserved this tested W=0 trajectory exactly.

The median logged step time after the first ten steps was 2.4 s for baseline
and 2.0 s for both fast-path runs (about 17% faster); full job elapsed times
were 13:43, 11:57, and 11:52. The sampled GPU profile attributes most of the
per-call device saving to `density_copy_A` (about 140 ms on the baseline
sample versus 109–112 ms on the repeated fast-path samples). CPU-side packing
and CSR construction remain visible costs, and the RK4 update still rebuilds
the geometry four times per accepted step. The improvement helps, but by
itself is far too small to make a direct `t=600` run practical at the current
step cadence.

The tested fast path is enabled only after an exact global check that every
packed weight, including ghost/padding entries, is zero; any nonzero weight
falls back to the Laguerre expression. In this project, `W=0` therefore means
all cell weights `w_i=0`, which gives the Voronoi tessellation. The longer
q1 continuation to `t=2.221` is still running the baseline binary (job
**437212**, latest inspected `t=1.81356`,
`dEtot/|Etot0|≈-1.45e-10`, `n_ie_le0=0`, and `floor_inj=0`). Its dependent `t=6.5` continuation, job **437213**,
is staged to use the separately hashed W=0 fast-path executable after the
parent completes. The q1 conditional checkpoint-resume script is present in
`halo28_from2221_to65_resume1`; its input hashes pass, and q1 fallback job
**438008** is confirmed PENDING with dependency `afterany:437213`. The earlier
assumed ID **437772** and job **437779**
are actually q2 420² jobs. Keep q1 and q2 job chains distinct when checking or
resuming them. No replacement or duplicate evolution was started.

### Corrected global orbital-cap diagnosis (2026-10-28)

The q1 `t=0.825` checkpoint has a generator essentially at the softened
potential center (`r=1.08e-11`, `rho=0.01`, `mass=2.1143e-6`). Since
`kepler_dt_limit` minimizes over all 147,456 cells, this one floor-density
cell sets the `eta=0.06` cap to `7.6420e-5`; the value matches the nominal
accepted q1 step. By contrast, the cap at the active annulus inner edge
`r=0.5` is `2.1222e-2`. Thus the cap does **not** explain dips below
`7.642e-5`, but it does explain the unusually small nominal step and most of
the long-run cost. The prior claim that the central cap could not explain the
q1 timestep was incorrect because it used only `r=0.5` and missed the
near-center floor cells.

The matched A10 pair from the same `t=0.825` checkpoint and fast-path binary
isolates `EUNHA_KEPLER_DTETA`: with `eta=0.06`, the endpoint step at `t=0.845`
was `4.50e-5`; with `eta=0`, it was `2.052e-4`. Active-annulus orbit and
angular-momentum metrics were effectively identical at this short endpoint
(`active_radius_error_rms≈6.38e-10` in both; relative active `Lz` change
`≈8.42e-13`). Total-energy drift grew from `-1.58e-12` to `-9.40e-10`, so
the no-cap branch is faster but has measurably larger integration error.

No-cap continuation job **438040** is still running from the eta=0 branch at
`t=0.845`; the latest log is at `t=1.25082`, with total-energy drift
`-1.427e-8`, `n_ie_le0=0`, and `floor_inj=0`. Over its latest 100 steps, `dt`
has median `2.063e-4`, minimum `7.840e-6`, maximum `4.914e-4`, with one value
below `1e-5`. Its `t=1.5` orbital/conservation audit has not run yet.
Keep the already-running capped q1 chain intact. Use the `t=1.5` audit to
decide whether a no-cap continuation from a verified shared checkpoint can
replace the cap for later long-horizon branches; do not infer long-time
accuracy from the `t=0.845` pair alone.

### Same-executable dteta control (2026-10-28 05:50 KST)

The latest `438040` log marker is `t=1.25728` (`dt=3.3863e-4`); its latest
energy ledger reports `dEtot/|Etot0|=-1.443e-8`, with no nonpositive internal
energy and no floor injection. This is still an intermediate result.

An executable-provenance check found that `438040` uses SHA-256
`24ad8026...9dd86ce3`, while the running H200 reference `437212` uses
`75844269...67c0284`. Their `params.dat` hashes match, but the executable
hashes do not. Therefore the existing eta=0.06 H200 endpoint is useful as a
reference, but comparing it directly with `438040` cannot isolate the
timestep-cap effect.

To make that comparison controlled, job **438083** was submitted from the
same t=0.845 checkpoint `khout.010920.dat`, with the same `24ad8026...`
executable and identical `params.dat` as `438040`; its only physics-setting
change is `EUNHA_KEPLER_DTETA=0.06`. It requests exact stop at `t=1.5`, with
the same output cadence and A10 partition. The run directory contains a
SHA-256 manifest for the executable, parameters, checkpoint, and job script.
At the submission check, both `438040` and `438083` were RUNNING on `syn01`;
the H200 q1 baseline `437212` was also still RUNNING. Compare the two A10
endpoints by particle ID and conservation metrics before deciding whether
removing the orbital cap is acceptable. Keep the earlier H200 comparison
marked as cross-binary evidence.

### Matched no-cap t=1.1 audit (2026-10-28)

The already-completed pair **438003/438004** provides a clean short-horizon
comparison with the same executable SHA-256
`75844269...67c0284`, identical `params.dat`, and identical restart
`khout.013200.dat` at `t=1.0000125425`. Both reached exactly `t=1.1` with
fixed IDs and masses. The eta=0.06 run took 4,139.427 s; eta=0 took 1,126.546
s, a 3.67x wall-time reduction (1,334 versus 371 reported RK4 steps).

The established 15-check Kepler gate passes for both. Their active-disk
mass-weighted absolute density changes are effectively identical
(`4.7819574e-4`); active-disk position errors are `2.0847e-8`; GIZMO-annulus
velocity L1 is `1.256e-8`; and active-disk angular-momentum drift is about
`-8.5e-13`. A direct ID-matched endpoint comparison gives RMS position
difference `2.83e-12` relative to initial radius, RMS velocity difference
`1.25e-10` relative to circular speed, and mass-weighted relative density
difference `1.35e-10`. Per-ID mass is bitwise identical.

The cost of removing the cap is a larger energy error: `-3.7473e-9` at eta=0
versus `-6.6298e-11` at eta=0.06, still far below the current `1e-3` energy
gate. This supports no-cap as a promising performance candidate through
`t=1.1`; it is not evidence of long-time stability. The same-executable
continuations **438040/438083** to `t=1.5` remain the next controlled check.

### OMP RHS candidate and t>6 no-cap continuation (2026-10-28)

The previous cached-RHS candidate **438007** reached the same `t=1.1` endpoint
as capped baseline **438003**, and the saved physical fields (`x`, `y`, both
velocity components, mass, volume, internal energy, density, pressure, and
`w2`) are bitwise identical. Its full run time was 4,190.382 s versus
4,139.427 s for the baseline, but it ran on a different A10 node, so this is
inconclusive rather than evidence of a speedup or regression. The instrumented
CPU phase-sum median was 2.4 s for both.

A new isolated candidate was built from the exact 24ad W=0 source overlay. It
adds `schedule(static)` OpenMP parallel loops only to the four independent
per-cell K1–K4 RK4 RHS loops and initializes the lazy SEDOV mode before K1.
Its executable SHA-256 is `affdc7bf...f3f20dc7`; the user's working `src/`
tree remains untouched. Job **438089** runs the 24ad baseline and this
candidate sequentially on one A10 allocation, both at OMP=2 and eta=0, from
the same `t=0.845` checkpoint to `t=0.855`. At the latest check the baseline
case had reached the endpoint and the OpenMP candidate was running; no
speedup is claimed until both logs and native fields are compared.

The no-cap `t=1.5→6.5` continuation, job **438092**, is staged with dependency
`afterok:438040`. At job start it runs the native 15-check gate and radial
audit on the exact no-cap t=1.5 checkpoint; only if the gate passes does it
continue eta=0 to exact `t=6.5`, writing checkpoints every 0.5. It has a
two-day A10 limit and runs independently of same-executable cap control
**438083**. The purpose is to reach the requested `t>6` gate promptly while
preserving the matched cap comparison as a separate validation.

At the latest read, no-cap job **438040** was at `t=1.32284`, with relative
total-energy drift `-1.609e-8`, zero nonpositive internal-energy cells, and
zero floor injection. Cap control **438083** and H200 baseline **437212**
were still running; the original capped continuation **437213** remained
dependency-pending.

### Static-RHS thread-count result and current continuation state (2026-10-28 06:17 KST)

The follow-up A10 job **438093** ran the same static-RHS candidate executable
from the same `t=0.845` checkpoint to exact `t=0.855`, first with
`OMP_NUM_THREADS=2`, then with 8, on one allocation. Both endpoints have the
same sorted particle IDs and are bitwise identical in `x`, `y`, `vx`, `vy`,
mass, volume, internal energy, density, pressure, and `w2`. The recorded case
times are 130.377 s and 127.034 s respectively, only a 2.6% difference; this
does not establish a useful OpenMP speedup. The job log records exact endpoint
`t=0.855` and `EXIT:0` for both cases.

At the 06:17 KST progress read, the controlled no-cap run **438040** had
reached `t=1.35350`, with relative energy drift `-1.675e-8`, no nonpositive
internal-energy cells, and no floor injection. Its matched `eta=0.06` control
**438083** was at `t=0.877361`, with drift `-4.584e-12` and the global cap
`dt=7.64202e-5`; this slower branch is still progressing. The original capped
H200 baseline **437212** was at `t=1.85252`, with drift `-1.534e-10` and no
nonpositive internal-energy cells or floor injection. The no-cap `t=1.5→6.5`
continuation **438092** remains staged behind `afterok:438040`; the original
capped continuation **437213** remains staged behind `afterok:437212`.

The completed W=0 GPU profile in `tess_profile_w0global_fast_t0825_a10_20261028`
sampled eight tessellation calls. Median per call was 21.317 ms for cell-CSR
construction, 18.961 ms for AoS-to-SoA packing, and 22.970 ms for writeback.
On device, Kernel A (logged as `density_copy_A`) took 104.159 ms; pressure
kernel B took 0.034 ms, scan 0.011 ms, host round-trip 0.216 ms, and face
scatter 1.661 ms. The total tessellation path median was 173.709 ms. Together
with the completed 2-vs-8 thread result, this points to GPU Kernel A and the
repeated geometry work as the next optimization target; adding threads to the
RK4 RHS alone is not a demonstrated throughput improvement.

The pending no-cap continuation **438092** was re-audited at 06:20 KST: its
input/source SHA-256 manifest passes, the submitted script retains
`afterok:438040`, and it enforces the native t=1.5 gate before selecting the
restart and evolving to exact t=6.5. No staging change was needed.

### W=0 seed-bin pruning candidate result (2026-10-28 06:35 KST)

An isolated CUDA candidate added geometric lower-bound pruning to the W=0
nearest-support prepass and 3×3 seed-bin traversal. Job **438105** ran the
control (`affdc7bf…f3f20dc7`) and candidate (`3fbac30a…b69b12ea`)
sequentially on one A10, with the same t=0.845 restart, OMP=8, eta=0, profile
instrumentation, and exact t=0.855 stop. Slurm reports COMPLETED, exit 0.
Sorted-ID endpoints are bitwise identical for position, velocity, mass,
volume, internal energy, density, pressure, and weight.

The candidate is not a speedup: whole-case time was 127.810 s versus 126.442
s for control. Across the eight profiled tessellation calls, median Kernel A
time rose from 110.111 ms to 124.677 ms, and total tessellation time rose
from 179.958 ms to 189.992 ms. Reject this pruning candidate for production;
keep its code isolated in the archived test directory. This confirms that
the next optimization must target Kernel A without adding enough per-bin
bookkeeping to outweigh saved candidate work.

At the 06:35 KST follow-up, no-cap **438040** was at `t=1.40954`
(`dEtot/|Etot0|=-1.832e-8`), matched capped control **438083** at
`t=0.898091` (`-3.217e-12`), and H200 baseline **437212** at `t=1.87424`
(`-1.582e-10`). All three still report zero nonpositive internal-energy cells
and zero floor injection. The two t=6.5 continuations remain dependency-pending.

### 420² refined W=0 geometry-cost growth (2026-10-28 11:36 KST)

The running refined A10 job **437771** is a separate `420²`, fixed-mass,
zero-weight Voronoi branch. Its restart manifest pins binary SHA-256
`7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284`.
The source archive used for that binary already streams candidate halfplanes
from CellCSR and avoids the former large per-thread candidate arrays and
insertion sort; repeating that optimization would be redundant. The workspace
`src/Exam/exam_gpu.cu` has a different hash and was not edited.

The call-indexed profile shows a large, evolving cost inside Kernel A
(`density_copy_A`), while host-side CSR, packing, and writeback costs stay
roughly steady:

| tessellation call | Kernel A | GPU section | whole tessellation |
|---:|---:|---:|---:|
| 2,048 | 34.9 ms | 41.5 ms | 82.9 ms |
| 66,560 | 162.4 ms | 168.8 ms | 213.7 ms |
| 67,584 | 324.3 ms | 330.9 ms | 373.0 ms |
| 77,824 | 407.9 ms | 414.6 ms | 457.6 ms |
| 78,848 | 434.0 ms | 440.6 ms | 483.3 ms |

At the last sample, the RK4 step took about 3.6 s and used `dt=7.64202e-5`.
The same call reported only 276 stencil-expanded cells out of 176,400 and no
halo-clipped cells. Thus broad stencil expansion does not explain the growth;
the number of candidate planes actually clipped and polygon-edge visits is not
currently recorded. The next useful profiling sample is one representative
late-state call with per-cell candidate, accepted-plane, edge-work, and peak
corner counts. Keep the production geometry unchanged until those counts
identify a safe optimization.

At the 11:36 KST scheduler check, **437771** was still RUNNING; its log had
advanced to `t=1.47674`, step 19,604, with no nonpositive internal-energy
cells or floor injection in the nearby ledger. Its exact `t=2.221` gate is
incomplete. The `afterany` continuation **437772** remains the existing
recovery path if the one-day allocation expires before the endpoint.

At the 11:49 KST check, **437771** was still RUNNING at `t=1.48620`, step
19,729. Profile call 78,848 confirms the geometry-kernel growth persists;
CSR (17.3 ms), packing (15.9 ms), and writeback (8.2 ms) remain small beside
Kernel A (434.0 ms). The latest complete snapshot is still before the
scheduled `t=1.5` checkpoint, so the late-state Nsight Compute sample must
wait for that checkpoint before it can be staged without interfering with the
active run.

At 11:57 KST, the same run was still advancing at `t=1.49416`, step 19,838.
The nearby RK4 ledger at step 19,794 reports total-energy relative drift
`-1.693e-10`, `dEhyd/Ehyd0=2.128e-4`, no nonpositive internal-energy cells,
and zero floor injection. Short steps at `t=1.49401–1.49402` (`1.72e-5` and
`5.00e-6`) recovered to `7.63e-5` within two more steps. The original run has
`GFS_VARPRESS_GEOM_AUDIT` disabled, so those dips do not identify their
limiting cell. The isolated one-step Nsight Compute script also enables this
read-only audit and is syntax-checked; submit it after the exact `t=1.5`
checkpoint is complete. That checkpoint had not yet appeared at this read.

### q1 uncapped continuation live check (2026-10-28 11:59 KST)

The matched q1 no-cap continuation **438092** remains RUNNING on A10 `syn01`
with a 48-hour allocation. At elapsed 4:47:59 its latest log reached
`t=2.06784`, step 16,959. The recent accepted steps ranged from
`1.04e-4` to `2.37e-4`, versus the `7.642e-5` nominal eta=0.06 cap in the
matched capped branch. Per-step wall time was about 4.0 s. The nearby ledger
reports `dEtot/|Etot0|=-6.018e-9`, `dEhyd/Ehyd0=2.238e-4`, `n_ie_le0=0`, and
zero floor injection. This is promising throughput for reaching the
intermediate `t=6.5` gate, not yet a long-time accuracy result; the endpoint
audit remains required before the no-cap branch is extended further.

### 420² late-state Nsight diagnostic submitted (2026-10-28 12:05 KST)

The refined W=0 A10 trajectory **437771** produced `khout.019919.dat` at
`t=1.500020939799`. Its native header is complete and consistent: 176,400
cells, 576-byte cell stride, stress fields present, and exact file size
101,606,468 bytes. The isolated single-step Nsight Compute plus geometric-CFL
audit was submitted as **438386** from that checkpoint. It uses a separate
diagnostic directory, profiles the first `voronoi_tessellate_kernel` launch,
and does not modify the advancing trajectory. The one-step run completed with
exit 0 and left the production trajectory untouched, but Nsight reported
`No kernels were profiled`; the geometry audit did run. It measured a maximum
swept-volume rate of 212.084 per code-time unit at local array index 176258.
The implemented rate is `sum(|Vdot_face|)/cell_volume`, with a timestep bound
`C/rate`; for the run's `C=0.1`, this is `4.72e-4`. Nearby accepted steps are
about `7.64e-5`, so this bound was not the active timestep limit at t≈1.5.
The audit used a forced `1e-5` one-step stop, giving `h*rate=0.00212`; that
number is diagnostic scaling, not the run's adaptive step.

Changing the Nsight filter to the mangled CUDA symbol did not resolve capture:
**438387** also completed in 13 seconds with `No kernels were profiled` and no
`.ncu-rep`. The CUDA kernel does exist in the executable, so the remaining
suspect is the MPI launch boundary: Nsight wrapped `mpirun`, while Intel
Hydra/Slurm starts the rank separately. A third one-step attempt, **438388**,
places Nsight inside the one-rank MPI launch and uses a fresh output directory.
The capture result was `ERR_NVGPUCTRPERM`: Nsight connected to `eunha2.exe`
but the account cannot access GPU performance counters on the A10. The one-step
simulation completed its geometry audit, while the Slurm job exited 1 and
produced no `.ncu-rep`. Stop further Nsight attempts unless counter access is
enabled; use CUDA event timing or add explicitly gated application counters.

At the same scheduler check, the no-cap q1 continuation **438092** had reached
`t=2.08162`, step 17,082, with `dEtot/|Etot0|=-6.051e-9`, no nonpositive
internal-energy cells, and zero floor injection. It remains short of the
required t=6.5 gate. The 420² W=0 trajectory **437771** also remains active;
its checkpoint-near-t=1.5 diagnostic did not alter it. MPI one-rank fix job
**438363** has started on H200; its four-rank dependent job remains queued.

At the follow-up log check, no-cap q1 **438092** had advanced to `t=2.08542`,
step 17,109, with accepted `dt=1.80614e-5`; it has not reached t=6.5. H200
one-rank fix job **438363** remained RUNNING and its four-rank dependent job
**438364** remained pending. The refined 420² W=0 job **437771** was still
RUNNING near `t=1.5024`. No existing job was changed.

### q1 t=2.0 geometry-CFL audit and CUDA event profile (2026-10-28)

The complete q1 checkpoint `khout.016489.dat` is at `t=2.000012530836` and
stores a preceding step of `2.29852e-5`. Isolated one-step job **438389**
completed with exit 0. At its four RK4 RHS stages, the maximum swept-volume
rate ranged from 920 to 997 per code-time unit (local array index 138943).
With `GFS_VARPRESS_GEOM_CFL=0.1`, the corresponding geometric bound is about
`1.00e-4–1.09e-4`, over four times the stored preceding step. Thus the
geometric swept-volume bound was not setting that small step at t≈2; identify
the other active timestep constraint before changing this limiter. The
diagnostic's own `1e-5` stop was forced and is not an adaptive-step result.

Application CUDA events work without Nsight permissions. In this q1 build,
`density_copy_A` took 396–437 ms per 384² tessellation call; B and scan each
took under 0.05 ms, host roundtrip under 0.5 ms, and scatter about 1.7 ms.
The whole tessellation call took 464–502 ms. The 388 stencil-expanded cells
and zero halo-clipped cells are not enough by themselves to explain Kernel A's
cost. This confirms Kernel A geometry clipping remains the optimization
target, while GPU hardware-counter analysis is blocked by account permissions.

At the latest log read, no-cap q1 **438092** had advanced to `t=2.09195`,
step 17,153, with `dt=1.57775e-4`; the t=6.5 endpoint and its strict audit
remain incomplete.

### Latest live branches (2026-10-28 12:18 KST)

No-cap q1 **438092** advanced to `t=2.09359`, step 17,169. Its nearby ledger
has `dEtot/|Etot0|=-6.068e-9`, no nonpositive internal-energy cells, and zero
floor injection. The refined 420² W=0 branch **437771** reached `t=1.51418`;
its nearby ledger has `dEtot/|Etot0|=-1.663e-10`, no nonpositive internal
energy, and zero floor injection. Both remain RUNNING. H200 one-rank fix
**438363** is RUNNING and its four-rank dependent **438364** remains pending.
The q1 t=2 geometry audit **438389** is COMPLETED, exit 0; it wrote only to
its isolated diagnostic directory.

### Exact-output CUDA clipping A/B submitted (2026-10-28 12:40 KST)

The isolated candidate changes only the exact Sutherland-Hodgman path in
`gpu_clip_voronoi_halfplane`: it swaps polygon/scratch array views after each
accepted halfplane instead of copying the completed polygon back. The clipping
order and arithmetic are unchanged. Candidate executable SHA-256 is
`6718551061d70e24e09f0fe000a35b5cf1ac8c1a15b579f1e97e6d907b42baa0`; the
q1 baseline executable is
`a6f16e8b1e143f70e8b8c32155008077bd84e813f7f4e4aff2f04b23dfd04572`.

The candidate has slightly higher static register/stack use than baseline on
the available GPU cubins (A100 `sm_80`: 128/2976 B versus 127/2960 B). The
runtime comparison is therefore necessary; static resources alone do not show
whether the removed copy traffic wins. H200's test allocation forecast was
weeks away, while the A100 PCIe node had an immediate slot, so the paired test
uses one `gpu:A100:1` allocation. Job **438400** ran three paired one-step
restarts from the identical 384² native checkpoint at `t=2.000012530836` and
recorded CUDA event timings. This diagnostic did not modify the ongoing q1
trajectory. The initial raw-byte equality gate was invalid and is described
below; the numerical comparison uses the stored fields and ignores native
record pointers and padding.

The six one-step executions in **438400** all returned exit 0; the batch
wrapper then exited 4 because its raw `cmp` gate compared pointer-bearing
native records byte for byte. `treevorostressrk4particletype` stores linked-list
pointers at offsets 0 and 568, so raw snapshot equality is not a valid physics
comparison. A structured comparison that ignores those pointers and ABI
padding found all primary state fields (positions, velocities, masses,
density, pressure, volume, internal energy, and W values) exactly equal
between baseline and candidate. Repeated runs of each binary also agree in all
numeric fields. Only derivative/stress scratch fields differ, at no more than
`2.85e-13` absolute (`lap_vy`) and below `4.8e-17` relative to each field's
maximum magnitude.

Across 12 RK4 tessellation calls per build, `density_copy_A` averaged
`263.97 ms` baseline and `241.57 ms` candidate (9.27% shorter); whole
tessellation averaged `326.73 ms` and `302.69 ms` (7.94% shorter). However,
the first tessellation call consistently reports two more faces in the
candidate (`884778` versus `884776`); the other three stage counts match.
The one-step primary state and cell volumes are unchanged, but this topology
count difference still needs a short evolution check before the optimization
is considered safe. The original raw-byte comparator is being replaced with a
structured numerical comparison.

`compare_kepler_native_state_20261028.py` now models and validates the complete
576-byte record layout, checks IDs and metadata, and excludes the serialized
`next`/`bp` pointers and ABI padding. It verified all three one-step pairs and
the same-binary repeats. The matched `t=2.000012530836` → `t=2.001` evolution
check was submitted as **438414** on `a100_pcie` (`gpu:A100:1`); the scheduler's
test allocation forecast was 13:35 KST. Both variants use the same restart,
one rank, and the q1 W=0 fixed-mass settings. The long q1 trajectory and all
other existing jobs remain untouched.

The short evolution **438414** completed successfully on A100. Baseline and
candidate both advanced five accepted RK4 steps from `t=2.000012530836` to
`t=2.001`, with the same step sizes, endpoint metadata, and exact primary
state fields. Both had `n_ie_le0=0` and zero floor injection; the final energy
ledger was `dEtot/|Etot0|=-1.267e-11`. Across the 20 profiled tessellation
calls, only the first face count differed (candidate +2); the remaining face
counts matched. Mean full-tessellation time was `323.01 ms` baseline and
`306.24 ms` candidate (5.47% shorter). Non-primary derivative/stress fields
differed only at roundoff scale (maximum absolute `5.26e-13`, maximum relative
to field maximum `4.73e-16`). This supports using the candidate for an isolated
longer Kepler branch, while the first-call degeneracy remains noted.

The primary-paper endpoint was rechecked against Hopkins et al.'s GIZMO paper:
Figure 7 is the 2D Keplerian disk density map at `t=600`, and the inner
`r≈0.1` region has completed over 250 orbits without decaying or disrupting.
The final acceptance therefore still requires a density-field comparison at
that horizon; the present t=6.5 candidate branch is only an intermediate gate.
[Source: original GIZMO paper PDF](https://www.tapir.caltech.edu/~phopkins/Site/GIZMO_files/gizmo.pdf).

The isolated candidate continuation from `khout.016489.dat` (`t=2.00001253`)
to the same t=6.5 gate was submitted as **438416** on A100 PCIe, one GPU, with
a 48-hour limit. It uses the same q1 W=0, fixed-mass, one-rank controls and
will run the existing long-time and radial audits at completion. The
test-only scheduler estimate was 13:40 KST, but actual allocation began at
12:55:55 KST on `syn103`. After 33 seconds it had reached `t=2.00159`, step
16,499; the nearby energy ledger is `-1.747e-11`, with no nonpositive
internal-energy cells or floor injection. The main q1 baseline **438092**
continues independently and had reached `t=2.1444` at the latest log read.

### Direct shared-log check (2026-10-28 14:52 KST)

The current grammar Slurm queue does not contain the syntax jobs 438092/438416,
so it cannot provide their scheduler states. The shared run logs nevertheless
advanced at this check and contain no terminal `EXIT` marker or t=6.5 audits:

- **438416**, isolated A100 clip-swap continuation: `t=2.21970`, step 18,286;
  log updated 14:52:25 KST. The last logged candidate tessellations are about
  325 ms. From the 14:45 log sample, its short-window progress was about
  `0.10 t/hour`; this is only a rough projection, not a stable throughput
  measurement.
- **438092**, A10 no-cap baseline continuation: `t=2.28828`, step 18,929;
  log updated 14:52:25 KST. Its last accepted timesteps vary from roughly
  `3e-6` to `1.4e-4`, so a short-window extrapolation is especially noisy.

Neither run has reached the t=6.5 physics gate. Keep both outputs intact and
resume from their latest complete checkpoint if a wall-time limit interrupts
them; do not replay either trajectory from its original checkpoint. The
current grammar Slurm view is not authoritative for these syntax job IDs.

**Follow-up at 14:55 KST:** both logs advanced again. The A100 candidate is at
`t=2.22610`, step 18,337; the A10 baseline is at `t=2.29242`, step 18,968.
Using the observed job start times and constant average `Δt_sim/Δt_wall` only
as a planning estimate, the A100 job has about 45.0 wall hours left and would
need about 37.8 hours to t=6.5. The A10 job has about 40.3 hours left and
would need about 41.1 hours, so it may stop just short of the gate; the
half-unit snapshots allow a continuation from the last complete checkpoint.
Both projections are uncertain because accepted `dt` varies substantially.

### Sparse buffered-area continuation from the last complete checkpoint

The 54,312-site W0 buffered-area branch passed its exact `t=0.825` gate and
was continued by two older binaries. Both stopped at `t=2.50986` with
`[GPU_CLIP_CAPACITY] nc=48 capacity=48`, followed by
`[GPU_FINAL_VOLUME] invalid=1` and `EXIT:1`. The last complete checkpoint is
`khout.001318.dat` at `t=2.4753977144699952`; it was audited against the
exact `t=0` state in
[the checkpoint receipt](KEPLER_Q2_SPARSE_PROFILE_T2P475_CHECKPOINT_AUDIT_2026-10-28.json).

At that checkpoint, fixed masses are bitwise equal by ID, all weights are
exactly zero, every cell is finite with positive volume and internal energy,
the periodic volume sums to 36, and whole-box relative energy drift is
`1.84e-8`. For the active `0.5<=r<2` disk, mass-weighted density change is
`1.44e-4`, absolute entropy-log change `4.92e-6`, mass-weighted orbit-position
error `3.62e-8`, and disk-cohort angular-momentum change `-1.41e-13`. The
whole-box Lz diagnostic changes by `-0.94%`; it is not a conservation gate on
this Cartesian periodic domain. Cells initially beyond `r=2.35` already show
large edge motion, so the checkpoint does not prove stability beyond the
active disk or beyond this time.

Job **438533** was submitted to H200:1 to resume from this exact checkpoint
with clip-swap binary SHA-256
`6718551061d70e24e09f0fe000a35b5cf1ac8c1a15b579f1e97e6d907b42baa0`, retain
W0/fixed mass and the sparse-area state, cross the previous `t=2.50986` failure,
and reach exact `t=10` with the existing strict audit. It writes into a new
run directory and does not alter the failed runs. At 15:05:58 KST, syntax
Slurm reported it `PENDING (QOSGrpGRES)`; the H200 partition has one eight-GPU
node in mixed use. Jobs 438092 and 438416 were still `RUNNING` at that same
syntax query. The local grammar queue is a separate cluster; use the syntax
login alias for scheduler control and the shared logs for scientific output.

**15:06:50 KST syntax query:** 438533 is still `PENDING (QOSGrpGRES)` and has
not started. The two original jobs remain `RUNNING`; their shared logs advanced
to `t=2.24549` (438416) and `t=2.30452` (438092) by 15:06:53 KST. Neither has
produced the t=6.5 audit yet. The sparse continuation is staged and queued,
not counted as an active simulation until Slurm allocates its H200.

**15:12 KST syntax query:** 438092 (`KepEta0To65`) remains `RUNNING` on
`syn01` at elapsed 8:01; its log had reached `t=2.31122`, step 19,166.
438416 (`KepClipSwapTo65`) remains `RUNNING` on `syn103` at elapsed 2:16;
its log had reached `t=2.25602`, step 18,599. Neither has reached the t=6.5
audit. The H200 sparse continuation 438533 remains `PENDING (QOSGrpGRES)` and
has not started. Keep all three run directories intact; do not submit duplicate
continuations while these states persist.

**15:16 KST syntax query:** all three scheduler states are unchanged. Job 438092
has reached `t=2.31448`, step 19,206; its latest RK4 ledger reports relative
total-energy drift `-6.392e-9`, no nonpositive internal energy, and no floor
injection. The sampled step took 3.9 s and accepted `dt=4.46e-5`. Job 438416
has reached `t=2.26147`, step 18,652; its latest ledger reports drift
`-6.807e-10`, with no nonpositive internal energy or floor injection. Its
sampled step took 2.7 s with `dt=1.39e-4`; the candidate-stream marker accounts
for all candidates. Neither run has reached the t=6.5 gate. The last complete
scheduled checkpoint in each run directory is still `t≈2.0`; both are expected
to write their next `t=2.5` checkpoint after reaching that time. The A10
baseline's elapsed-average projection is borderline for the remaining wall
limit; the A100 candidate has more wall-time margin, but variable `dt` makes
both projections provisional. Job 438533 remains queued for `QOSGrpGRES`.

**15:18 KST syntax query:** 438092, 438416, and 438533 retain their prior
scheduler states. The latest log times are `t=2.31757`, step 19,238, for 438092
and `t=2.26559`, step 18,694, for 438416. The A100 log briefly reached
`dt=1.40e-6` at step 18,690, then recovered to `1.72e-4` two steps later; this
sample is a transient timestep dip, not a sustained collapse. Both trajectories
remain below their t=2.5 checkpoint and have no t=6.5 audit yet. Job 438533 is
still `PENDING (QOSGrpGRES)`.

### Exact q1 t=2.000 checkpoint physics audit (2026-10-28)

The native q1 W=0, fixed-mass checkpoint `khout.016489.dat` at
`t=2.000012530835539` passes all 15 checks in the existing long-gate audit.
For the initialized disk `0.5<=r0<2`, mass-weighted absolute density change is
`5.0531e-4`, absolute entropy-log change `1.7088e-4`, orbit-position error
`7.1412e-8`, and relative angular-momentum drift `9.34e-13`. Whole-domain
relative energy change is `-2.6942e-8`; volume sums are `36.00000000` at both
ends, `rho=m/V` agrees to `5.93e-8` relative, and mass and W=0 are exact. In
the GIZMO-style `0.8<r<1.2` annulus, L1 velocity error is `2.95e-8` and
plateau-corrected density error is `5.93e-4`.

The companion [radial audit](KEPLER_Q1_NOCAP_T2_CHECKPOINT_RADIAL_AUDIT_2026-10-28.json)
shows the active disk's inner `0.5<=r0<0.75` density error is `1.046e-3` and
entropy error `4.13e-5`. It also shows large excursions outside the initialized
disk: for `2.75<=r0<5`, density error is `0.312`, radial speed `0.166`, and
periodic-aware ballistic position error `0.03936`. The long-gate receipt also
lists `1.995` for this exterior bin, but that diagnostic compares unwrapped
coordinates; it is not the periodic-domain error. Density and radial-speed
metrics agree between the two audits. Keep the exterior density/velocity
excursion separate from the passing active-disk gate and investigate its
effect before making a long-time stability claim. The `r<0.1` cells are a
small core cohort (134 cells) inside the initial `r=0.5` inner edge, so this
run is not yet the paper's `r≈0.1`, `t=600` benchmark.

The [15-check receipt](KEPLER_Q1_NOCAP_T2_LONG_GATE_AUDIT_2026-10-28.json)
uses the exact t=2 snapshot for all state and orbit metrics. Its cumulative
injection fields are read from the live continuation log; a separate scan of
all 5,687 `[RK4E]` records through `t=2.3211` found `floor_cum`, `de_cum`, and
`sfl_cum` identically zero. Thus the snapshot gate is valid and no injection
has occurred through that later log time.

At the 15:22 KST syntax query, **438092** and **438416** remained `RUNNING` at
`t=2.32110` and `t=2.27130`, respectively; neither had written its t=2.5
checkpoint. **438533** remained `PENDING (QOSGrpGRES)`. The t=6.5 audits remain
the next long-horizon gates.

**15:23 KST syntax/H200 check:** 438092 and 438416 remain `RUNNING`; their
logs advanced to `t=2.32212` (step 19,293) and `t=2.27238` (step 18,769).
438533 remains `PENDING (QOSGrpGRES)`. `sinfo` reports the H200 partition up
with one eight-GPU node (`syn104`) in mixed use. The pending reason is the
Slurm GRES group limit, rather than an offline H200 partition.

**15:28 KST follow-up:** 438092 advanced to `t=2.32772`, step 19,341; 438416
advanced to `t=2.27899`, step 18,832. Both remain `RUNNING` and below the t=2.5
checkpoint. The H200 limit is now identified: partition `h200` applies QoS
`part_h200`, whose GRES group cap is six GPUs. Node `syn104` has eight H200s
with seven allocated, so one physical card is free while the partition QoS
group is saturated. The alternate `sjjlee_hgpu` partition allows only QOS
`sjjlee_unlim`; an `sbatch --test-only` request using that QOS was rejected as
`Invalid qos specification`, and the user's association lists `normal` only.
No alternate H200 submission was made. The existing 438533 request remains
queued for the H200 partition's QoS capacity.

### Area-weighted CFL continuation result (2026-10-28)

The exact `t=2.0001076440463654` checkpoint from area-CFL job 438384 was
audited against the branch's own t=0 snapshot with the existing 15-check
long-gate audit; all checks passed. Its active-disk density, entropy, orbit,
and angular-momentum errors agree closely with the no-cap baseline at t=2.
The machine-readable receipt is
[the area-CFL t=2 audit](KEPLER_Q1_AREA_CFL_T2_LONG_GATE_AUDIT_2026-10-28.json).
The outer-bin orbit-position field in that generic receipt uses unwrapped
coordinates; the [periodic-aware radial receipt](KEPLER_Q1_AREA_CFL_T2_PERIODIC_RADIAL_AUDIT_2026-10-28.json)
gives the correct wrapped error `0.03936` for `2.75<=r0<5` and confirms the
active-disk position error remains `7.14e-8`.

That does not make the long branch viable: it aborted at `t=2.1450139886`,
step 15412, during the RK4 k3 predictor when outer-transition cell 92210 had
`U=-5.94420e-10` (`EXIT:92`). Its accepted start had `U=9.77513e-10`,
`rho=0.06404`, and `P=5.61969e-6`; its initial radius was 2.328, outside the
initialized `0.5<=r<2` disk. A native snapshot read shows the same cell had
rho 0.02143 at t=2.0001 and 0.02020 at t=1.515. The failed-stage snapshot is
not restartable; `khout.014793.dat` at t=2.0001 is the last valid state.

Source inspection confirms `GFS_AREA_CFL_APPLY` replaces the minimum
face-acoustic timestep by an area-weighted acoustic crossing time, while
retaining non-face bounds. This likely explains the speed increase but is not
positivity preserving by itself. A near-coincident generator was found about
`9.2e-5` away from cell 92210 at failure start; the available checkpoint does
not contain enough face topology to prove that pair caused the negative
predictor. Treat this branch as a failed t>2.145 candidate, not as a long-time
stability improvement.

The last syntax Slurm query showed 438092 (A10 no-cap) and 438416 (A100
clip-swap) RUNNING; 438533 (H200 sparse clip-swap), 438394 (MPI2), and their
dependent jobs remain pending. The failed area-CFL job is absent from the
queue. No jobs were cancelled or duplicated.

### Area-CFL timestep reduction follow-up

To distinguish a too-large area-averaged timestep from a broader method
failure, job **438549** was submitted on one A40. It restarts from the area
branch's exact last-valid t=2.000107644 checkpoint and targets t=2.2, just past
the prior t=2.145 k3 positivity failure. It reuses the same area-CFL binary,
source, W=0/fixed-mass state, and force settings; the only parameter change is
`Courant_hydro=0.4→0.2`. The runner aborts on nonpositive internal energy and
will write the existing 15-check audit on successful arrival. Its immutable
inputs and source overlay passed SHA-256 preflight, and `bash -n` accepted the
runner. Run receipt and script are in the new GPFS run directory
`area_cfl_c020_from_t2000_to220_a40_20261028`; the audit will be written to
[the planned t=2.2 receipt](KEPLER_Q1_AREA_CFL_C020_T2P2_LONG_GATE_AUDIT_2026-10-28.json).

The active no-cap baseline log subsequently reached `t=2.34139` with
`dEtot/|Etot0|=-6.417e-9`, zero nonpositive internal energies, and zero floor
injection; its last complete listed checkpoint is still t=2.0. This confirms
that the unchanged branch crossed the area-CFL failure time, although its
trajectory is not an exact step-by-step match after t=1.515.

At the next live log comparison, the first 85 common global step IDs
(14794–14878, through `t=2.02367`) in 438549 and the C=0.4 run had exactly
identical accepted times and dt values. The Courant change has therefore not
altered the accepted timestep in this prefix; a non-face bound likely wins,
but its identity is not recorded. Keep 438549 running toward the prior failure
time to see whether the acoustic area bound becomes active there. Any claim
that C=0.2 stabilizes the branch requires a timestep divergence before the
previous failure point.

The follow-up 438549 remained exactly matched to the C=0.4 branch through
187 common steps (`t=2.05742`) in its last comparison; it had not yet reached
the prior failure. Because the Courant change was inactive over that prefix,
a second independent diagnostic branch, **438626**, was submitted from the
same exact t=2.000107644 checkpoint. It keeps `Courant_hydro=0.4`, lowers only
`GFS_VARPRESS_GEOM_CFL=0.1→0.05`, and enables `[GPU_AREA_CFL_APPLIED]` logging
to separate the area-acoustic and non-face bounds. It targets t=2.2 with the
same positivity abort and long-gate audit. Static input/source hashes and
`bash -n` passed; no existing job was cancelled. Its receipt and runner are
in `area_cfl_geomcfl005_from_t2000_to220_a40_20261028`.

At the failed t=2.145 k3 stage, cell 92210's predicted volume is 1.3875 times
its accepted start volume while integrated internal energy changes from
`9.775e-10` to `-5.944e-10`. A gamma=5/3 adiabatic expansion by this volume
ratio would instead leave roughly `7.858e-10`. This indicates a substantial
stage-predictor positivity violation, not roundoff-level negativity; the
adiabatic estimate is diagnostic, not a claim about the exact local solution.

### Limiter telemetry and C=0.2 branch correction

The first telemetry sample from 438626 at t≈2.00074 reports
`area_dt_min=0.015431` and `nonface_dt_min=0.00020456`, with
`selector=nonface`. In the same first steps, changing only
`GFS_VARPRESS_GEOM_CFL` from 0.1 to 0.05 reduces accepted dt from about
`4.11e-4` to `2.10e-4`. Thus the area-acoustic candidate is not controlling
near t=2; the geometry/non-face limit is active and its coefficient directly
changes the accepted step.

Correction to the preceding status note: job 438549 was subsequently
cancelled at t=2.06582 after 213 common global step IDs had identical accepted
times and dt values to the C=0.4 path. It produced no t=2.2 endpoint and no
additional restart; the t=2.0001 checkpoint remains preserved. This was our
own redundant job, stopped after evidence showed its Courant change did not
affect the accepted timestep. Job 438626 is the active timestep-reduction
branch; the original no-cap and A100 clip-swap jobs were not cancelled.

### Latest live Q1 states

The last syntax scheduler poll shows baseline 438092, clip-swap 438416, and
geometry-CFL diagnostic 438626 RUNNING; 438549 is absent after our cancellation.
Their latest shared logs are at t=2.36968, t=2.33576, and t=2.01120,
respectively. The baseline and clip-swap ledgers have relative energy drifts
`-6.431e-9` and `-7.425e-10`, with no nonpositive internal energy and zero
floor injection. Job 438626's latest ledger is `-3.637e-12`, also with no
nonpositive internal energy or floor injection; its area/non-face telemetry
still selects non-face (`0.01523` versus `0.0002063`). It has not yet reached
the prior t=2.145 failure point. H200 MPI2 438394 and sparse continuation
438533 remain `QOSGrpGRES`-pending; dependency jobs 438168, 438364, 438172, and
438420 remain pending. No baseline or clip-swap trajectory was changed.

### Execution recovery and resumed evolution (2026-10-28 evening)

The default shell failure was traced to a missing Codex sandbox executable;
the reviewed host-execution path works. The [recovery record](KEPLER_EXECUTION_RECOVERY_2026-10-28.md)
supersedes the preceding live status: 438626 completed exact t=2.2 and passed
all 15 checks; 438092 and 438416 subsequently stopped at the identical
t=2.654918 halo-only guard. Their accepted step-start snapshots were retained
and independently audited. Recovery 438957 is running from that accepted
state on A100 with only halo 28→40 changed, and has crossed the failed stage.
Validated segments submit their next t=3→4→5→6→6.5→7→8→9→10 continuation
automatically. Independent A100 MPI pair 438959 initially waited on the normal
QoS total eight-GPU user cap, then started on syn102. At 22:35 both jobs were
RUNNING; recovery 438957 reached t=2.66689 with zero halo-clipped cells and
zero recorded energy injection. Neither t>6 nor the GIZMO long endpoint is
yet verified.

### User-requested pending-job cleanup (2026-10-28 evening)

The [verified cleanup receipt](SLURM_CLEANUP_2026-10-28_225639/README.md)
records 34 obsolete/deferred pending Kepler jobs cancelled with zero elapsed
allocation time. These include the failed-parent q2 chain, old MPI/KDK
variants and unstarted cap16 pair/audit. Preserve their files and revisit a
specific experiment only when current results require it. Remaining Kepler
jobs are 438957, 438959 and 437772 RUNNING, plus the immediate refined-control
continuation 437779 PENDING. Other projects and in-flight calculations were
preserved. The 384² long trajectory still submits its next segment only after
the native and radial acceptance gates pass.
