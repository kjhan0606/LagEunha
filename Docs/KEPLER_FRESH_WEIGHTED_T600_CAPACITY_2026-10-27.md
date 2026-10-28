# Fresh weighted Kepler t=600 capacity from measured native A10 throughput

The accepted [t=6.5 native run](KEPLER_Q2_FRESH_WEIGHTED_NATIVE_T65_REVIEW_2026-10-27.md)
advanced from `t=0.825` in **2,926 RK4 steps**. Its own Bash timing was
`4,815.824 s` (Slurm elapsed `1:20:22`), so the measured mean is
`1.646 wall s/accepted step` at mean `dt=0.0019395`. This is the actual
55,351-site CUDA+MPI+OpenMP weighted solver on one A10, not a Python gas
simulation. A current RK4 step tessellates four times; the recent native
log reports roughly `0.16–0.18 s` per tessellation call. The remaining
wall time includes host work, forces, communication, logging, and I/O;
those costs are not isolated by the tessellation timer.

During monitoring of the live `t=10→20` continuation **437531**, at Slurm
elapsed `47:04` the largest logged time was `t=13.3365`, step 1721. The
corresponding largest RK4 ledger step was 1720, with relative total-energy
error `−2.69e−9` and zero floor or dual-energy injection. This is a progress
sample, not an endpoint audit or an isolated benchmark; another A10 job was
running on the same node, and merged stdout/stderr lines can appear out of
order.

If both the step size and measured throughput stayed fixed from `t=6.5`
to the [GIZMO Figure 7 time `t=600`](https://arxiv.org/pdf/1409.7395),
approximately **306,000 more steps** and **140 A10 wall hours** would be
required. This is a planning extrapolation, not a scheduled completion
time: cell topology, timestep, hardware contention, dumps, physics failure,
and code changes can alter it. The A10 partition's 48-hour job limit means
one segment cannot cover the whole remaining interval at this measured
speed; restarting from audited checkpoints is essential.

The isolated [v5/v10 same-GPU short pair](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_wallselect_cache_pair_v1/speed_audit.json)
kept ten physical fields bitwise equal for both W=0 and weighted controls.
For its weighted 13-step replay, wall time fell from `25.185` to `21.579 s`
(ratio `0.857`, about 14% faster). The later matched weighted replay from
the accepted `t=6.5` checkpoint also passed: [job 437538](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_wallselect_weighted_late_pair_v1/pair_audit.json)
completed both 21-check audits and found v5/v10 bitwise identical across ten
fields through `t=6.525`; whole-run time was `24.159→20.867 s` (13.6% lower
for v10). These are short A10 results, not proof of long-time equivalence.

The accepted v5 weighted chain completed `t=10→20` on A10 as job **437531**,
then `t=20→30` on H200 as job **437613**. The t=30 receipt passes all 21
configured checks. In the matched v10/v11 `t=30→30.025` comparison, job
**437794** completed both native runs and both 21-check audits, but its final
pair auditor failed to import a helper. The saved outputs were preserved; the
comparison was rerun with the helper directory on `PYTHONPATH`. All four pair
checks passed, including bitwise equality of the audited fields at the initial,
parent, and endpoint states. Its one short runtime sample was 15.019 s for v10
and 11.451 s for v11 (23.8% lower wall time), which is not yet sustained
throughput evidence. The recovery reused the saved native outputs and did not
repeat either solver run.

The guarded v11 `t=30.025→40` continuation is staged as job **437812**; it
started on H200 at 21:05 KST. Its matched 1-rank/2-rank probe is job **437813**,
dependent on 437812. Neither binary should be projected to `t=600` as a
sustained speedup until longer native intervals pass. Even if the previously measured 13.6% v10
gain persisted, it would reduce the A10 estimate only from about 140 to 121
hours; additional measured acceleration is still needed. The H200 pair and
MPI comparison are the next throughput evidence.
No direct same-IC AREPO timing is available, so this record does not claim
that either code meets or misses the user's AREPO speed target.

The measured A10 baseline would need about **2.9×** throughput to finish
`6.5→600` in one 48-hour segment, though multiple guarded segments are
allowed. Reducing global
RK4 tessellation count or using individual timesteps could change the cost
substantially, but those integrator changes require their own conservation
and restart gates and should not be inferred safe from the present short
wall-select cache test. The next prerequisites are the v10/v11 gate through
`t=30.025`, the v11 gate through `t=40`, and the matched MPI rank comparison.

**Measured H200 scale from job 437613:** the accepted v5 state advanced from
`t=20` to `t=30` in 6,471 s, or 647.1 wall seconds per simulated time unit.
Holding that rate fixed from `t=40.25` to `t=600` gives about **100.6 H200
hours**. Three guarded segments of `200`, `200`, and `159.75` time units
would each fit under the 48-hour job limit at that rate. This is a v5
throughput scale estimate; the t>40 step size, v11 timing, H200 contention,
and long-time physics may change it. The v10/v11 and MPI comparisons must
provide current rates before choosing the production rank count.

**Execution correction and recovery (2026-10-27):** the first t=30 v10/v11
pair job **437623** stopped in its restart-staging Python before either
simulation because a local `root` variable was undefined. The fresh v2 pair
job **437794** then completed both simulations but failed only in the final
comparison import. The comparison was recovered from the saved receipts and
snapshots. Its blocked t=40 and MPI dependents (**437796**, **437798**) were
canceled before allocation; replacements **437812** and **437813** preserve
the endpoint gates, with 437813 after 437812. These packaging repairs did not
change solver code or the physical acceptance criteria.

For the future segmented continuation, `EUNHA_DUMP_DT` is a physical-time
interval (confirmed in `src/Exam/KH/kh.c`). A completed 55,351-cell snapshot
is 31,882,244 bytes. Using a 5-unit interval over planned segments
`40.25→240.25→440.25→600`, including segment-end snapshots, would produce
about 114 files or **3.63 GB** of particle snapshots. Extending the current
0.825-unit cadence would produce roughly 680 files and about **21.7 GB**.
The current GPFS scratch mount reports 222 TB available, so the 5-unit cadence
keeps restart points roughly 54 minutes apart at the measured v5 throughput
while reducing snapshot volume by about sixfold. Apply this cadence when the
post-MPI continuation scripts are staged; the current 0.825-unit runs remain
unchanged.

**Preliminary v11 H200 throughput sample (2026-10-27 21:11 KST):** job
**437812** advanced from its accepted `t=30.025` restart to `t=31.1714` after
6:34 Slurm elapsed, or about **344 wall seconds per simulated time unit** over
this early prefix, including restart checks and startup. The live RK4 ledger
at step 596 reports `dt=0.00193939`, relative total-energy drift
`−9.007e−10`, and zero floor/source-floor injection. If that wall rate held,
`40.25→600` would take about **53.5 hours**, with `200`, `200`, and `159.75`
unit segments taking roughly `19.1`, `19.1`, and `15.3` hours. Across 2,454
GPU weight-solve stages, the median was 4.881 ms, p99 8.039 ms, and maximum
0.381 s (15 stages exceeded 0.1 s). Thus weight solve time is usually small
but has rare long stalls. These are early-prefix measurements, not production
estimates; the completed `t=40` receipt and MPI rank comparison must supply
sustained throughput and parallel-efficiency evidence before choosing the
long-run layout.

**Live v11 H200 checkpoint (2026-10-27 21:17 KST):** job **437812** was still
RUNNING after 12:16 Slurm elapsed and reached `t=32.1643`, global step 1,119,
from its accepted `t=30.025` restart. The latest RK4 ledger reports
`dEtot/|Etot0|=-1.695e-9`, `dEhyd/Ehyd0=-7.318e-6`, no nonpositive internal
energy, and zero floor/source-floor injection. The early-prefix rate is about
344 wall seconds per simulated time unit, including startup. If it held over
`40.25→600`, that interval would take about 53.5 hours (three segments of
about 19.1, 19.1, and 15.3 hours). This remains a rough projection; the `t=40`
physics receipt and dependent MPI rank comparison **437813** are still
required before selecting the long-run layout. Across 4,428 logged stages,
GPU weight-solve time has median 5.078 ms, p95 6.603 ms, p99 8.377 ms, and
maximum 380.723 ms (24 stages exceeded 100 ms). Tessellation time has median
68.064 ms and p99 72.890 ms. These timings do not yet establish sustained
throughput or long-time physical accuracy.

**Live v11 H200 update (2026-10-27 21:23 KST):** at the 21:23:04 Slurm
check, job **437812** remained RUNNING with 18:01 elapsed and its log at
`t=33.13`, step 1,617. The latest ledger reports
`dEtot/|Etot0|=-2.471e-9`, `dEhyd/Ehyd0=7.326e-7`, and zero
floor/source-floor injection. The full-prefix measured rate is now about
348 wall seconds per simulated time unit, so a fixed-rate projection from
`40.25→600` is about **54.1 hours** (approximately 19.3, 19.3, and 15.5 hours
for segments of 200, 200, and 159.75 units). This remains only a capacity
estimate; the native `t=40` receipt and MPI comparison **437813** have not
completed.

Across 7,570 logged v11 GPU stages, weight-solve time has median 4.989 ms,
p95 6.450 ms, p99 8.362 ms, and maximum 380.723 ms (41 stages exceeded
100 ms). Tessellation time has median 69.858 ms and p99 89.365 ms. The latest
100 CPU-step records show a 0.50 s median step with `update=0.30 s` and
`force=0.10 s`; these single-GPU measurements identify costs but do not
establish the performance of the future MPI layout.

**Live v11 H200 update (2026-10-27 21:31 KST):** job **437812** remained
RUNNING at the 21:31:25 queue check, with 26:22 elapsed and the log at
`t=34.4211`, step 2,325. Its latest ledger reports
`dEtot/|Etot0|=-3.383e-9`, `dEhyd/Ehyd0=-5.395e-6`, and zero floor/source-floor
injection. The full-prefix rate is now about 360 wall seconds per simulated
time unit. If it held, `40.25→600` would take about **56.0 hours**, or roughly
20.0, 20.0, and 16.0 hours for the planned segments. This still requires the
exact `t=40` gate and dependent MPI scaling job **437813**; it is not a
validated production runtime or proof of `t=600` stability.

**Live v11 H200 update (2026-10-27 21:34 KST):** the 21:34:51 poll still
showed **437812** RUNNING after 29:48 elapsed; its log had reached `t=34.8677`,
step 2,617. The latest ledger reports `dEtot/|Etot0|=-3.591e-9`,
`dEhyd/Ehyd0=-1.169e-7`, and zero floor/source-floor injection. Dividing the
full elapsed time by the advance from `t=30.025` gives about **369 wall
seconds per simulated time unit**; a fixed-rate estimate for `40.25→600` is
about **57.4 hours** (three segments of about 20.5, 20.5, and 16.4 hours).
This is still an evolving one-rank sample; the t=40 gate and MPI comparison
remain required.

## Actual t=40 gate result

The v11 continuation **437812** subsequently reached exact `t=40` but ended
`FAILED 1:0`: the independent Laguerre audit found a hidden active disk
power cell. Initial site ID 30046 at `r0=0.508232999` has positive stored
volume `1.40004e-4` and weight `7.99359e-5`, but no lower-hull power-cell
facets. A second inner-ring site has a 35.3% stored-volume/power-area
mismatch, and the maximum neighboring `|Δw|/d²=1.06252`. The disk density
equilibrium error is `2.48189e-3`, also above its `1e-3` gate. Full diagnosis
and the preserved endpoint path are recorded in
[the segmented-run report](KEPLER_Q2_HALO28_SEGMENTED_T65_2026-10-27.md).

Uniformly scaling the saved weights is not a valid repair: a 0.75 scale
brings the measured face ratio below 0.8 but leaves the cell hidden; a 0.5
scale restores visibility but creates an 8.91 maximum stored-volume/power-area
mismatch. Thus the previous t=600 timing projections are capacity estimates
only. The accepted continuation still ends at `t=30.025`; do not promote the
failed `t=40` snapshot as a restart parent. A corrected candidate must keep
every active power cell visible and its stored volume consistent with the
power geometry, then pass the same exact t=40 audit before extending the
long-time projection. No solver source was changed in this update.
