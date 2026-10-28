# Full q2 halo-28 continuation through the historical t>6 interval

## Why the pending second interval was split

The accepted full-q2 H200:1, CUDA+MPI+OpenMP job **437076** is evolving
from exact `t=.825` to `t=2.221` with a 28-cell periodic halo. Near
`t=.835`, its measured accepted step is `7.64202e-5`; a live 21-step
window used 53 seconds, or `2.52 s/step`. Keeping both values constant
would make the original `2.221→6.5` single interval about **39.2 hours**
against Slurm's 48-hour limit. This is a capacity estimate, not a
performance promise: evolving geometry or a smaller timestep could erase
the margin. Its `0.825` dump interval could also force several hours of
replay after a timeout. The present `t=2.221` job was left untouched.

The not-yet-started `t=6.5` job **437077** was replaced by two physically
gated intervals. This adds a useful check at `t=4.5` and reduces the
projected compute time of each new job to roughly 21 and 18 hours at the
sampled rate. Both have a 36-hour cap and write intermediate checkpoints
every `0.5` time unit. The dump schedule is an output setting; the EOS,
force, RK4, fixed mass, W=0, halo 28 and timestep controls remain those
of the accepted q2 run. These new jobs have not yet produced physical
results.

The exact stop at `t=4.5` clips one RK4 step and a restart can change the
subsequent step sequence, so this segmented path is not claimed to be
bitwise identical to one uninterrupted `2.221→6.5` run. The predecessor
and endpoint audits must establish physical accuracy on the actual path.

The original runtime projection used a 2.52 s/step sample. The latest
100-step H200 profile gives 2.895 s/step at `dt=7.64202e-5`. Holding those
values fixed estimates about 24.0 h for `2.221→4.5` and 21.1 h for
`4.5→6.5`, within the 36 h caps. This excludes checkpoint overhead and any
later timestep reduction, so each segment still needs live walltime review.

| Job | Exact interval | Required predecessor | Audit outputs |
|---|---|---|---|
| **437250** | `2.221→4.5` | `afterok:437076` plus passed `audit_t2p221.json` | `audit_t4p5.json`, `audit_radial_t4p5.json` |
| **437251** | `4.5→6.5` | `afterok:437250` plus passed `audit_t4p5.json` | `audit_t6p5.json`, `audit_radial_t6p5.json` |

The independent sparse `t=10` and q1 arithmetic probes have no artificial
dependency on this full-q2 chain. They remain eligible to run concurrently
if GPU quota permits; the present `QOSGrpGRES` wait does not justify
serializing independent science jobs in their Slurm dependency graph.

The [single segment runner](kepler_polar_q2_halo28_segment_20261027.slurm)
uses a fixed two-target argument contract. The
[bundle stager](stage_kepler_polar_q2_halo28_segment_20261027.py)
copied and SHA-verified its script and four auditor/restart dependencies
into new directories beneath the current `to2p221_h28` work directory.
Both bundles passed `bash -n`, full `sha256sum -c`, and Python import
checks before submission; the runner SHA-256 is
`506966cbb4a598c0c116b2253f47d1893281c58d6c149991278d50f326f8e899`.
At run time each job checks the predecessor's exact-time complete
checkpoint, all preceding physics checks, zero injection, accepted binary
SHA-256, and `params.dat` SHA-256
`2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43`
before copying state. The native endpoint must pass the full 15-check
orbit/energy/entropy/area gate and the radial state audit before the next
Slurm dependency can succeed. No failed RK4 stage is reused.

After both new jobs and their `afterok` dependencies were verified in
Slurm, the superseded pending job **437077** was canceled without
execution (`CANCELLED`, elapsed `00:00:00`). The live `437076` remained
`RUNNING`; no completed interval was discarded or recomputed.

These are full-count periodic q2 controls. A successful `t=6.5` endpoint
would cross the historical time barrier, but the GIZMO `t=600` image and the
paper's separate statement that MFM was evolved for approximately 1,000
orbits remain separate gates. The paper does not specify the radius or exact
time for that latter statement in the sentence, so it must not be silently
equated with the `t=600` figure.

## Latest live checkpoint

At the latest verified Slurm read, job **437076** was still RUNNING on
`syn104`; its log had reached `t=1.95547`, global step 26248. The nearby
complete RK4 ledger at `t=1.95290` reports relative total-energy drift
`-2.523e-10`, `n_ie_le0=0`, and zero floor injection. A low-rate live snapshot
audit of the already written `t=1.500026676` checkpoint (CPU job **437570**,
`khout.019926.dat`) passes fixed mass, W=0, positivity, initial and final
area tiling, `rho=m/V`, and no-injection checks. Its disk mass-weighted
density and absolute entropy-log changes are `1.33690e-4` and `2.39716e-4`;
orbit-position error is `9.4571e-9`, disk angular-momentum drift
`-4.7701e-13`, and total-energy drift `-2.13734e-10`. This is an accepted
intermediate checkpoint, not a `t>6` result. The [density-field figure](KEPLER_Q2_H28_DENSITY_T150_2026-10-27_corrected.png)
shows native surface density and relative density change on the common GIZMO
`[-2,2]` view and color scales. The `2.221` endpoint has not yet been written
or audited; jobs 437250 and 437251 remain gated by that receipt and the
subsequent `t=4.5` receipt respectively.

The published-annulus audit for `0.8<R<1.2` gives unweighted
`L1_velocity=5.0849e-9`, literal-paper `L1_density=0.01000005`, and
plateau-corrected `L1_density=1.8124e-4` at this endpoint. The density values
are reported both against the paper's unit reference and the run's actual
initial plateau `1.01`; these early values are below the Figure 8 thresholds,
but this is not the paper's `10`-orbit measurement or its `t=600` Figure 7
comparison. The full [metric receipt](KEPLER_Q2_H28_GIZMO_METRICS_T150_2026-10-27.json)
keeps the unweighted paper norm separate from the image's mass-weighted disk
change.

## Current live observation (2026-10-27)

Job **437076** is still RUNNING on H200. Its latest log sample reached
`t=2.09442`, global step 28,238, with `dt=7.64202e-5`. The latest printed
RK4 energy ledger has `dEtot/|Etot0|=-2.975e-10`,
`dEhyd/Ehyd0=3.917e-4`, `n_ie_le0=0`, and zero
floor, dual-energy, or source-floor injection. The exact `t=2.221` endpoint
and its audit are not complete. The last Slurm check still showed 437076
RUNNING; q2 continuation jobs 437250/437251 and later t=10/20 jobs 437639/
437640 remain gated on the endpoint receipts.

## t=2.000016 intermediate checkpoint audit (2026-10-27)

The already-complete native snapshot `khout.026880.dat` was audited while
437076 continued toward t=2.221. The strict long-gate receipt
[`audit.json`](/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/halo12_after_failure_to0825/to2p221_h28/live_t200_audit_cpu_v1/audit.json)
passes all 15 checks: finite state, exact fixed mass and W=0, positive
internal energy and volume, both 6x6 area tilings, `rho=m/V`, no floor/dual
injection, orbit and conservation tolerances, and the two paper-norm bounds.
It reports disk (`0.5<=r0<2`) mass-weighted absolute density change
`1.48529e-4`, absolute/signed entropy-log changes `3.19862e-4`/`-7.41898e-5`,
radial-velocity error `1.48406e-8`, orbit-position error `1.89678e-8`, disk
angular-momentum drift `-5.46148e-13`, and total-energy drift `-3.27844e-10`.
The paper annulus has `L1_velocity=7.07984e-9` and plateau-corrected
`L1_density=1.74547e-4` (literal unit-reference norm `0.01000004`). The
[radial receipt](/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/halo12_after_failure_to0825/to2p221_h28/live_t200_audit_radial_repair_v1/radial.json)
also passes its five checks.

Job 437790 wrote the passing long-gate receipt but exited 1 when its optional
radial step lacked the imported `audit_kepler_seam_checkpoint_20261026.py`.
The preserved checkpoint was not rerun; CPU job 437791 completed only the
radial audit from a separate directory with the missing helper staged. This
was a packaging failure, not a simulation or physics-gate failure.

The same radial receipt localizes much larger distortion to the exterior
background: `r0>=2.75` has mass-weighted density change `0.3120` and absolute
radial speed `0.1658`; cells currently within 0.2 of a periodic edge have
angular-momentum change `-0.6013`. The `0.5<=r0<2` disk remains accurate at
this checkpoint, but the outer seam remains a separate possible source of
later contamination. This t=2 receipt is not the pending t=2.221 endpoint,
does not cross t=6, and says nothing yet about survival to t=600.

Compared with the accepted t=1.500026676 snapshot, disk density error rises
from `1.33690e-4` to `1.48529e-4` and absolute entropy-log error from
`2.39716e-4` to `3.19862e-4`; signed entropy change remains negative
(`-4.90122e-5` to `-7.41898e-5`). This is a gradual increase across these two
checkpoints, not evidence of the earlier runaway heating failure. The t=4.5
and t=6.5 audits must establish whether this trend remains controlled.

## GIZMO figure metrics and shell spacing at t=2 (2026-10-27)

The [paper Fig. 8 metric receipt](KEPLER_POLAR_Q2_GIZMO_METRICS_T200_2026-10-27.json)
recomputes its `0.8<R<1.2` annulus norms on `khout.026880.dat` from the q2
`t=0` snapshot. It gives `L1_velocity=7.07984e-9`, literal-paper
`L1_density=0.01000004`, and density error relative to the actual `1.01`
plateau of `1.74547e-4`. These are below the Fig. 8 thresholds (`0.01` for
velocity and `0.25` for density), but `t=2` is only about 0.32 orbital times
at `R=1` and 0.90 at `R=0.5`. It is an early checkpoint, not a long-time pass.

The [ring-spacing audit](KEPLER_Q2_RING_SPACING_AUDIT_2026-10-27.json) traces
the initial 68 radial shells through `t=2`. In the scientific disk the median
nearest same-shell generator distance is `0.01087755`; the median nearest
adjacent-shell distance is `0.02218841`. No disk cell has a closer adjacent
shell candidate. This does not support an initial inter-shell crowding cause
through `t=2`; the next matched check is required at `t>=6.5`.

Hopkins' [GIZMO paper, §4.2.4 and Figures 6–8](https://arxiv.org/abs/1409.7395)
uses an effective `256^2` spacing across `[-2,2]^2`; q2 uses `384^2` across
`[-3,3]^2`, giving the same nominal spacing `0.015625` in the active disk.
The initial layouts are not identical: the paper removes particle-method
sites inside `r<0.5` and outside `r>2`, while its mesh comparison uses a
low-density floor and smoothed disk edges. Thus the current density images
are resolution-matched comparisons, not exact same-site reproductions. The
paper's Figure 7 is the exact visual endpoint `t=600`; its caption says the
`r≈0.5` inner disk has completed more than 250 orbits by then. The current q2
[t=2 density field](/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/halo12_after_failure_to0825/to2p221_h28/kepler_density_t2.000016_w0_h200_20261027.png)
is retained for visual tracking, while the `t=600` density comparison remains
the required figure-level gate.

## Primary-source GIZMO endpoint audit (2026-10-27)

The paper's §4.2.4 defines a 2D `γ=5/3` disk with `Σ=0.01+(r/0.5)^3`
inside `r=0.5`, a unit plateau for `0.5≤r≤2`, and an outer profile
`Σ=0.01+[1+(r-2)/0.1]^-3`; the pressure is `10^-6`, self-gravity is off,
and the softened Kepler potential and circular speed use the same softening
parameter. Its text says the qualitative behavior is unchanged for any small
softening and edge values, but does not give the numeric softening in this
section. This supports the present low-pressure, fixed-mass W=0 baseline as a
physics comparison, while the exact site distribution and softening still
need to be recorded when claiming an exact setup match.

The paper's Figure 8 uses the annulus `0.8<R<1.2`, measures velocity error
at `t=10` orbits at `R=1` (`t=20π` in the paper's units), and plots the first
crossings of `L1_v=0.01` and `L1_ρ=0.25`. The current early q2 receipts compute
these same annulus norms, but their `t≈2` values are not the Figure 8
measurement. The paper's text additionally reports MFM integration to about
1,000 orbits without specifying the radius or corresponding exact time in
that sentence. Therefore the audit plan keeps distinct gates: the exact
`t=600` Figure 7 density image, the Figure 8 norm checkpoint at `t=20π`, and
resolution of the approximately 1,000-orbit statement before claiming the
paper's longest reported duration has been matched.

## Latest live q2 observation (2026-10-27)

At 20:50 KST, job **437076** was still RUNNING on `syn104`; its live log had
reached `t=2.15906`, global step 29,168, with `dt=1.91481e-5`. The latest RK4
ledger had `dEtot/|Etot0|=-3.041e-10`, `dEhyd/Ehyd0=4.246e-4`, no nonpositive
internal energy, and zero floor/source-floor injection. The preceding global
time marker at `t=2.15904` had recovered to `dt=7.64202e-5`; intermittent
one-step minima recur, but this sample does not show sustained timestep
collapse. The exact `t=2.221`
endpoint has not yet been written; dependent `t=4.5` and `t=6.5` jobs remain
gated on its complete physics and radial receipts.

The active run's parameter file sets `KH Kappa=0`; in this source,
`getw2forHydroParticle()` returns exactly `w2=0` when `GAS_Kappa==0`
(`src/Exam/exam.c`). The run script also leaves dynamic `w2` controls and
`GFS_LAGUERRE_ROTATION` unset. Thus this is the zero-weight Voronoi baseline:
the Laguerre power weights are zero, and no weight evolution or cell-spin
state is active. The preserved `t=2.000016` native audit receipt independently
passes both `mass_exact` and `W0_exact`; the cell masses are held fixed by the
mass-update path, which is a separate condition from setting `w2=0`.

## Live status and execution profile (2026-10-27 20:54 KST)

The latest log sample for **437076** reached `t=2.16398`, global step 29,242.
Its accepted `dt` recovered to `7.64202e-5` after the isolated smaller step
seen near `t=2.1614`; this supports intermittent low-`dt` events, not a
sustained timestep collapse. The accompanying RK4 ledger reports
`dEtot/|Etot0|=-3.053e-10`, `dEhyd/Ehyd0=4.271e-4`, no nonpositive internal
energy, and zero floor/source-floor injection. The last Slurm query at 20:51
showed 437076 RUNNING; the log continued updating afterward. Its exact
`t=2.221` endpoint and gate are still outstanding. The dependent H200 chain is
already staged as `437076 → 437250 (t=4.5) → 437251 (t=6.5) → 437639 (t=10)
→ 437640 (t=20)`, with `afterok` dependencies.

The same step telemetry gives a concrete performance lead. On 437076, a recent
full step took about 3.0 s in the CPU step path (`update=1.6 s`, `force=0.9 s`);
the corresponding CUDA telemetry reported `SoA=116 ms`, `kernel=1.7 ms`, and
`writeback=8 ms`, while recent tessellation calls were about 300 ms each.
This points to the host-side mesh/update and data-movement path as the first
performance target; it does not support spending effort first on the CUDA
force kernel alone. A separate A10 W0 refined-420 run, **437771**, reached
`t=0.239584` at 20:54; it is a separate branch and is not a substitute for
the full-q2 H200 chain.

## Live checkpoint (2026-10-27 21:11 KST)

Job **437076** advanced to `t=2.18322`, global step 29,520. Its latest RK4
ledger reports `dEtot/|Etot0|=-3.119e-10`, `dEhyd/Ehyd0=4.372e-4`, no
nonpositive internal energy, and zero floor/source-floor injection. A nearby
smaller step recovered to `7.64202e-5`; the `t=2.221` stop and audit remain
pending. The separate A10 refined-420 W0 branch **437771** reached
`t=0.264115` at its last confirmed log sample;
it is not a continuation of the full-q2 H200 state.

## Live checkpoint (2026-10-27 21:17 KST)

Slurm confirms the W=0 full-q2 H200 job **437076** is still RUNNING. Its log
has reached `t=2.18991`, global step 29,614. The latest RK4 ledger reports
`dEtot/|Etot0|=-3.134e-10`, `dEhyd/Ehyd0=4.409e-4`, no nonpositive internal
energy, and zero floor/source-floor injection. At the observed step, the CPU
step path is about 3.0 s and the four GPU tessellations are about 0.31–0.33 s
each; the force kernel remains about 1.7 ms. The exact `t=2.221` endpoint and
native audit remain outstanding, so the dependent `t=4.5→6.5→10→20` chain has
not started.

The separate W=0 q1 H200 control **437212** is RUNNING at `t=1.1244`, global
step 14,861; its latest ledger has `dEtot/|Etot0|=-4.578e-11`,
`dEhyd/Ehyd0=5.542e-5`, and zero floor/source-floor injection. The refined
420² A10 W=0 branch **437771** is RUNNING at `t=0.279475`; its latest ledger
has `dEtot/|Etot0|=1.325e-11`, `dEhyd/Ehyd0=6.629e-6`, and zero floor/source
floor injection. These are live-prefix conservation checks, not endpoint
gates or evidence of `t>6` stability.

## Live checkpoint (2026-10-27 21:23 KST)

At the 21:23:04 Slurm check, **437076**, **437212**, **437771**, and **437812**
were all RUNNING. The full-q2 W=0 job **437076** reached `t=2.19632`, global
step 29,708. Its latest ledger reports `dEtot/|Etot0|=-3.167e-10`,
`dEhyd/Ehyd0=4.442e-4`, no nonpositive internal energy, and zero
floor/source-floor injection. Over its latest 100 CPU-step records, total
step time has median 3.00 s and mean 3.039 s (`update` 1.60 s, `force`
0.90 s); 400 corresponding GPU tessellation calls have median 314.3 ms.
This run remains OMP=2. The force kernel is not the dominant measured GPU
cost.

At the same log check, q1 W=0 job **437212** reached `t=1.13408`, refined
420² A10 W=0 **437771** reached `t=0.28918`, and weighted v11 H200 **437812**
reached `t=33.13`. These branches are still running and their exact endpoint
receipts have not completed. The q2 `2.221→4.5` and `4.5→6.5` segments each
have a 36-hour limit; using the current 3.0 s/step and `dt=7.64202e-5`
sample projects about 25.1 h and 22.0 h respectively, before startup,
checkpoint, or timestep-change overhead. The 4-rank/4-GPU job **436751** is
still pending QOS; it uses a 256² grid on a 24-by-24 domain with `Rin=2`,
`Rout=8`, and `epsilon=.047`, and stops at `t=.003`. This differs from the
384² W=0 run (`Rin=.5`, `Rout=2`, `epsilon=.01175`), so even if it runs it
cannot establish scaling for the full-q2 production branch.

## Matched OMP scaling probe staged (2026-10-27 21:30 KST)

To measure the actual full-q2 bottleneck, a short 384² W=0 H200 pair was
staged from the exact audited `t=2.221` parent. Both runs target `t=2.231`,
use the same binary, parameters, one MPI rank, one H200, and the same physics
environment; only OpenMP threads and the matching CPU allocation differ:
**437835** uses OMP=2/2 CPUs and is `afterok:437076`, while **437836** uses
OMP=8/8 CPUs and is `afterok:437835`. Each runner checks the native endpoint
and radial audits. Both jobs have a low Slurm priority (`nice=10000`), so the
already-staged `t=4.5` continuation remains ahead of this optimization probe.
No solver source or active run was changed.

## Live endpoint approach (2026-10-27 21:31 KST)

The 21:31:25 Slurm check still showed **437076** RUNNING. Its log reached
`t=2.20564`, global step 29,845; the latest ledger reports
`dEtot/|Etot0|=-3.177e-10`, `dEhyd/Ehyd0=4.493e-4`, and zero floor,
dual-energy, or source-floor injection. About `0.0154` simulation time remains
to the exact `t=2.221` endpoint. Jobs **437250** (t=4.5) and **437251**
(t=6.5) remain dependency-gated. The new OMP pair **437835→437836** also
remains gated behind this accepted parent and does not replace or alter the
long continuation.

## Short timestep excursion at t≈2.208 (2026-10-27 21:33 KST)

A recent 100-marker window contains a short burst with `dt<5e-5`: the minimum
was `7.90e-7`, and the longest consecutive burst was nine steps. The log
printed `Time=2.20799` repeatedly because time is displayed to five decimals;
the accepted steps continued to advance. By `t=2.20827`, the logged timestep
had returned above `5e-5` and was approaching the usual `7.642e-5`. The RK4
ledger remained near `dEtot/|Etot0|=-3.16e-10`, `dEhyd/Ehyd0=4.505e-4`, with no
floor, dual-energy, or source-floor injection. This is a finite timestep
excursion, not a demonstrated stall or sustained collapse. Since the active
runner does not enable `GFS_VARPRESS_GEOM_AUDIT`, the limiting cell/cause is
not identified; no CFL setting was changed.

## Live checkpoint (2026-10-27 21:34 KST)

The last Slurm poll at 21:34:51 showed the W=0 q2 job **437076** still
RUNNING; its latest log reached `t=2.20936`, global step 29,903. The ledger
reports `dEtot/|Etot0|=-3.152e-10`, `dEhyd/Ehyd0=4.512e-4`, and zero floor,
dual-energy, and source-floor injection. The q1 control **437212** reached
`t=1.15404` with `dEtot/|Etot0|=-4.225e-11`, `dEhyd/Ehyd0=6.227e-5`;
the refined A10 W=0 branch **437771** reached `t=0.309126` with
`dEtot/|Etot0|=1.192e-11`, `dEhyd/Ehyd0=8.191e-6`. They all still require
their exact endpoint audits. The t=2.221 q2 gate and its dependent t=4.5,
t=6.5, and OMP-pair jobs remain unfinished.

## Verified continuation and weighted-geometry failure (2026-10-27 22:23 KST)

At the 22:23:46 KST Slurm poll on `syntax`:

- Q2 W=0 H200 segment **437076** completed at exact `t=2.221` with exit 0. Its endpoint audit is
  `/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/halo12_after_failure_to0825/to2p221_h28/audit_t2p221.json`.
  It confirms exact W=0 and fixed cell mass. On the initial-radius disk `0.5<=r<2`, mass-weighted density change is `1.40808e-4`, entropy-log change `3.56052e-4`, position error `2.46283e-8`, angular-momentum drift `-6.39371e-13`, and total-energy drift `-3.77607e-10`; floor/dual-energy/source-floor injections are all zero. This is only the `t=2.221` gate, not long-time stability.
- The continuation to `t=4.5`, job **437250**, remains PENDING with `QOSGrpGRES`; `t=6.5` job **437251** remains dependency-pending.
- Q1 W=0 H200 **437212** is RUNNING at `t=1.23665`, step 16368. Its latest ledger has `dEtot/|Etot0|=-6.148e-11`, `dEhyd/Ehyd0=8.135e-5`, and zero floor, dual-energy, and source-floor injections. A brief timestep dip (`3.62e-5`, then `2.74e-5`) recovered to `7.642e-5` on the next step.
- Refined 420² A10 W=0 **437771** is RUNNING at `t=0.388809`, step 5078; its latest ledger has `dEtot/|Etot0|=2.040e-11`, `dEhyd/Ehyd0=1.296e-5`, with zero injection.

### W64 time-dependent-weight endpoint at t=40

H200 job **437812** reached exact `t=40`, but exited FAILED because the independent Laguerre geometry gate raised `Hidden or nonpositive disk power cell`; see its preserved `audit.log` and snapshot `khout.005387.dat` in
`/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_thermalhoist_v11_fresh_to40_h200_v2/`.
The failing site is initial ID 30046 at `r0=0.508232999`, with `w2=7.99359e-5` and stored positive volume `1.40004e-4`, but it has no lower-hull power-cell facets. A second inner-ring site (ID 31092, same initial radius) has a 35.3% mismatch between stored volume and independently reconstructed power area. The maximum neighboring `|Delta w|/d^2=1.06252` exceeds the configured 0.8 guard.

Other endpoint diagnostics, computed separately because the geometry gate aborts before writing its receipt, are: disk density-equilibrium error `2.48189e-3` (above the `1e-3` gate), entropy-log error `5.10481e-5`, orbit-position error `1.67216e-5`, radial-speed RMS/circular speed `4.44651e-6`, angular-momentum drift `-1.44713e-8`, and total-energy relative drift `-1.66614e-8`. The favorable orbit/energy values do not make the invalid power geometry acceptable.

A static endpoint sensitivity calculation scaling all stored weights showed that a 0.75 scale brought the measured face ratio to 0.7969 but **did not restore the hidden inner-ring cell**. At 0.5 scale the cell reappeared, but the maximum stored-volume/power-area mismatch rose to 8.91. These are geometry-only counterfactuals, not evolved solutions. Therefore a face-ratio cap or uniform weight damping alone is insufficient; the update must also preserve every active power cell and its area/volume consistency. No solver source was changed in this session.

The dependent MPI probe **437813** was CANCELLED after its parent failed, so it will not run. A scheduler-only H100 dry run predicted 2026-10-29 23:50 KST; A40/A100 predictions were later. No alternative-GPU job was submitted. The requested `t>6` and GIZMO `t=600` stability gates remain open.


### Static target-volume re-solves at the failed t=40 positions

Two CPU-only Newton/CG geometry re-solves were run on the exact preserved
`khout.005387.dat` positions, using the same 49,038 active disk sites and the
independent all-site regular triangulation. The first targets each site's
original t=0 Voronoi volume; the second targets its stored t=40 volume. Both
start from zero weights and use a line search that rejects hidden or
nonpositive active power cells.

| Static target | Newton corrections | Hidden cells | Maximum relative area error | Weight range with inactive sites fixed at zero | Maximum face `|Δw|/d²` |
|---|---:|---:|---:|---:|---:|
| Initial t=0 volumes | 3 | 0 | `5.43e-10` | `−5.26e-5 … 3.39e-5` | `1.64` (1 face above 0.8) |
| Stored t=40 volumes | 3 | 0 | `2.20e-8` | `−5.61e-5 … 5.99e-5` | `2.31` (4 faces above 0.8) |

Both target volume sets therefore admit valid static power diagrams at the
failed generator positions. The hidden cell does not prove either target is
geometrically impossible; the evolved weight field or its native update path
failed to reach a valid representation. These geometry-only re-solves do not
repair the failed endpoint receipt or establish conservative weight
evolution. With inactive weights fixed at zero, both solutions use mixed
signs and exceed the heuristic 0.8 face-ratio guard. The earlier uniform
0.75 scaling also remained hidden despite bringing that ratio below 0.8, so
the ratio alone is neither a visibility test nor a sufficient area test.

Adding one common constant to **all** site weights leaves every power-cell
boundary unchanged. The driver should check whether a nonnegative-weight
storage convention can use this global gauge freedom instead of clipping
individual active weights. The same shift must cover active, core and
exterior sites, and the solver must depend only on weight differences. A
corrected native run still needs to restart from the accepted t=30.025 parent,
preserve stage-wise visibility and area/volume consistency, and pass the
exact t=40 physics audit.

### Snapshot replay identifies the t=40 clipped-face artifact

A snapshot-only replay was performed on the preserved t=40 state
khout.005387.dat. It used the runtime grid reported in log
(cellsize=0.0625, mx=my=98), the CUDA kernel's 3x3 candidate stencil,
the run's boxsize=5*BOXSIZE/NX=0.0868056, and the same candidate sorting
and Sutherland-Hodgman clipping. No evolution or solver source was changed.

Across 3,838,169 unique candidate pairs, the 0.05--0.95 boundary-fraction
clamp activated for three pairs touching the active disk; the largest stencil
had 158 neighbors, below the kernel's 256-neighbor cap. One clamped pair is
exactly the problematic hidden pair (30046, 31092): its signed
(w_i-w_j)/d^2=-27.52089, so its true fraction is -13.26045, but the
tessellator replaces it with 0.05. The exact power diagram has no cell for
site 30046.

For this same pair, the clipped polygon face is placed at fraction 0.05,
while the face-state interpolation and the Laguerre face-velocity anchor use
the unclamped fraction -13.26045. Their locations differ by 13.31045
generator separations, about 0.01795 in code length for this pair. Thus the
stored face area and the pressure/density/velocity and ALE face-speed
evaluation do not refer to the same geometric face; the issue is also a
moving-face conservation inconsistency, not only an area error.

Replaying those two cells with the clamp removed versus enabled gives:

| Site | Independent exact power area | Replay without clamp | Replay with clamp | Stored native volume |
|---|---:|---:|---:|---:|
| 30046 | hidden | 0 | 1.400043858e-4 | 1.400043839e-4 |
| 31092 | 4.816732569e-4 | 4.816732569e-4 | 3.118563741e-4 | 3.118563618e-4 |

The clamped replay reproduces both stored native volumes to about 4e-8
relative or better; the unclamped replay agrees with the independent power
geometry and makes the hidden cell empty. This directly identifies the
fraction clamp as the cause of these two t=40 geometry failures, rather than
merely a correlation with a large weight ratio. The independent lower-hull
audit found only two regular edges beyond the clamp threshold; the hidden
pair was not a valid power edge, which is why counting regular edges alone
missed the artifact.

The solver source should preserve the exact power half-plane during clipping
and reject or backtrack a weight update when an active site becomes hidden.
It should not manufacture a positive cell by moving an out-of-range
bisector fraction into [0.05, 0.95]. Before another long weighted run, the
driver should verify all-site area tiling and exact power-area agreement on
this saved snapshot. The run manifest records binary SHA-256
adeca5d738bce85a740b966fbbef04e3aaaba898ac7978eb5db797a620eaf5cd, but
does not record source-file hashes; the near-exact match to native volumes
strongly supports the replay, while source-to-binary provenance remains
unverified.

The common additive weight gauge discussed above preserves static geometry
only. The current face-speed path uses square roots of absolute weights and
sound-speed clamps, while the weight controls include absolute floors and
rate limits. Gauge recentering therefore needs a transformed, revalidated
weight-evolution and face-speed scheme before it can be used in production.

## Long-run W=0 continuation re-segmented and queued (2026-10-27 22:33 KST)

The W=0 control remains the zero-weight Voronoi case, with the same RK4/GPU binary and no Laguerre weight or cell-spin operator enabled. No solver source or physics settings were changed. The pending post-`t=6.5` work was split into 36-hour H200 jobs so each Slurm allocation has a realistic chance to finish:

| Job | Segment | Dependency |
|---|---|---|
| 437868 | 6.5 → 8.5 | afterok:437251 |
| 437869 | 8.5 → 10 | afterok:437868 |
| 437870 | 10 → 12 | afterok:437869 |
| 437871 | 12 → 14 | afterok:437870 |
| 437872 | 14 → 16 | afterok:437871 |
| 437873 | 16 → 18 | afterok:437872 |
| 437874 | 18 → 20 | afterok:437873 |

All seven stage directories passed `bash -n` and their `scripts.sha256` checks before submission; the runner preserves the exact-endpoint and independent radial audit gates. Each job requests one H200, two CPUs, 32 GiB, and 36 hours. The submission record is `chain_jobs_20261027.txt` beneath the `to6p5_h28_guarded_20261027` directory. The two prior unstarted jobs, 437639 and 437640, were cancelled before resubmission; neither had run or produced a checkpoint.

At submission, jobs 437868–437874 were all PENDING behind their required predecessors. The chain root 437250 (2.221 → 4.5) was still PENDING with `QOSGrpGRES`; 437251 (4.5 → 6.5) remained dependent on it. Therefore no post-6.5 evolution has started yet. Current independent controls were still RUNNING: q1 H200 job 437212 had log time `t=1.25235` (step 16577), and refined A10 W=0 job 437771 had `t=0.404524` (step 5298). These live prefixes do not pass any long-time endpoint gate.


## RTX 5000 Ada throughput screen and current queue (2026-10-27 22:43 KST)

A separate, foreground OMP=2 screen used the same accepted binary, parameters, initial reference, and audited q2 `t=2.221` restart as the pending H200 OMP probe. The binary contains `sm_89` code and ran on the RTX 5000 Ada (`compute capability 8.9`). After eight accepted steps it reached `t=2.22153` (step 30078). Its recent instrumented CPU-phase sum was 15.25 s/step median (`update=6.8 s`, `force=6.0 s`); GPU tessellation was 1.006 s/call median across 33 calls. The corresponding full-q2 H200 profile near `t=2` recorded about 2.90 s/step phase sum and 0.314 s/tessellation. These timer scopes omit RK4 stage arithmetic, so this is a throughput screen rather than a fully isolated hardware benchmark; it nonetheless shows the Ada path is roughly five times slower in the measured phase sum.

The short Ada screen showed no immediate energy failure (`dEtot/|Etot0|=-1.207e-13`, `dEhyd/Ehyd0=2.680e-7`, zero floor/dual/source-floor injection at step 8), but it did **not** reach its planned `t=2.231` endpoint and produced no endpoint receipt. It was intentionally interrupted after the timing sample because continuing the 0.01 interval was projected to take roughly half an hour and Ada is not a competitive replacement for the H200 long-run chain. The partial log and restart are preserved at
`/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/halo12_after_failure_to0825/to2p221_h28/ada_probe_t2221_to2231_20261027_omp2/`; do not count this as a physics pass.

At the 22:43:17 KST Slurm check, the `t=2.221→4.5` root **437250** remained PENDING (`QOSGrpGRES`); 437251 and 437868–437874 were dependency-pending, so no `t>6` segment had started. The q1 control **437212** remained RUNNING at `t=1.26771` and refined A10 W=0 **437771** at `t=0.420343`. The H200 queue is still the path to the long-time W=0 gate; the staged OMP=2/8 probe also remains pending and will measure thread scaling on the exact H200 branch when scheduled.


## H200 group-GRES queue diagnosis (2026-10-27 22:45 KST)

At the 22:44:20 KST check, the H200 partition consisted of `syn104` with eight H200 cards. Six one-card jobs were RUNNING on that node, leaving nominal physical capacity, while q2 root **437250** still reported `QOSGrpGRES`. The user-visible `sacctmgr` output for QOS `normal` and the `kjhan` associations showed no numeric `GrpTRES`/`MaxTRES` values, so the precise group cap is not exposed by these queries. The queue reason is therefore a group GRES limit rather than proof that all eight physical cards are allocated. No unrelated running jobs were changed; the q2 chain remains queued with its existing dependencies.


### QoS limit and alternate-partition check (2026-10-27 22:48 KST)

Follow-up Slurm accounting exposed the exact limit: partition `h200` applies QoS `part_h200`, whose `GrpTRES` is `gres/gpu=6`. At the 22:44 check, exactly six one-GPU jobs were RUNNING on `syn104`; the root job 437250 is waiting for one of those group allocations to be released. The node has eight cards, but the two spare physical cards are outside the partition QoS allowance. The same node also belongs to `sjjlee_hgpu`, but that partition only allows QoS `sjjlee_unlim`; the `kjhan` associations expose only QoS `normal`, so this is not an authorized alternate queue for this account. No queue workaround or unrelated-job cancellation was made.


### H100 alternate dry-run and q1 slot progress (2026-10-27 22:51 KST)

A `sbatch --test-only` request for one H100, two CPUs, 32 GiB, and 36 hours (no job submitted) predicted start **2026-10-30 00:20 KST**. This is substantially later than the likely release of the currently running q1 H200 control; keep the main continuation on its existing H200 dependency chain. A fresh read of q1 job **437212** reached `t=1.27982`, step 16951. This is live progress toward its `t=2.221` endpoint, not an endpoint pass. The previous scheduler poll still had root **437250** at `QOSGrpGRES`.


### Latest live check (2026-10-27 22:52 KST)

At 22:51:40 KST, root **437250** and descendants **437251**, **437868–437874** remained pending; 437250 was still blocked by `QOSGrpGRES`. q1 control **437212** remained RUNNING at `t=1.28107`, step 16968. The initially recorded remaining-time estimate was later corrected after checking the `t=0.825` restart; see the correction below. No endpoint or slot release time is guaranteed. The H200 chain remains preferable to the H100 dry-run's predicted start on 2026-10-30.


### Correction: q1 release-time estimate (2026-10-27 22:53 KST)

The previous ~3.3-hour estimate incorrectly treated the q1 job as starting at `t=0`. Its authoritative `restart_plan.json` shows `expected_time=0.825`, target `t=2.221`. At the latest log time `t=1.28298` after 4:32:34 elapsed, the average continuation rate projects about **9.3 more hours** to the q1 endpoint, roughly 08:11 KST on Oct 28 if the rate holds. This is a rough projection, not a scheduler reservation; any of the other five H200 jobs could release the group slot earlier.


## Matched W=0 MPI scaling probes staged (2026-10-27 22:57 KST)

The prior 4-GPU startup probe used a different domain and annulus, so it could not measure scaling of the full q2 branch. Two low-priority probes now use the exact accepted q2 W=0 `t=2.221` parent, binary and parameters, target `t=2.231`, OMP=2, and the same native plus radial endpoint audits: **437892** uses two MPI ranks/two H200 GPUs; **437893** uses four ranks/four GPUs and waits `afterok:437892`. Both are limited to one hour and set `nice=10000`; neither changes the long continuation.

The staging helper verified the t=2.221 receipt, fixed-mass parent gate, checkpoint, executable SHA-256 `7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284`, and parameter hash before copying inputs. Slurm rejected an `afterok` dependency on already-completed job 437076, so 437892 was submitted directly from this verified parent; the 4-rank job remains gated on 437892. Both were PENDING at submission with priority 1, and the 4-rank branch showed the expected dependency. They will not outrank or replace the t>6 H200 continuation.


## Standalone CUDA spin-accounting checks completed (2026-10-27 23:01 KST)

A10 CUDA jobs **437559** and **437550** completed with exit 0. Job 437559 reports `PASS: GPU affine spin projection, energy split and invalid rollback`; job 437550 reports `PASS: GPU nonrigid/spin/thermal pair conservation and rollback`. Their staged source and input checksums passed. These are isolated projection and two-cell ledger kernels; they do not exercise native RK4 evolution, physical face traction, CUDA/MPI halo exchange, or Kepler stability. The active W=0 q1/q2 runs still contain no intrinsic cell-spin state; the existing geometric face-speed correction is a different mechanism.

At 23:00:09 KST, q1 job **437212** remained RUNNING at `t=1.29458`; the q2 continuation and low-priority MPI probes **437892/437893** remained pending behind H200 group GRES availability. No `t>6` result is available yet.

## Live check (2026-10-27 23:07 KST)

At 23:07:13 KST, Slurm still showed q1 W=0 H200 job **437212** RUNNING on
`syn104` (elapsed 4:46:58). Its log reached `t=1.30563`, step 17,302, with
`dt=7.64202e-5`; the latest RK4 ledger reports
`dEtot/|Etot0|=-8.673e-11`, `dEhyd/Ehyd0=9.880e-5`, no nonpositive internal
energy, and zero floor/dual-energy/source-floor injection. This remains a
live prefix toward `t=2.221`, not an endpoint pass. Based on the observed
average progress since its `t=0.825` restart, roughly nine more hours may be
needed if the rate holds; this is only a queue-planning estimate.

The q2 W=0 continuation root **437250** (`2.221→4.5`) remains PENDING with
`QOSGrpGRES`; **437251** (`4.5→6.5`) and the dependent `6.5→20` chain
**437868–437874** remain dependency-pending. The exact-q2 OMP and MPI probes
also remain queued. No `t>6` evolution has started, and the q2 endpoint
receipts remain the gates for advancing the chain.

## Live check (2026-10-27 23:12 KST)

At the 23:12:22 KST Slurm check, q1 W=0 job **437212** was still RUNNING on
`syn104`, with elapsed time 4:52:07. Its log reached `t=1.314`, step 17,415;
the latest RK4 ledger reports `dEtot/|Etot0|=-9.025e-11`,
`dEhyd/Ehyd0=1.010e-4`, `n_ie_le0=0`, and zero floor, dual-energy, and
source-floor injection. This is live-prefix evidence only; its `t=2.221`
endpoint audit remains pending.

The q2 continuation root **437250** is still PENDING with `QOSGrpGRES`, and
its `t=4.5→6.5` continuation **437251**, later `t=6.5→20` chain, and exact-q2
OMP/MPI probes have not started. Thus there is still no `t>6` q2 result.

## Post-20 continuation staged (2026-10-27 23:22 KST)

To keep the long-run chain advancing toward the GIZMO comparison after its
current `t=20` endpoint, the guarded `t=20→22` W=0/fixed-mass segment was
staged and submitted as H200 job **437897**, with `afterok:437874`. It requests
one H200, two CPUs, 32 GiB, and 36 hours. The runner uses the existing
accepted binary/parameter checks, requires the exact `audit_t20.json`
predecessor receipt, then performs both the native endpoint and radial audits
at `t=22`. The staging helper verified the parent script hashes and passed
shell syntax, SHA-256, and Python import checks. The job is PENDING on its
dependency; it has not started and produces no new physics result yet. Its
submission is recorded in `chain_jobs_20261027.txt`.

At the same check, q1 W=0 job **437212** was RUNNING at `t=1.32939`, step
17,618, with `dEtot/|Etot0|=-9.449e-11`, `dEhyd/Ehyd0=1.051e-4`, no
nonpositive internal energy, and zero floor/dual-energy/source-floor
injection. q2 root **437250** remained PENDING with `QOSGrpGRES`; no q2
segment beyond `t=2.221` has begun, so the `t>6` gate remains untested.

## Continuation extended to t=24 (2026-10-27 23:25 KST)

The generic next-segment stager in
[`stage_kepler_polar_q2_halo28_next_20261027.py`](stage_kepler_polar_q2_halo28_next_20261027.py)
validated the staged t=20→22 runner as the exact predecessor, then prepared
and submitted **437909**, `t=22→24`, with `afterok:437897`. It requests one
H200, two CPUs, 32 GiB, and 36 hours; its runner requires `audit_t22.json`
and writes native plus radial audits at `t=24`. Both **437897** and **437909**
are PENDING on their `afterok` dependencies, so neither has started.

At the 23:24:59 KST check, q1 W=0 job **437212** remained RUNNING and had
advanced to `t=1.334`, step 17,682. Its latest RK4 ledger reports
`dEtot/|Etot0|=-9.898e-11`, `dEhyd/Ehyd0=1.063e-4`, `n_ie_le0=0`, and zero
floor, dual-energy, and source-floor injection. The q2 root **437250** still
waits with `QOSGrpGRES`; no q2 result beyond `t=2.221` or at `t>6` exists.

## Guarded continuation staged through t=30 (2026-10-27 23:27 KST)

Job **437912** (`t=24→26`, `afterok:437909`) had been accepted by Slurm but
was missing from the local chain ledger because the submission-ID parser
failed after `sbatch` returned success. The ledger is now repaired. Using
437912 as the predecessor dependency, the generic stager prepared
`t=26→28` and `t=28→30`; both passed parent hash checks, shell syntax,
SHA-256, and Python import validation. Slurm accepted **437915** (`t=26→28`,
`afterok:437912`) and **437916** (`t=28→30`, `afterok:437915`). The three
segments are recorded in `chain_jobs_20261027.txt`. They remain
dependency-gated; staging does not mean a segment has started or an endpoint
has passed.

At 23:27 KST, q1 W=0 job **437212** was still RUNNING on `syn104`; its latest
log time was `t=1.33768`, step 17,733. q2 root **437250** remained PENDING with
`QOSGrpGRES`; descendants through **437916** were pending on dependencies.
There is still no q2 result beyond `t=2.221`, so the required `t>6` stability
gate has not been tested.

## Live state and throughput checkpoint (2026-10-27 23:29 KST)

The next Slurm check confirmed q1 **437212** RUNNING and q2 root **437250**
still blocked by `QOSGrpGRES`. The queued continuation **437251**, **437868–
437874**, and **437897**, **437909**, **437912**, **437915**, **437916** remain
dependency-pending. q1 reached `t=1.34117`, step 17,780. Its latest RK4
ledger reports `dEtot/|Etot0|=-1.025e-10`, `dEhyd/Ehyd0=1.082e-4`,
`n_ie_le0=0`, and zero floor, dual-energy, and source-floor injection.

Between the 23:27 and 23:29 log reads, q1 advanced about `0.00349` simulation
time in roughly 127 seconds, or about `0.10` simulation-time units per wall
hour. A straight-line extrapolation from the current `t=1.341` to the GIZMO
Fig. 7 time `t=600` is roughly 250 days on this observed single-GPU rate.
This is only a scale estimate, not a forecast, but it identifies solver
throughput as a major remaining gate alongside stability. The timing record
also reports about 2.2 seconds per step and two GPU tessellation calls near
130 ms each; this suggests profiling the full step, since tessellation time
alone does not explain total wall time.

## q1 resource layout and revised throughput checkpoint (2026-10-27 23:31 KST)

The active q1 job's Slurm record confirms one H200, one MPI task, and two
CPUs; its staged runner invokes `mpirun -np 1` with `OMP_NUM_THREADS=2`.
Thus q1 is a single-rank control and provides no MPI scaling result. At the
23:30:47 log read it had reached `t=1.34328`, step 17,808, with the current
timestep reduced to `5.53e-5`; its RK4 ledger still showed
`dEtot/|Etot0|=-9.883e-11`, `dEhyd/Ehyd0=1.088e-4`, and zero floor or
dual-energy injections. Slurm confirmed it remained RUNNING at 23:31:55
(elapsed 5:11:53).

The latest measured interval advanced about `0.00211` in `t` over 89 seconds,
approximately `0.085` simulation-time units per wall hour. A linear extension
to `t=600` is about 290 days at that rate. The timestep is changing, so this
remains a rough scale estimate; it strengthens the case for obtaining the
matched multi-rank scaling measurements before deciding which code path to
optimize. The exact-q2 scaling probes are separate from this q1 control.

## Refined-resolution W=0 branch is advancing (2026-10-27 23:33 KST)

The separate A10 job **437771** is a 420² W=0 periodic-Voronoi run designed
to test the proposed inner-resolution increase: relative to the 384² parent,
it adds 536 angular sites across the 84 inner rings while preserving the
original coordinates at `r>=1`. Its exact-count `t=0.001` startup gate passed
earlier; that gate checked fixed masses, persistent W=0, periodic area tiling,
positive state, and orbital initialization. The live continuation to
`t=2.221` is now at `t=0.501194`, step 6,575. The latest RK4 ledger reports
`dEtot/|Etot0|=-8.480e-12`, `dEhyd/Ehyd0=2.176e-5`, and zero nonpositive
internal energy or floor injection. This is still a prefix, not the endpoint
audit.

The prepared A10 chain includes conditional timeout recovery at `t=2.221`
and a separately gated continuation to `t=6.5` with OMP=8. The latter is the
first planned `t>6` run for the refined mesh; neither it nor the H200 384² q2
branch has passed that gate yet. The refined run retains W=0 and fixed masses
and has no intrinsic cell spin, so it isolates resolution effects rather than
testing rotation or nonzero Laguerre weights.

At the 23:32:58 Slurm check, **437771** was RUNNING on one A10 with one MPI
task/two CPUs. H200 q1 **437212** was also RUNNING, while q2 root **437250**,
matched H200 OMP2 probe **437835**, and MPI2 probe **437892** remained blocked
by `QOSGrpGRES`; their dependent jobs were still held. Thus the A10 branch can
provide the next resolution evidence while H200 comparison and scaling jobs
wait for group GPU capacity.

## Follow-up live checkpoint (2026-10-27 23:35 KST)

The Slurm queue is unchanged: H200 q1 **437212** and refined A10 **437771**
are RUNNING; H200 q2 root **437250**, exact-q2 OMP2 **437835**, and MPI2
**437892** are still waiting on `QOSGrpGRES`. The matched four-GPU and OMP8
jobs remain dependency-held. No q2 t>6 result has been produced.

q1 advanced to `t=1.34962`, step 17,901, with `dt=7.64202e-5` restored.
From its `t=0.825` restart at 18:20:15 to the 23:35:09 log read, its mean
rate is about `0.100 t/hour`; a linear `t=600` extrapolation is about 250
days. This whole-segment average supersedes the noisier 290-day estimate from
the short timestep-dip interval. The A10 refined branch reached `t=0.504098`,
step 6,613; its latest ledger reports `dEtot/|Etot0|=-1.050e-11`,
`dEhyd/Ehyd0=2.202e-5`, and no floor or dual-energy injection. Its observed
mean since `t=0.001` is also about `0.101 t/hour`, but the distinct mesh and
GPU make that unsuitable as a hardware-scaling comparison. These are still
short prefixes and establish neither long-time physics nor MPI speedup.

## Current logs after queue check (2026-10-27 23:36 KST)

The 23:36:16 Slurm query still showed q1 **437212** and refined A10
**437771** RUNNING. The H200 q2 root **437250**, exact-q2 OMP2 **437835**,
and MPI2 **437892** remained at `QOSGrpGRES`; their descendants were held.
q1 reached `t=1.35165`, step 17,928. Its newest RK4 ledger has
`dEtot/|Etot0|=-1.018e-10`, `dEhyd/Ehyd0=1.111e-4`, no nonpositive internal
energy, and zero floor/dual-energy/source-floor injection. The latest `dt`
briefly dipped to `5.10e-5` after a `7.642e-5` step; this short event is not
evidence of sustained timestep collapse.

The 420² refined W=0 branch reached `t=0.506238`, step 6,641. Its latest RK4
ledger has `dEtot/|Etot0|=-1.256e-11`, `dEhyd/Ehyd0=2.221e-5`, and zero
nonpositive internal energy or floor injection. Both runs remain well before
their `t=2.221` endpoint audits; neither supplies a `t>6` result.

## Live checkpoint (2026-10-27 23:37 KST)

The queue remains unchanged: q1 **437212** and refined A10 **437771** are
RUNNING; q2 root **437250** and matched OMP2/MPI2 probes remain
`QOSGrpGRES`; all gated continuations remain pending. q1 reached `t=1.35306`,
step 17,947, with `dt=7.64202e-5`; its latest ledger reports
`dEtot/|Etot0|=-1.004e-10`, `dEhyd/Ehyd0=1.115e-4`, and no nonpositive
internal-energy cells or floor/dual-energy/source-floor injection. The
refined A10 run reached `t=0.50769`, step 6,660; its ledger reports
`dEtot/|Etot0|=-1.263e-11`, `dEhyd/Ehyd0=2.234e-5`, and zero floor injection.
These are short-time prefixes; the `t>6` stability and MPI scaling gates
remain open.

## Live checkpoint (2026-10-27 23:38 KST)

Slurm still reports q1 **437212** and the refined A10 run **437771** RUNNING.
H200 q2 root **437250** and the OMP2/MPI2 probes remain at `QOSGrpGRES`; all
dependent segments are still held. q1 reached `t=1.35512`, step 17,974,
with `dEtot/|Etot0|=-9.902e-11`, `dEhyd/Ehyd0=1.120e-4`,
`n_ie_le0=0`, and no floor injection. The A10 420² branch reached
`t=0.509753`, step 6,687, with `dEtot/|Etot0|=-1.117e-11`,
`dEhyd/Ehyd0=2.253e-5`, and no nonpositive-energy cells or floor injection.
Both logs have advanced, but neither has reached its `t=2.221` endpoint gate.

## Live checkpoint (2026-10-27 23:39 KST)

At 23:39:30, the same two jobs remained RUNNING and all H200 q2 and scaling
jobs retained their previous dependency/QoS waits. The q1 log advanced to
`t=1.35633`, step 17,995; its latest ledger has
`dEtot/|Etot0|=-9.860e-11`, `dEhyd/Ehyd0=1.124e-4`, zero nonpositive
internal-energy cells, and zero floor/dual-energy/source-floor injection.
The refined 420² A10 log advanced to `t=0.511358`, step 6,708; its ledger has
`dEtot/|Etot0|=-1.498e-11`, `dEhyd/Ehyd0=2.268e-5`, and zero floor injection.
No endpoint or `t>6` gate has passed.

## Live checkpoint and A10 chain record correction (2026-10-27 23:41 KST)

The queue remains unchanged. q1 **437212** reached `t=1.35916`, step 18,032;
its latest ledger reports `dEtot/|Etot0|=-9.768e-11`,
`dEhyd/Ehyd0=1.132e-4`, `n_ie_le0=0`, and zero floor/dual-energy/source-floor
injection. The 420² refined A10 job **437771** reached `t=0.514262`, step
6,746; its ledger reports `dEtot/|Etot0|=-1.216e-11`,
`dEhyd/Ehyd0=2.294e-5`, and zero nonpositive-energy cells or floor injection.
Both remain before the `t=2.221` gate.

Review of the A10 continuation script confirms job **437779** waits
`afterany:437772`, then refuses to evolve to `t=6.5` unless either the
original or conditional-resume run has passing exact-time global and radial
`t=2.221` receipts and a matching 420² checkpoint. It uses OMP=8 and audits
both global and radial behavior at `t=6.5`. The A10 README had retained the
superseded job number 437773; that record now correctly names 437779.

## Live checkpoint (2026-10-27 23:42 KST)

The queue remains unchanged: q1 **437212** and refined A10 **437771** are
RUNNING; H200 q2 root **437250** and the OMP2/MPI2 probes remain blocked by
`QOSGrpGRES`. q1 reached `t=1.36153`, step 18,063, with
`dEtot/|Etot0|=-1.003e-10`, `dEhyd/Ehyd0=1.138e-4`, `n_ie_le0=0`, and no
floor, dual-energy, or source-floor injection. The 420² A10 branch reached
`t=0.516701`, step 6,778; its latest ledger reports
`dEtot/|Etot0|=-1.317e-11`, `dEhyd/Ehyd0=2.316e-5`, and zero floor
injection. Neither branch has reached `t=2.221`, and no `t>6` result is
available.

Using the full observed segment averages through 23:43 KST, q1 has roughly
8.6 wall hours to its `t=2.221` endpoint, while the refined A10 run has about
16.9 hours. These are rate-based planning estimates only. The A10 estimate
leaves roughly two hours before its 24-hour limit; its conditional recovery
job is already staged if the original run times out before passing the gate.

## Live checkpoint (2026-10-27 23:43 KST)

At the latest Slurm check, q1 **437212** and refined A10 **437771** were still
RUNNING. H200 q2 root **437250** and OMP2/MPI2 probes remained
`QOSGrpGRES`; their continuations remained dependency-held. q1 advanced to
`t=1.36313`, step 18,084, with `dEtot/|Etot0|=-9.867e-11`,
`dEhyd/Ehyd0=1.142e-4`, zero nonpositive-energy cells, and no floor,
dual-energy, or source-floor injection. A10 advanced to `t=0.518305`, step
6,799, with `dEtot/|Etot0|=-1.176e-11`, `dEhyd/Ehyd0=2.331e-5`, and zero
floor injection. No endpoint receipt or `t>6` result exists yet.


## Live checkpoint (2026-10-27 23:44 KST)

The q1 H200 and refined A10 runs remained RUNNING; H200 q2 and the OMP2/MPI2
probes remained `QOSGrpGRES`, with their dependent chains held. q1 reached
`t=1.36517`, step 18,112, with `dEtot/|Etot0|=-1.008e-10`,
`dEhyd/Ehyd0=1.148e-4`, and zero nonpositive-energy cells or floor injection.
The 420² A10 branch reached `t=0.520598`, step 6,829, with
`dEtot/|Etot0|=-1.085e-11`, `dEhyd/Ehyd0=2.352e-5`, and zero nonpositive
internal energy or floor injection. Neither has passed the `t=2.221` gate.


## Live checkpoint (2026-10-27 23:45 KST)

The q1 H200 job **437212** and refined A10 job **437771** remain RUNNING.
H200 q2 root **437250** and the matched OMP2/MPI2 probes remain at
`QOSGrpGRES`; dependent continuations are still held. q1 reached `t=1.36663`,
step 18,131, with `dEtot/|Etot0|=-1.009e-10`, `dEhyd/Ehyd0=1.152e-4`,
`n_ie_le0=0`, and zero floor/dual-energy/source-floor injection. The 420² A10
branch reached `t=0.522126`, step 6,849, with
`dEtot/|Etot0|=-1.068e-11`, `dEhyd/Ehyd0=2.366e-5`, and zero nonpositive
internal-energy cells or floor injection. Both are still short of the
`t=2.221` endpoint gate; no `t>6` result exists.


## Live checkpoint (2026-10-27 23:46 KST)

At the latest queue check, q1 **437212** and refined A10 **437771** were
RUNNING. H200 q2 root **437250** and OMP2/MPI2 probes remained
`QOSGrpGRES`; dependent jobs remained held. q1 reached `t=1.36863`, step
18,159, with `dEtot/|Etot0|=-1.038e-10`, `dEhyd/Ehyd0=1.158e-4`,
`n_ie_le0=0`, and no floor/dual-energy/source-floor injection. The 420² A10
branch reached `t=0.524266`, step 6,877, with
`dEtot/|Etot0|=-1.173e-11`, `dEhyd/Ehyd0=2.386e-5`, and zero
nonpositive-energy cells or floor injection. Both remain before `t=2.221`;
no `t>6` endpoint has been reached.


## Live checkpoint (2026-10-27 23:47 KST)

q1 **437212** and the refined 420² A10 job **437771** remain RUNNING. H200
q2 root **437250** and the OMP2/MPI2 probes are still at `QOSGrpGRES`; their
dependent continuations remain held. q1 reached `t=1.36984`, step 18,176,
with `dEtot/|Etot0|=-1.035e-10`, `dEhyd/Ehyd0=1.161e-4`, no nonpositive
internal energy, and zero floor injection. A10 reached `t=0.525642`, step
6,895, with `dEtot/|Etot0|=-1.128e-11`, `dEhyd/Ehyd0=2.399e-5`, and zero
nonpositive-energy or floor-injection counts. Neither has passed the
`t=2.221` gate; there is no `t>6` result.


## Extension beyond t=30 staged (2026-10-27 23:49 KST)

The generic stager validated the t=30 parent and prepared the guarded
t=30→32 segment. Hash checks, shell syntax, and Python imports passed. Slurm
accepted **437932** with `afterok:437916`; a queue check confirms it is
PENDING on that dependency. It is recorded in `chain_jobs_20261027.txt`, so
the q2 chain now extends through t=32 if the preceding exact t=30 audits pass.

At the same check, q1 **437212** remained RUNNING at `t=1.37202`, step
18,205, with `dEtot/|Etot0|=-1.032e-10`, `dEhyd/Ehyd0=1.167e-4`, and zero
nonpositive-energy cells or floor injection. The refined 420² A10 job **437771**
remained RUNNING at `t=0.527934`, step 6,925, with
`dEtot/|Etot0|=-9.201e-12`, `dEhyd/Ehyd0=2.420e-5`, and zero floor
injection. H200 q2 root **437250** and the OMP2/MPI2 probes remain at
`QOSGrpGRES`; neither active branch has reached `t=2.221`, so no `t>6`
result is available.


## Guarded q2 continuation extended through t=40 (2026-10-27 23:49 KST)

Following the t=30→32 segment **437932**, four more two-unit W=0/fixed-mass
segments have been validated and submitted sequentially: **437933** (32→34,
afterok:437932), **437934** (34→36, afterok:437933), **437935** (36→38,
afterok:437934), and **437936** (38→40, afterok:437935). Each stager
verified the exact predecessor manifest and input hashes, passed shell syntax,
SHA-256, and Python import checks, and refused to copy restart data before the
predecessor physics gate. The chain ledger now contains all five additions
(437932–437936). They are queued behind their exact-time global and radial
audits; submission is not evidence that they have run or that any new endpoint
passed. The executed horizon remains unchanged and still requires the planned
t=600 comparison.


## Continuation queue and live-run check (2026-10-27 23:50 KST)

Slurm confirms jobs **437932–437936** are PENDING on their sequential
dependencies; the existing q2 chain through **437916** is also dependency-held.
The q2 root **437250** and exact-q2 OMP2/MPI2 probes remain `QOSGrpGRES`. q1
**437212** is RUNNING at `t=1.37401`, step 18,231; its latest ledger reports
`dEtot/|Etot0|=-1.058e-10`, `dEhyd/Ehyd0=1.173e-4`, and zero nonpositive
internal-energy or floor-injection counts. Refined A10 **437771** is RUNNING at
`t=0.529998`, step 6,952, with `dEtot/|Etot0|=-9.618e-12`,
`dEhyd/Ehyd0=2.440e-5`, and zero nonpositive-energy or floor-injection counts.
Neither active run has reached its `t=2.221` endpoint; the newly staged chain
has not begun executing.


## Live checkpoint (2026-10-27 23:51 KST)

The queue confirms q1 **437212** and refined A10 **437771** remain RUNNING.
H200 q2 root **437250** and OMP2/MPI2 probes are still `QOSGrpGRES`; all
continuations, including the new chain through **437936** (t=40), remain
held on dependencies. q1 reached `t=1.37551`, step 18,254, with
`dEtot/|Etot0|=-1.059e-10`, `dEhyd/Ehyd0=1.177e-4`, `n_ie_le0=0`, and
zero floor/dual-energy/source-floor injection. The 420² A10 branch reached
`t=0.531832`, step 6,976, with `dEtot/|Etot0|=-9.335e-12`,
`dEhyd/Ehyd0=2.457e-5`, and zero nonpositive-energy or floor-injection
counts. Neither has reached `t=2.221`; no `t>6` result is available.


## Live checkpoint (2026-10-27 23:54 KST)

Queried Slurm on `syntax` (the local `grammar` Slurm view is a different
cluster and must not be used for these jobs). q1 H200 **437212** and refined
420² A10 **437771** remain RUNNING. q1 reached `t=1.37951`, step 18,309,
with `dEtot/|Etot0|=-1.068e-10`, `dEhyd/Ehyd0=1.188e-4`, and zero
nonpositive-internal-energy cells or floor injection. A10 reached `t=0.536111`,
step 7,032; its latest energy audit at step 7,018 reports
`dEtot/|Etot0|=-1.028e-11`, `dEhyd/Ehyd0=2.498e-5`, and zero floor
injection. Neither has passed `t=2.221`.

H200 q2 root **437250** and exact OMP2/MPI2 probes **437835/437892** remain
`QOSGrpGRES`; all q2 continuations through t=40 remain dependency-held. The
A10 resume job **437772** is deliberately `afterany:437771` with a 48-hour
limit, and the t=6.5 A10 segment **437779** is `afterany:437772`, so the
planned chain can inspect/resume from the available checkpoint even if the
first 24-hour segment does not finish its t=2.221 audit. The q1 t=6.5 segment
**437213** remains `afterok:437212`, so it will only start after the q1 t=2.221
physics gate succeeds. No branch has yet established t>6 stability.

Recent profiler lines show approximately 0.9–1.0 s in CPU update/force and
about 0.12–0.14 s for an individual GPU tessellation call in q1; the A10
lines show roughly 0.8 s CPU update/force and 0.11–0.12 s tessellation.
This points to substantial non-tessellation cost in the current single-rank
runs and motivates reading the exact-q2 OMP/MPI probes once the H200 QOS allows
them; it does not yet prove which CPU-side kernel is the dominant optimization
target.


## Live checkpoint (2026-10-27 23:55 KST)

On `syntax`, q1 **437212** and refined A10 **437771** remain RUNNING;
q2 root **437250** and OMP2/MPI2 probes **437835/437892** remain held by
`QOSGrpGRES`. The q1 segment reached `t=1.38154`, step 18,336. A10 reached
`t=0.538515`, step 7,064; its latest energy audit (step 7,050) reports
`dEtot/|Etot0|=-1.226e-11`, `dEhyd/Ehyd0=2.521e-5`, and no nonpositive
internal energy or floor injection. The t=2.221 audits have not yet run.

The latest q1 log shows GPU tessellation calls around 0.236–0.242 s, above
the roughly 0.12–0.14 s seen at 23:54; this variation makes single-sample
GPU timing insufficient for attribution. The CPU update/force timings remain
around 0.8 s each in the latest visible step record. Preserve the planned
OMP/MPI comparison as the next useful scaling evidence once QOS releases the
probe jobs. The A10 run is progressing at about 0.10 simulation-time units
per wall hour averaged over this allocation; its `afterany` restart segment
**437772** remains queued as the recovery path if the 24-hour run reaches its
wall limit before finishing the t=2.221 audit.


## Live checkpoint (2026-10-27 23:56 KST)

q1 **437212** and A10 **437771** remain RUNNING. q1 reached `t=1.38299`,
step 18,355; A10 reached `t=0.540197`, step 7,086. The latest visible A10
energy record is from step 7,050 (`dEtot/|Etot0|=-1.226e-11`,
`dEhyd/Ehyd0=2.521e-5`, no floor injection). The t=2.221 gates are still
pending. q2 root **437250** and its OMP2/MPI2 probes remain `QOSGrpGRES`;
the t=40 continuation chain remains held behind the root and exact-time gates.

A scheduler query found 12 of this account's jobs RUNNING across A10, A40,
A100-PCIE, and H200 partitions at this check. The earlier accounting audit
in this document established the `part_h200` QoS cap as six GPUs. The current
H200 queue has exactly six one-GPU jobs RUNNING on `syn104`, including q1
**437212**; q2 root **437250** is therefore waiting on the partition QoS cap,
not on the node's physical eight-card capacity. The other two `kjhan` H200
allocations, **436963** and **437173**, belong to unrelated project paths and
were left untouched.


## Live checkpoint (2026-10-27 23:58 KST)

The H200 partition still has exactly six one-GPU jobs RUNNING on `syn104`,
matching the `part_h200` six-GPU group cap. q2 root **437250** and the two-rank
probe **437892** remain `QOSGrpGRES`; t=4.5, t=6.5 and later segments remain
dependency-held. q1 **437212** reached `t=1.38498`, step 18,381, and A10
**437771** reached `t=0.542566`, step 7,117. The latest A10 energy ledger at
step 7,103 reports `dEtot/|Etot0|=-1.410e-11`, `dEhyd/Ehyd0=2.560e-5`, and
zero nonpositive internal-energy cells or floor injection. q1 has not reached
its t=2.221 endpoint, and the active q2 chain has not advanced beyond the
already audited t=2.221 parent. No full-domain periodic t>6 result exists yet.


## Live checkpoint (2026-10-28 00:00 KST)

The H200 partition still has six RUNNING one-GPU jobs on `syn104`, exactly at
the `part_h200` group cap. q2 continuation root **437250** and the low-priority
two-rank probe **437892** remain `QOSGrpGRES`; the rest of the continuation is
dependency-held. q1 **437212** is RUNNING at `t=1.38774`, step 18,420. A10
**437771** remains RUNNING; its latest log is at step 7,149 with the most
recent tessellation timings near 0.12 s. Neither branch has reached its
current `t=2.221` audit. The periodic q2 chain's last accepted parent remains
the audited `t=2.221` state from job **437076**.


## Cross-history assessment against the long-term target

The repository already contains important `t>6` evidence, but it belongs to
distinct controls and must not be conflated with the current periodic q2 route:

- The bounded, aligned W=0 control passed its 17-check `t=6.5` gate and later
  reached `t=30`; however, its inner-disk native-minus-ballistic density
  residual grew 6.07x from `t=20` to `t=30`. Its static polygon boundary was
  eventually crossed at `t=48.3753`, so that branch cannot establish the
  periodic GIZMO `t=600` result. See
  [`KEPLER_BOUNDED_ALIGNED_T65_MILESTONE_2026-10-27.md`](KEPLER_BOUNDED_ALIGNED_T65_MILESTONE_2026-10-27.md),
  [`KEPLER_Q2_BOUNDED_W0_ALIGNED_T30_REVIEW_2026-10-27.md`](KEPLER_Q2_BOUNDED_W0_ALIGNED_T30_REVIEW_2026-10-27.md),
  and [`KEPLER_Q2_BOUNDED_WALLSELECT_LONG_CHAIN_2026-10-27.md`](KEPLER_Q2_BOUNDED_WALLSELECT_LONG_CHAIN_2026-10-27.md).
- Three bounded W=0 ICs pass at `t=6.5`, but independent angular site jitter
  has an 8.22x larger absolute density change and 91.6x larger native-minus-
  ballistic residual than the ring-coherent 5% IC. Within-ring spacing
  regularity is therefore a material IC sensitivity; aligned survival alone
  is not a general stability result. See
  [`KEPLER_Q2_THREE_IC_NATIVE_T65_REVIEW_2026-10-27.md`](KEPLER_Q2_THREE_IC_NATIVE_T65_REVIEW_2026-10-27.md).
- The full-domain periodic 384² W=0 route has passed its exact `t=2.221`
  parent audit. The active `2.221→4.5→6.5` H200 continuation and later
  segments are the outstanding route for resolving the historical `t≈6`
  failure under periodic boundaries. They are waiting on the known
  `part_h200` six-GPU group cap, currently occupied by six running jobs.
- A separate time-varying W64 run at `t=40` failed the independent Laguerre
  geometry audit: a boundary-fraction clamp created a positive stored volume
  for a hidden power cell and mismatched the face geometry used by fluxes.
  Therefore weight evolution is not ready to be used as a stability remedy;
  exact power-cell visibility and face-consistent geometry must be enforced
  first.

Cell spin is also absent from the live q1/q2 solver. The existing
`GFS_LAGUERRE_ROTATION` flag adjusts geometric face speed; it does not store or
evolve intrinsic cell angular momentum. The saved-cell affine analysis
estimates rigid spin as only about 10–11% of local linear Kepler subcell
kinetic energy, with most in non-rigid strain. Any native spin work therefore
needs an energy representation for that remainder and a consistent mass
centroid/first-moment and face-traction closure. It must preserve `Re=∞`
without adding friction. See
[`KEPLER_SPIN_NATIVE_INTEGRATION_GATES_2026-10-27.md`](KEPLER_SPIN_NATIVE_INTEGRATION_GATES_2026-10-27.md),
[`KEPLER_SPIN_CENTROID_CLOSURE_2026-10-27.md`](KEPLER_SPIN_CENTROID_CLOSURE_2026-10-27.md),
and [`KEPLER_SPIN_SHEAR_LOCAL_MODEL_2026-10-27.md`](KEPLER_SPIN_SHEAR_LOCAL_MODEL_2026-10-27.md).

Performance remains a gating requirement for `t=600`: the full-q2 profile
measured about 3.0 s per instrumented RK4 step, with four tessellation calls
near 0.315 s each and substantial CPU update/force work. A single-H200
extrapolation of roughly 250 days to `t=600` is only a rough estimate, but it
is too slow for an efficient verification campaign. Existing source audit
points to uninstrumented per-cell RK4 work and combined host staging/transfer
costs; the exact binary's archived source must be used before any code change,
and matched OMP/MPI timing evidence remains queued behind H200 QoS.

## Updated continuation status (2026-10-28 01:10 KST)

`sacct` confirms q2 parent **437076** COMPLETED after 16:30:41. Its exact
`t=2.221` global receipt passes all 15 checks and its radial receipt passes;
the scientific disk metrics and outer-background caveat are recorded in the
runtime profile. The `2.221→4.5` continuation **437250** is PENDING on
`QOSGrpGRES` (not on an unmet scientific dependency); **437251** and later
segments correctly remain dependency-held. The q2 chain has not yet crossed
`t=4.5` or `t=6.5`, so no full-domain periodic `t>6` claim is supported.
