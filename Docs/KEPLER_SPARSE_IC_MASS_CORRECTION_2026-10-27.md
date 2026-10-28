# Sparse q2 Kepler IC: initial-volume failure and correction (2026-10-27)

## Native failure identified before long evolution

The original 54,312-site narrow and 55,931-site wide preflights passed their
`t=.001` audits because those audits checked `m_i=rho_target(x_i)V_i(0)`,
positivity, fixed mass, disk orbits and disk pressure, but **did not check
that the initial periodic Voronoi areas tile the 6×6 box**. The assigned
initial masses were therefore internally consistent with incorrect areas.
The unmodified full q2 native initial snapshot is the comparison.

| Initial native IC | Sites | Sum of cell areas (required: 36) | Total fixed mass | Maximum cell area | Disk mass, 0.5≤r<2 |
|---|---:|---:|---:|---:|---:|
| Full q2 | 147,456 | 36.000000001 | 13.099331 | 0.0005749 | 11.904979 |
| Original sparse narrow | 54,312 | 30,469.710306 | 363.268457 | 44.444443 | 11.904979 |
| Original sparse wide | 55,931 | 30,761.608211 | 367.098226 | 44.444443 | 11.904979 |

At `t=.001` the original sparse native GPU tessellation gives total area
`≈36` (narrow `35.999999994`, wide `35.999999993`), but the already assigned
mass stays at 363/367. The narrow run then reaches density `33.415` even
though its intended initial profile is at most `1.01`. The initial total
kinetic energy of narrow is `71.6783`, versus full q2 `5.36253`; the
difference is mostly the incorrectly massive core and outer atmosphere.
The disk mass happens to agree, so annular checks alone concealed the error.

In the sparse-capable source, `Exam/KH/util.c::kh_findVol()` uses a local CPU
neighbor stencil to compute initial areas before the one-time Kepler mass
rebase. Sparse exterior sites leave some polygons nearly uncut, up to the
initial 6.667-unit square area `44.4444`. The GPU subsequently reconstructs
the actual periodic polygons, but fixed mass has already inherited the
incorrect areas. This is an **IC initialization defect**, not evidence that
the later native GPU tessellation has the same area error.

The four queued sparse continuations **437099–437102** were canceled before
execution (Slurm elapsed `00:00:00`). Their original preflight snapshots are
preserved for the regression audit. The [new auditor](audit_kepler_sparse_periodic_mass_preflight_20261027.py)
rejects the old narrow snapshot: initial area and mass differ from the
correct periodic IC by up to a factor `1876.78`, while its fixed-mass,
W=0, positivity, and disk-orbit checks still pass. The
[failed-audit receipt](KEPLER_OLD_SPARSE_MASS_GATE_FAILURE_2026-10-27.json)
preserves each check and the snapshot hashes.

The earlier periodic-seam taper used this same invalid narrow sparse glass.
Its preflight **436955** completed under the obsolete gates, and its
conditional continuation **436965** had started before the mass defect was
applied to this branch. Direct inspection of its `t=0` native snapshot gave
the same area sum `30,469.7103055` and mass `363.2684565`. Job 436965 was
canceled; its partial output cannot establish any boundary-physics result.
See the [taper correction](KEPLER_PERIODIC_SEAM_CONTROL_2026-10-26.md#superseding-sparse-mass-audit-2026-10-27).

## One-time periodic IC repair

The [area generator](make_kepler_periodic_area_sidecar_20261027.py) computes
W=0 Voronoi polygons for the exact original site positions with periodic
images at independent 1.5- and 2.5-unit margins. The maximum twice-cell
circumradius is `0.3281` (narrow) and `0.3105` (wide), safely below each
margin. Both independently computed area sums are 36 to about `1e-12`.
Their maximum cell-by-cell relative disagreements are `8.88e-11` and
`1.18e-10`. The immutable inputs are [narrow areas](KEPLER_Q2_PROFILE_BUFFERED_PERIODIC_AREAS_2026-10-27.json)
and [wide areas](KEPLER_Q2_PROFILE_WIDE_PERIODIC_AREAS_2026-10-27.json),
indexed by the original native glass site IDs. They imply total initial
masses `13.136045` (narrow) and `13.125032` (wide), with the same disk mass
`11.904979` as full q2. Different core/atmosphere sampling means their
whole-box masses need not equal full q2 exactly.

An isolated copy of the user's sparse-capable CUDA+MPI+OpenMP source at
`/tmp/Kepler_sparse_periodic_area_20261027` reads an opt-in
`EUNHA_KEPLER_INITIAL_AREAS_FILE` after the CPU initialization and uses its
site-ID area for the existing `m_i=rho_target(x_i)V_i(0)` rebase. It checks
the sidecar's count, box, exact length, positive finite entries, and total
area 36. It fails sparse Kepler startup if the sidecar or mass rebase is
missing. No RK4, force, EOS, GPU tessellation, or timestep formula changed.
The [source diff](KEPLER_SPARSE_PERIODIC_AREA_NATIVE_PATCH_2026-10-27.patch)
is against the previously accepted sparse glass source. The new binary
SHA-256 is `eff735174de1b487759f168d75e1e317a201692ef8d6cadba03cd03e49afb4b5`.
The complete changed and untracked native `Exam/` source overlay is preserved
in [KEPLER_NATIVE_SOURCE_OVERLAY_2026-10-27](KEPLER_NATIVE_SOURCE_OVERLAY_2026-10-27/README.md)
with per-file hashes and its base Git commit, independently of the temporary
`syntax` worktree.

Fresh H200 preflights **437156** (narrow) and **437155** (wide) use separate
run directories, byte-checked original glasses, their corresponding area
sidecars, the same new binary and physical settings, and the stronger audit.
The [narrow](kepler_q2_sparse_periodic_area_narrow_preflight_20261027.slurm)
and [wide](kepler_q2_sparse_periodic_area_wide_preflight_20261027.slurm)
runners differ only in job name, site count and input hashes. Both were
pending GPU quota at submission. Both have now **COMPLETED, exit 0:0** on
syntax H200 (7 seconds each), and their native `t=.001` audits passed all
14 checks. The audited binary SHA-256 is
`eff735174de1b487759f168d75e1e317a201692ef8d6cadba03cd03e49afb4b5`
in both runs. The authoritative receipts are
`/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_profile_buffered_area_t001/audit_t001.json`
and
`/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_profile_wide_area_t001/audit_t001.json`.

| Native `t=.001` gate | Narrow, 54,312 sites | Wide, 55,931 sites |
|---|---:|---:|
| Initial area sum | 36.0000001641 | 36.0000001535 |
| Evolved area sum | 36.0000000185 | 35.9999999947 |
| Initial fixed mass | 13.1360448619 | 13.1250321670 |
| Maximum relative initial area error against sidecar | 5.89e-8 | 5.70e-8 |

The gate also verified fixed mass across the step, W=0, positive state,
zero floor injection, disk pressure and generator circular speed. This
validates the corrected native initialization and first short evolution;
it does not establish sustained hydrostatic/orbital balance or a long-time
Kepler solution.

The [same-time native density-field image](KEPLER_FULL_SPARSE_NATIVE_DENSITY_T0001_2026-10-27.png)
compares full q2 and corrected sparse narrow at `t=.001` over the paper's
`[-2,2]²` view, using identical density and relative-change color scales.
It makes the sparse core and outer-profile sampling visibly coarser while
showing equal disk mass-weighted mean absolute density change to the displayed
precision (`2.04e-6`). The [plotter](plot_kepler_native_density_comparison_20261027.py)
requires complete same-time native snapshots, W=0, fixed mass, positive
state, area tiling, and `rho=m/V`; it can render later matched endpoints.
This early image is an IC check, not a long-time stability result.

The native orbital limiter `dt <= 0.06 (r²+0.01175²)^(3/4)` gives an initial
limit of `0.00193985` at the sparse first ring (`r=0.10080469`), versus
`0.0000764202` for the full-q2 central generator: about 25.4 times larger.
The exact-stop preflight deliberately truncated its only step to `.001`, so
this is a source/IC prediction, not a measured sustained-step speedup. The
`.825` runs will reveal the actual CFL and wall-time benefit.

The conditional native chain now continues only after each prior physical
audit exits successfully. Narrow **437160** and wide **437161** depend on
their respective corrected preflights and stop at exact `t=.825`; their
[generic early auditor](audit_kepler_sparse_area_t0825_20261027.py) checks
periodic area, fixed mass, W=0, energy, positivity, floor injection, disk
density, entropy, and orbital errors. Before either queued `.825` job began,
the auditor was strengthened to report the eight-site first ring separately
(`r<.2`) and the fixed-mass cell-centre crossings of `r=.5`. Promotion now
also requires first-ring mass-weighted absolute density change `<.05`,
absolute entropy change `<.01`, RMS radial speed/orbital speed `<.05`, and
total mass of inward plus outward generator crossings `<1e-4` of initial
disk mass. These are intentionally generous early-instability screens, not
proof of a balanced subcell force or long-time stability. The crossing
quantity is not a face-integrated Eulerian fluid flux. The new auditor SHA-256
`22691a957714845bf86bc63fd2acb85bb6280ad2a851718af7dba225746b6df4`
and updated `scripts.sha256` were checked in both held work directories;
jobs **437160/437161** were then released. On the two existing native
`t=.001` preflight snapshots, all 16 strengthened checks pass, including
first-ring density change `1.29019e-4`, absolute entropy change `5.6684e-8`,
and zero `r=.5` cell-centre crossing mass. This verifies the audit path only
at the preflight time.

Narrow **437174** and wide **437175**
depend on those respective `.825` jobs and stop at exact `t=10`. The
[checkpoint preparer](prepare_kepler_sparse_area_from0825_20261027.py)
requires every early gate, checkpoint hash, and the approved binary hash
before copying a restart. The [t=10 auditor](audit_kepler_sparse_area_t10_20261027.py)
adds whole-box periodic area and energy gates to the earlier annular GIZMO
metrics, and reports separate core, disk, guard, and outer radial bins.
Both [narrow](kepler_q2_sparse_area_narrow_to10_20261027.slurm) and
[wide](kepler_q2_sparse_area_wide_to10_20261027.slurm) runners request
H200:1, 32 GB, halo 32, and save intermediate snapshots through `t=10`.
Following an unrelated q1 final-auditor import packaging failure, the
actual Python import chains were checked successfully in all four staged
sparse `.825`/`t=10` work directories. Their `scripts.sha256` files were
also verified in place; no sparse job input or active binary was changed.
The `.001` preflights have passed. At the last 2026-10-27 check the `.825`
jobs **437160/437161** were pending `QOSGrpGRES`, and the `t=10` jobs
**437174/437175** were pending their `.825` dependencies. Neither stage has
a native result yet.

An additional [state-consistency auditor](audit_kepler_sparse_state_consistency_20261027.py)
checks complete IDs, fixed mass, W=0, periodic area, `rho=m/V`, and the
native ideal-gas relation `P=(gamma-1)U/V` at both endpoints. It complements
the time-specific orbit, energy, and entropy gates; the already staged GPU
jobs were not altered. On the existing narrow and wide native `t=.001`
preflights, both independent [narrow](KEPLER_SPARSE_STATE_CONSISTENCY_NARROW_T001_2026-10-27.json)
and [wide](KEPLER_SPARSE_STATE_CONSISTENCY_WIDE_T001_2026-10-27.json)
receipts pass all eight checks. Endpoint maximum relative `rho=m/V`
error is `5.90e-8` in both; EOS pressure errors are `9.16e-8` and
`9.17e-8`. Apply this same independent check to a completed `t=10`
checkpoint before promoting it toward `t=120` and `t=600`. These first-step
passes do not establish late physical accuracy.

### Prepared native continuation beyond t=10

The reusable [long-run H200 runner](kepler_sparse_area_long_20261027.slurm),
[audited restart preparer](prepare_kepler_sparse_long_restart_20261027.py),
and [bundle stager](stage_kepler_sparse_long_20261027.py) are prepared for
integer or decimal endpoints after `t=10`. The runner retains the accepted sparse CUDA
binary, `params.dat`, W=0 fixed mass, variational pressure, RK4, halo 32,
and no physical viscosity. It requests one H200, 32 GB and at most 48 hours,
with an explicit exact stop and configurable intermediate dump interval.
The prospective H200 runner now enables sparse host/device tessellation
timing at the first eight and every 1,024th call; its current SHA-256 is
`5e71f41ea3e0405348899e120974e97d4ea08156bb0e234ed74156a4cbab49e1`.
Each stage refuses a changed checkpoint, binary or physical input, and
requires both the predecessor's passed orbit/energy/entropy receipt and the
independent state-consistency receipt. The generalized
[native gate](audit_kepler_sparse_area_t10_20261027.py) accepts `--time`
for later endpoints and now also checks `rho=m/V` and the ideal-gas pressure
relation. Its default remains `t=10`; the currently queued `t=10` jobs use
their own immutable, previously staged auditor. The generalized gate passed
all checks on the actual corrected narrow native `t=.001` checkpoint;
shell syntax, Python compilation, staged imports and bundle SHA checks also
passed in an isolated temporary bundle.

No H200 job beyond `t=10` has been submitted: the required native `t=10`
checkpoint and its measured sustained timestep/wall cost do not yet exist.
After a `t=10` audit passes, run the independent state-consistency check on
its recorded endpoint, stage a new directory with the helper, and submit a
first interval to `t=120`. Later intervals can target `t=360` and `t=600`
if the previous physical gates and measured 48-hour wall budget permit;
intermediate checkpoints allow a shorter interval if the timestep shrinks.
For the separate longest-orbit objective, the current softened potential
gives `T(r=.5)=2.222361498867078`; a nominal 1000 inner-edge orbits ends at
`t=2222.361498867078`. The unsoftened convention gives
`t=2221.441469079183`. The runner now accepts either exact decimal target;
the paper's prose does not unambiguously specify the orbit reference radius
or its numerical softening, so the chosen convention must be stated in the
final comparison. A successful `t=600` image remains a distinct earlier gate.
This is a conditional calculation path, not a claim that the native solution
is stable at those times. The full-q2 and sparse `t>6` comparisons remain
the prerequisite for interpreting any sparse long-time result.

## Required promotion check

For each completed preflight, require initial and evolved area sums near 36,
initial areas matching the independently computed sidecar, fixed mass,
target-profile total mass, W=0, positive state, zero floor injection and disk
orbit accuracy. The submitted restart chain must refuse any corrected
checkpoint that does not pass these gates. The periodic outer velocity/potential seam and the sparse core's
point-source versus cell-average gravity discrepancy remain separate physical
questions for `t>6`, `t=120`, and `t=600`. A new
[initial core-force quadrature](KEPLER_SPARSE_CORE_GRAVITY_INTERPRETATION_2026-10-27.md)
measures the first-ring profile-weighted cell-mean gravity at about 1.308
times its generator value, versus about 0.993 in nearby full-q2 cells. It
confirms a core spatial-resolution limitation, without establishing that
this limitation causes disk heating or invalidates the native comparison.

## Whole-disk exact-orbit geometry control

Astra's IC review accepts full q2 and the corrected sparse layouts as valid
discrete Kepler test initial states, while reserving judgment on their native
long-time fluid stability. To distinguish mesh geometry from hydrodynamics,
the [whole-disk geometry auditor](audit_kepler_sparse_full_long_geometry_20261027.py)
uses each native `t=0` generator, area, and fixed mass, advances sites along
the same softened exact circular orbits, and rebuilds W=0 Voronoi cells at
`t=0, .825, 2.22144, 6.5, 10, 120, 600`. It compares density change across
the entire original `0.5 <= r < 2` disk, not only the first eight core cells.
Each disk polygon must be proven independent of omitted periodic images.

The independent [CPU Slurm runner](kepler_sparse_full_long_geometry_20261027.slurm)
completed on grammar as job **1223410**, exit `0:0`, in 61 s, with 1 CPU,
16 GB requested, and 470 MB peak batch RSS. The [result receipt](KEPLER_SPARSE_FULL_LONG_GEOMETRY_2026-10-27.json)
is copied from
`/gpfs/kjhan/LagForce/Kepler_Sparse_LongGeometry_20261027/result.json`.
The runner refuses to reuse an existing result directory and records SHA-256
for the three native initial snapshots, reader layout, and analysis scripts.

All three layouts contain **49,038 identical generator positions** in the
original `0.5 <= r < 2` disk. Their total original disk mass agrees to
`6e-11`; per-cell native initial areas differ by at most `5.96e-8` relative.
The following geometry-only disk mean absolute density changes are therefore
equal to the shown precision:

| Exact-orbit time | Full q2 | Corrected narrow | Corrected wide |
|---:|---:|---:|---:|
| 0.825 | 1.1654484e-4 | 1.1654484e-4 | 1.1654484e-4 |
| 6.5 | 1.8572658e-4 | 1.8572658e-4 | 1.8572658e-4 |
| 10 | 1.9429006e-4 | 1.9429005e-4 | 1.9429005e-4 |
| 120 | 2.0584731e-4 | 2.0584731e-4 | 2.0584731e-4 |
| 600 | 1.9698959e-4 | 1.9698959e-4 | 1.9698959e-4 |

Every disk polygon passed the omitted-periodic-image safety certificate;
its minimum edge-to-`2*circumradius` ratio exceeded 41. This result shows
that the sparse construction preserves the *disk's exact-orbit Voronoi
geometry*; removing sites elsewhere does not buy an intrinsic disk
geometric-error reduction. Any sparse/full difference in native evolution
must be interpreted through their different core/outer gas sampling, forces,
time stepping, or propagated boundary disturbances. This control does not
evolve gas or test entropy, pressure, forces, or the periodic seam. No
long-time stability claim follows from a favorable geometry result.

The [independent Eulerian inner-circle mass audit](KEPLER_EULERIAN_INNER_MASS_2026-10-27.md)
finds a further IC difference that the identical main-disk geometry does not
show: within `r<0.5`, sparse narrow/wide begin with 0.783%/0.768% more
piecewise-constant cell mass than full q2. Their short native runs create no
appreciable extra inner mass. Later comparisons must subtract each layout's
own initial enclosed mass rather than compare absolute sparse and full mass
as though the initial core sampling matched.

## Native same-time update (2026-10-27)

Both corrected sparse `.825` jobs **437160/437161** completed with all 16
early checks passing; their independent strict and state consistency audits
also passed. The [speed and same-disk comparison](KEPLER_SPARSE_NATIVE_T0825_SPEED_ACCURACY_2026-10-27.md)
records the measured 198-second sparse paths against the full-q2 H200 path,
the virtually unchanged initial-disk density field, and the much smaller
early entropy error. The earlier “pending `.825`” text above is historical.
The `t=10` native jobs **437174/437175** are queued behind H200 quota, and
their physical results remain necessary before any long-time promotion.

## Matched native disk and time-step control at t=0.825

The corrected narrow and wide H200 continuations **437160/437161** completed
at exact `t=0.825`, exit `0:0`, in 3:19 each. Both passed their staged early
gates. The independent [narrow strict](KEPLER_SPARSE_NARROW_T0825_STRICT_2026-10-27.json)
and [wide strict](KEPLER_SPARSE_WIDE_T0825_STRICT_2026-10-27.json)
audits pass 15/15 checks; the [narrow](KEPLER_SPARSE_NARROW_T0825_STATE_2026-10-27.json)
and [wide](KEPLER_SPARSE_WIDE_T0825_STATE_2026-10-27.json)
state audits pass 8/8. Neither run injected an energy floor. Whole-box
relative energy drift is below `3.7e-9` in both.

| Native layout | Disk mass-weighted mean absolute density change | Disk mean absolute log-entropy change | Signed log-entropy change | Disk orbit-position error |
|---|---:|---:|---:|---:|
| Full q2 | 1.16521e-4 | 1.29736e-4 | -1.89861e-5 | 2.50253e-9 |
| Sparse narrow | 1.16596e-4 | 2.50389e-6 | -6.30745e-7 | 4.12776e-9 |
| Sparse wide | 1.16521e-4 | 2.51848e-6 | -6.60358e-7 | 2.44020e-9 |

The [matched-site audit](KEPLER_FULL_SPARSE_MATCHED_DISK_T0825_2026-10-27.json)
matches all 49,038 original full-q2/narrow disk generators exactly at `t=0`.
Their masses differ by at most `6.06e-8` relative and initial velocities are
identical. At `t=0.825`, the mass-weighted mean absolute *pointwise difference*
in density change is `1.76e-7`; narrow's mean absolute entropy change is
about 1/52 of full q2's. The corresponding narrow/wide entropy magnitudes
differ by only 0.58%. Identical original disk sites and the near-equal
native disk geometry show that disk-site spacing alone cannot explain the
large full-versus-sparse entropy difference. The tested narrow/wide sampling
change has little early effect on this metric, but it does not isolate the
remaining changed variables: sparse core and outer gas sampling, their
propagated forces, and the global time step. Full q2 took 10,856 accepted
steps, versus 425 continuation steps for corrected narrow; a per-step
error hypothesis therefore needs a matched-IC time-step test.

The same narrow preflight and accepted binary are staged in the immutable
`/gpfs/kjhan/LagForce/Kepler_SparseNarrow_DTeta0006_T0825_20261027`
directory. The [control runner](kepler_sparse_narrow_dteta0006_to0825_20261027.slurm)
changes only `EUNHA_KEPLER_DTETA` from `0.06` to `0.006`, requests one H200,
and will audit the exact `t=0.825` endpoint with the early, strict, and state
gates. Slurm job **437240** was accepted with low scheduling priority and is
pending `QOSGrpGRES` at the 2026-10-27 04:47 KST check. This controls the
time-step contribution without assuming that the shorter early run proves
survival at `t>6`. The previously submitted narrow/wide `t=10` jobs
**437174/437175** remain pending the same GPU quota.

### Stored native entropy history without a new GPU run

The [CPU Slurm auditor](audit_kepler_q2_sparse_entropy_series_20261027.py)
and [runner](kepler_q2_sparse_entropy_series_20261027.slurm) analyzed four
complete saved checkpoint pairs. Grammar job **1223809** completed, exit
`0:0`, in 5 s with 569,432 KB peak batch RSS. Its [receipt](KEPLER_Q2_SPARSE_NATIVE_ENTROPY_SERIES_T0825_2026-10-27.json)
preserves all input hashes, exact native checkpoint times and IDs, and
periodic-area, fixed-mass, W=0, `rho=m/V`, and ideal-gas EOS checks.
Each full/sparse pair differs in physical time by less than `0.002`.

| Full q2 time | Full steps | Sparse steps | Full disk `mean |Δln K|` | Sparse disk `mean |Δln K|` | Full/sparse |
|---:|---:|---:|---:|---:|---:|
| 0.206265 | 2,700 | 107 | 3.74917e-5 | 1.34981e-6 | 27.78 |
| 0.412548 | 5,403 | 214 | 6.58035e-5 | 1.86804e-6 | 35.23 |
| 0.618764 | 8,119 | 320 | 9.66021e-5 | 2.11698e-6 | 45.63 |
| 0.825000 | 10,856 | 426 | 1.29736e-4 | 2.50389e-6 | 51.81 |

The mass-weighted disk density-change magnitudes agree across each pair to
within `7.5e-8` at the final endpoint; both are close to the same
exact-orbit tessellation prediction. Absolute entropy error grows with
accepted steps in both layouts, much faster in full q2. A straight-line fit
through only these four saved points gives `1.13e-8` per full step and
`3.49e-9` per sparse step; this is a descriptive trend, not a convergence
order or a per-step source measurement. Full q2's *signed* mean entropy
change is negative at every saved point, reaching `-1.90e-5` at `.825`;
the series does not demonstrate net irreversible heating. Core/exterior
sampling differs, and the checkpoint times are only approximately matched.
The same-IC `EUNHA_KEPLER_DTETA=0.006` control **437240** is the direct
next discriminator for accumulated-step error.

An ID-ordered direct read of the `t=0` and `.825` snapshots named and hashed
in the series receipt decomposes
`Δln K = Δln U + (gamma-1)Δln V` over the initial disk cells:

| Native layout | `mean |Δln U|` | `mean |(gamma-1)Δln V|` | `mean |Δln K|` | Unweighted correlation of the first two signed terms |
|---|---:|---:|---:|---:|
| Full q2 | 1.65798e-4 | 7.76812e-5 | 1.29736e-4 | -0.60249 |
| Sparse narrow | 7.77302e-5 | 7.77310e-5 | 2.50389e-6 | -0.99977 |

The two signed terms nearly cancel cell by cell in sparse and substantially
less well in full q2, despite their nearly identical disk density evolution.
This identifies the *thermal-work/volume cancellation* as the immediate
source of the measured entropy difference. It does not show whether the
underlying error comes from more RK4 steps, finite precision, different
core/exterior forces, or native stage-to-final geometry consistency. The
independent small-mesh adjoint identity remains valid; this snapshot
decomposition does not contradict it. The pending same-IC timestep control
and same-checkpoint energy-versus-fixed-K probe are needed before changing
the operator or activating a thermodynamic closure.

The complementary **full-q2 same-IC** control is staged at
`/gpfs/kjhan/LagForce/Kepler_Q2_DTeta030_T020626_20261027` and submitted
as H200 job **437246**, at scheduling nice 1000. Its
[runner](kepler_q2_dteta030_to020626_20261027.slurm) restarts the byte-checked
accepted full-q2 `t=.001` checkpoint with the same accepted binary,
fixed masses, W=0, pressure operator, no physical viscosity, and the same
8-cell early-time halo. It changes only `EUNHA_KEPLER_DTETA=0.06` to `0.3`
and stops at the *actual* saved baseline time
`t=0.20626460064652694`. The [restart preparer](prepare_kepler_polar_restart_20261026.py)
now records an optional exact target time while retaining its original
default `.825`; existing staged jobs remain immutable. The new run's
strict auditor checks complete native state, fixed mass, W=0, area tiling,
zero injections, disk density, entropy, orbital error, angular momentum,
total energy, and paper-annulus metrics. The larger step is a diagnostic,
not an accepted speedup: if it reduces entropy error but worsens orbital or
energy accuracy, the tradeoff must be reported. Along with the sparse
smaller-step job **437240**, it tests the step-count hypothesis in both
directions without changing either IC's site sampling. No new native
outcome exists at submission.
