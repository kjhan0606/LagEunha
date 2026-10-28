# Native polar Kepler control, prepared 2026-10-26

## Why this experiment exists

The native C010 Cartesian fixed-mass Voronoi disk had first-annulus density
drift almost entirely predicted by exact-orbit Voronoi geometry. The accepted
t=6.5 checkpoint still had 14.69% mass-weighted absolute drift at
`2<=r0<2.25`; the subsequent whole-domain RK4 failure at t=6.599 was a
separate inner-atmosphere cell. A polar generator layout directly tests
whether the grid's orientation under differential rotation is a major cause
of the disk drift. It does not by itself implement a general Laguerre fix.

## Matched paper-unit initial problem

This test uses the user's 384² CUDA+MPI+OpenMP executable and the Hopkins
profile on a 6×6 domain: inner disk radius 0.5, outer disk radius 2,
`Sigma=0.01+profile`, `P=1e-6`, gamma=5/3, and softened potential epsilon
0.01175. The background remains Cartesian. The generators with initial
radius below 2.35 (71,032 of 147,456) are reassigned to regular polar rings;
the total particle count and box stay fixed. Both q=1 and q=2 are independent
cases, where q is target radial spacing divided by azimuthal spacing and
nominal local area remains `h²`, `h=6/384`.

This matches the paper's disk profile, pressure, potential *form*, and nominal
disk spacing (`4/256 = 6/384 = 0.015625`). It is a geometry control rather
than an exact replay of the published MFM particle IC: the paper describes an
effective 256² distribution over `[-2,2]²` with particles selected to follow
the density profile, while this run retains all 384² cells over `[-3,3]²`
and explicitly imposes polar rings inside `r=2.35`. The paper does not give
the numerical softening length in its test description. Treat comparison to
its `t=600` density image as a physical benchmark with these IC differences,
not a bitwise or strictly matched-method replication. Primary source:
[Hopkins (2015), §4.2.4](https://arxiv.org/html/1409.7395).
The paper's resolution footnote also describes thinning particle sites in
the inner/outer profile regions and notes that particle methods can use
vacuum edges. It does not establish our populated 6×6 periodic torus as
the paper's boundary. The separate [boundary-gap audit](KEPLER_ISOLATED_BOUNDARY_GAP_2026-10-27.md)
records why the current native reflecting-wall paths and CUDA ghost-face
operator cannot simply be selected as a matched free boundary.

The [generator](make_kepler_hybrid_polar_glass_20261026.py) writes a native
`LAGMFM_GLASS_FILE`; nearest separation is 0.01496 for q=1 and 0.01058 for
q=2, with no collision or out-of-box position. The [isolated source
patch](KEPLER_POLAR_NATIVE_PATCH_2026-10-26.patch) adds two opt-in initial
condition operations: (1) after native initial tessellation, set each cell's
mass once to `rho_target(x_i) * V_i(0)`, keeping density and pressure exact on
the non-Cartesian layout; (2) save a genuine t=0 native snapshot. With W=0
the subsequent cell mass stays fixed. The patch does not change the frozen
C010 source or active binary; [build provenance](KEPLER_POLAR_NATIVE_BUILD_2026-10-26.json)
records source and executable hashes.

## Geometry-only prediction on the full 384² IC

Independent SciPy Voronoi areas were computed for exact softened Kepler
generator orbits with fixed initial cell masses. Mass-weighted absolute
density changes are below. They are predictions for an orbit-only control,
not results of the gas solver.

| Initial layout | t=0.825, r=0.5–0.75 | t=2.221, r=0.5–0.75 | t=2.221, r=0.5–2 | t=600, r=0.5–2 |
|---|---:|---:|---:|---:|
| Cartesian | 0.03137 | 0.39379 | 0.05004 | not evaluated |
| Polar q=1 | 0.000914 | 0.001131 | 0.000525 | 0.000605 |
| Polar q=2 | 0.000361 | 0.000436 | 0.000141 | 0.000197 |

The [full machine-readable control](KEPLER_HYBRID_POLAR_BALLISTIC_2026-10-26.json)
also records t=120 and annular p95/max errors. The polar error remains bounded
in this geometry-only test through the paper's t=600 image time. Native
pressure forces, central atmosphere motion, boundary transitions, and RK4
positivity may invalidate this prediction; those are the purpose of the
native tests.

The [density-field comparison at t=2.221](KEPLER_POLAR_DENSITY_GEOMETRY_T2221_2026-10-26.png)
shows the exact-orbit Cartesian and polar q=2c controls with identical
color scales. Its lower panels show relative density change; the figure
is explicitly a Voronoi geometry prediction, not a native gas snapshot.

An [independent long-time exact-orbit Voronoi audit](KEPLER_POLAR_Q1_Q2_LONG_GEOMETRY_2026-10-26.json)
uses the actual native q=1 and q=2 `t=0` generators and their fixed initial
cell masses. All evaluated disk polygons pass a periodic-image independence
check by a minimum factor greater than 41. The initial `0.5<=r<2` disk's
mass-weighted absolute density change, `|V(0)/V(t)-1|`, is:

| Time | q=1 | q=2 | q=2/q=1 |
|---:|---:|---:|---:|
| 6.5 | 5.59457e-4 | 1.85727e-4 | 0.332 |
| 17.7789 | 5.79088e-4 | 2.16460e-4 | 0.374 |
| 120 | 6.46423e-4 | 2.05847e-4 | 0.318 |
| 600 | 6.04586e-4 | 1.96990e-4 | 0.326 |
| 2221.4415 (nominal 1000 inner-edge periods) | 5.95912e-4 | 2.12966e-4 | 0.357 |

At t=0, both errors are about `1.5e-8`, consistent with float32 stored
volumes. Thus q=2's geometric advantage persists through the paper's
plotted `t=600` interval and the longer nominal orbit count, without
secular growth in this *geometry-only* metric. It does not establish
native gas stability, heating, or orbital accuracy at those times.

## Same-shell versus adjacent-shell spacing

For a ring of radius `r` with `nphi` generators, the same-shell nearest
distance is `2*r*sin(pi/nphi)`. The distance between generators on adjacent
shells is at least their radial separation `dr`, regardless of their changing
relative phase. The actual generated 384² layouts give, near `r=1`:

| Layout | Same-shell chord | Adjacent-shell lower bound | `dr/chord` |
|---|---:|---:|---:|
| q=1 | 0.01543765 | 0.01562500 | 1.012 |
| q=2 | 0.01088462 | 0.02209709 | 2.030 |

Across every ring within `0.5<=r<=2`, the radial lower bound is greater
than the same-shell chord in both layouts. Thus initially closer
adjacent-shell generator pairs cannot explain their disk errors. However,
different shells still share Voronoi faces and rotate at different angular
speeds; face shape, neighbor identity, and numerical flux may change under
shear. The q=2 exact-orbit geometry reduces disk density error at `t=2.221`
by 3.7-fold relative to q=1, but the native `t>6` heating comparison has not
been made. Changing `q` also changes the angular-velocity jump across each
radial interface, so these two runs alone cannot isolate a flux mechanism.

The next discriminating native diagnostic is a face-level audit of neighbor
changes, radial versus tangential face orientation, and irreversible internal
energy production, grouped by radial separation. This tests whether the
remaining error tracks shear across shell boundaries rather than initial
generator proximity.

An independent [face auditor](audit_kepler_polar_faces_20261026.py) now reads
the complete native `t=0` and `t=0.001` snapshots. Among faces whose two
cells begin inside the disk, q=1 has 96,251 adjacent-shell faces out of
144,709, and q=2 has 96,593 out of 145,631. Thus both layouts have about
two-thirds of their disk face *count* between shells, despite those generator
pairs being farther apart than same-shell neighbors. The adjacent-shell
*total face length* is 799.83 in q=1 and 535.95 in q=2, a 33.0% reduction;
as fractions of all disk-disk face length these are 55.8% and 34.0%.
No cross-shell face has a generator separation shorter than either endpoint's
same-shell chord. At `t=0.001`, the native and exact-orbit Voronoi disk edge
graphs match exactly in both layouts; the entropy-change absolute mean is
only `3.81e-7` (q=1) and `2.76e-7` (q=2). These observations identify
cross-shell *contact length*, rather than nearest-center ordering, as a
plausible geometric control. They do not establish the cause of late native
heating: q=2 also increases the angular-velocity difference across each
adjacent-shell interface. [q=1 receipt](KEPLER_POLAR_Q1_FACES_T001_2026-10-26.json),
[q=2 receipt](KEPLER_POLAR_Q2_FACES_T001_2026-10-26.json).

A [native-velocity face diagnostic](KEPLER_POLAR_CROSS_SHELL_SHEAR_T0_2026-10-26.json)
adds the local co-rotating shear speed `r_mid |Omega_i-Omega_j|` and the
approaching component of the generator-pair velocity on every initial disk
face. It confirms that radial spacing alone cannot predict late numerical
heating. Relative to q=1, q=2 has **0.670 times** the cross-shell contact
length but **1.414 times** the contact-length-weighted local shear speed.
The integrated length-times-shear-speed is **0.947 times** q=1, whereas the
length-times-shear-speed-squared integral is **1.337 times** q=1. The
integrated positive approaching speed across these faces is **0.497 times**
q=1. Same-shell shear and approach are zero to roundoff initially. These are
kinematic proxies from the [independent script](audit_kepler_cross_shell_shear_20261026.py),
not an AV work integral or measured entropy production. A late native
q=1/q=2 comparison remains necessary to identify which proxy, if any,
predicts the actual gas error.

An independent [face-motion gap audit](KEPLER_FACE_MOTION_GAP_INTERPRETATION_2026-10-27.md)
compares analytic Kepler material motion with exact W=0 Voronoi face motion
on the same native t=0 ICs. The integrated absolute gap is larger on
**same-shell** than cross-shell disk faces in both layouts; q=2 cuts the
mass-weighted absolute geometric cell-area rate by about half. This is a
kinematic material-volume diagnostic, not the native HLLC contact speed or
measured heat production. The active variational W0 CUDA work law uses the
geometric face speed, so this gap must not be called a thermal-work defect.
It makes a cross-shell-only geometric explanation insufficient; native
entropy and face-force data at t>6 are still needed to attribute heating.

## Native q=1 first continuation checkpoint (2026-10-26)

The first complete continuation snapshot is native RK4 step 2700 at
`t=0.20626460064652694`, job 436889. The
[checkpoint audit](KEPLER_POLAR_Q1_T020625_AUDIT_2026-10-26.json) passes finite
state, exact fixed mass, W=0, positive volume, and positive internal energy.
For initial `0.5<=r<2`, mass-weighted absolute native density change is
`4.9132123e-4`, versus `4.9133441e-4` from the independent exact-orbit
Voronoi control. Their per-cell density difference has a mass-weighted
absolute mean of `2.634e-8`. The first annulus `0.5<=r<0.75` similarly gives
`0.0012121699` native versus `0.0012122019` geometry. The
[native density image](KEPLER_POLAR_Q1_NATIVE_DENSITY_T020626_2026-10-26.png)
shows the same ring/spoke error pattern in the native field and orbit-only
prediction.

The [face audit](KEPLER_POLAR_Q1_FACES_T020626_2026-10-26.json) finds 11,328
initial disk-disk faces gone and 11,418 new after differential rotation, but
the full native and exact-orbit disk face graphs at this time are identical.
No cross-shell generator pair on a disk face is nearer than either endpoint's
same-shell nearest neighbor. The disk's mass-weighted signed logarithmic
entropy change is `-3.09e-7`; its absolute mean is `5.74e-5`. The latter is
an error/noise scale, **not evidence of net irreversible heating**. A longer
native evolution and per-face flux or entropy accounting are still required.

A [per-cell face-turnover audit](KEPLER_Q1_FACE_FLIP_ENTROPY_T020626_2026-10-26.json)
finds 17,759 of 48,680 disk cells touched by a lost or new disk-disk face
at this checkpoint. Across the disk, changed-face count has Spearman
association `+0.195` with absolute density change but `-0.104` with
absolute entropy-log error. Within four radial bands, the entropy
associations range only from `-0.054` to `+0.015`. Thus early face turnover
is not a positive predictor of native entropy error in this run, although
it accompanies some of the geometric density pattern. This is a cell-level
association, not a per-face flux-work calculation or a result for `t>6`.

## Native q=1 second continuation checkpoint (2026-10-26)

The complete native step-5403 snapshot is at `t=0.41254840604704807`.
The [checkpoint audit](KEPLER_POLAR_Q1_T041254_AUDIT_2026-10-26.json)
passes finite state, exact fixed mass, W=0, positive volume, and positive
internal energy. For the initial `0.5<=r<2` disk, mass-weighted absolute
density change is `5.7094817e-4`, versus `5.7101075e-4` on exact circular
orbits. Native density differs from the independent orbit-only Voronoi
prediction by only `1.07086e-7` mass-weighted absolute mean. The native
disk orbit displacement is `2.226e-9` mass-weighted, and its angular
momentum drift is `-4.19e-13`. The [density field](KEPLER_POLAR_Q1_NATIVE_DENSITY_T041254_2026-10-26.png)
shows the ring/spoke residual without visible disk disruption. In the
published `0.8<r<1.2` annulus, [L1 velocity](KEPLER_POLAR_Q1_GIZMO_METRICS_T041254_2026-10-26.json)
is `3.031e-9`, literal L1 density `0.01000058`, and the initial plateau
offset-corrected L1 density `7.302e-4`.

The [face audit](KEPLER_POLAR_Q1_FACES_T041254_2026-10-26.json) finds
22,706 initial disk-disk faces gone and 22,801 new, but native and exact
orbit face graphs still agree exactly. No cross-shell disk face has
generator centers closer than a same-shell nearest neighbor. Mass-weighted
signed disk entropy-log change is `-1.90e-6` and its absolute mean is
`9.585e-5`, which does not demonstrate net irreversible heating. The
[per-cell turnover audit](KEPLER_Q1_FACE_FLIP_ENTROPY_T041254_2026-10-26.json)
finds 29,062 disk cells touched by a face change, with Spearman correlation
`-0.175` between change count and absolute entropy-log error. The four
radial-band correlations range from `-0.137` to `+0.001`; face turnover
does not positively predict this early entropy error.

The [periodic pair-distance audit](KEPLER_Q1_PERIODIC_PAIRS_T041254_2026-10-26.json)
finds exactly 44 native generator pairs below the pair-pressure distance
cutoff, involving 88 distinct cells. Their **current** radii are
`2.9819`–`3.3389`, entirely outside the disk. The nearby `[RK4E]` log
reports `n_pair=88` face ends and zero cumulative floor/dual-energy
injection. This match places the observed short-pair **criterion** in the
outer background at this checkpoint; it does not establish that the same
will hold later. The current disk has no pair below the cutoff. The q=1
run continues toward `t=.825`; neither this checkpoint nor the pair
location proves `t>6` stability.

A [dedicated periodic-seam receipt](KEPLER_Q1_SEAM_T041254_2026-10-26.json)
locates all 88 short-pair endpoint cells only `0.0458`–`0.1146` from a
box edge. The 44 pairs do not themselves straddle the periodic seam at
this snapshot, but are separated by only `0.00136`–`0.00227` while their
relative speeds are `0.944`–`1.142` (median `0.989`). At the same time,
cells within `0.2` of an edge have mass-weighted absolute density and
pressure changes `0.0681` and `0.1298`, versus `0.000571` and `0.000944`
in the disk. The seam area is already strongly perturbed, although its
signed entropy-log change is still tiny. This supports tracking the
propagation of outer disturbances in any t=120/600 comparison; location
and correlation alone do not prove the periodic seam is the only cause.

During the continuing q=1 run, the native `[RK4E]` log first reports
`n_pair=8` at `t=0.287729`, step 3753. This counts face ends satisfying the
short-pair criterion, not all geometric neighbors. In the active variational
W0 branch, the computed pair pressure does not contribute to the force or
work; see the [operator audit](KEPLER_POLAR_OPERATOR_IDENTITY_2026-10-27.md).
At an inspected complete log line near `t=0.337606`, `n_pair=96`,
while `floor_cum=de_cum=sfl_cum=0`. The log does not identify the radii of
those face ends, so this is an early-crowding signal, not evidence that the
disk itself is heating. The next complete checkpoint should be used to
locate the close pairs and measure their entropy change.

The [periodic nearest-generator screen](KEPLER_POLAR_PERIODIC_PAIR_SCREEN_2026-10-26.json)
checks the same distance cutoff with minimum-image distances. At the complete
native `t=0.2062646` checkpoint it finds no cell whose nearest generator is
inside the pair-pressure cutoff; the minimum separation/cutoff ratio is
`2.291`. A geometry-only exact-orbit screen using **initial** cell volumes
finds no such candidate at `t=0.287729`, but at `t=0.337606` finds 18
short generator pairs involving 36 cells, all starting at `r>=2.35`
(24 cells inside `r<3`, 12 outside). This supports an
outer-background origin for the later pair encounters, but does not locate
the *native* active pair ends. The screen tests all generator pairs within
the maximum cutoff, but not whether each pair shares a Voronoi face or
passes the pressure-work limits. Periodic wrapping is
essential here: an unwrapped outer orbit creates spurious near-coincident
points. Exact circular orbits for initial `r>=3` are not physical controls
after crossing the periodic boundary. The
[script](audit_kepler_polar_periodic_pairs_20261026.py) records these limits.

At this same time, a [matched-time geometry-only comparison](KEPLER_POLAR_Q1_Q2_BALLISTIC_T020626_2026-10-26.json)
predicts disk density error `1.3681e-4` for q=2 versus `4.9133e-4` for q=1,
or 3.59-fold smaller. This is not a native q=2 result; its queued H200
continuation must validate the prediction.

The [paper-metric reader](audit_kepler_gizmo_metrics_20261026.py) also
reproduces Hopkins (2015) Fig. 8 definitions in the current-radius annulus
`0.8<R<1.2`: mean `|v_phi-v_K|/v_K` and mean `|Sigma-1|`.
At this q=1 native checkpoint they are `9.60e-10` and `0.01000059`.
The paper's literal density metric already starts at `0.01` here because
its plateau profile is `Sigma=1.01`; the baseline-corrected mean
`|Sigma-1.01|/1.01` is `7.2465e-4`. This is a reproducible early benchmark,
not evidence for the paper's t=600 result. See the
[metric receipt](KEPLER_POLAR_Q1_GIZMO_METRICS_T020626_2026-10-26.json).

## Native gates and next action

Independent H200 preflights q=1 job 436885 and q=2 job 436886 use the same
patched native binary, exact t=0 and t=0.001 output, no generator steering,
W=0, Re=0, RK4, and one H200 each. At submission both were pending
`QOSGrpGRES`. Each [audit](audit_kepler_polar_preflight_20261026.py) checks
the t=0 geometry and target mass-density relation, exact fixed mass across
the step, positive volume and internal energy, W=0, and orbital motion.
Promote only passing cases to longer native checkpoints at equal orbital
counts. The current full-RK4 throughput must be measured before committing
to the paper's t=600 and ~1000-orbit objectives.

## 2026-10-26 21:10 KST: native preflight and centerless variant

Both q=1 and q=2 native t=0.001 preflights passed every check: exact fixed
mass, W=0, positive V/U, the prescribed t=0 mass-density relation, and
disk ballistic position error below `3.1e-14`. Each used 14 RK4 steps,
with normal initial `dt=7.64202e-5` and a shortened last step. q=1
continuation 436889 had begun and was at t=0.01155 after about 2m47s;
q=2 continuation 436890 was queued for GPU quota. These are interim
observations, not final t=0.825 audits.

The native orbital limiter is `dt <= eta*(r²+eps²)^(3/4)`. The central
generator at r=0 sets `dt=7.64202e-5` in both original polar layouts.
The Cartesian grid's innermost generator is at r=0.0110485, giving
`dt=1.22899e-4`. A third layout, **q=2c**, omits the central generator
and starts with six sites on a ring at the latter radius. It retains the
same 147,456 total cells and nominal q=2 area. Its predicted initial
limiter recovers `1.22899e-4`, a 1.61-fold larger step; actual long-run
speed has not yet been measured.

The [q=2c exact-orbit control](KEPLER_Q2C_BALLISTIC_2026-10-26.json)
predicts disk density change 0.000200 at t=2.221 and 0.000243 at t=600,
versus q=2 center layout 0.000141 and 0.000197. The
[central atmosphere control](KEPLER_POLAR_CORE_BALLISTIC_2026-10-26.json)
stays below 0.001 mass-weighted absolute density change at t=600 for both
q=2 layouts, though q=2c has larger individual extrema. Native q=2c
preflight job 436893 and its afterok t=0.825 continuation 436894 test
whether the speed gain survives full RK4/CUDA.

## Time-limit preservation

Each t=0.825 job has an `afternotok` fallback in a new directory: 436895
after q=1 job 436889, 436896 after q=2 job 436890, and 436897 after q=2c
job 436894. The fallback script reads the actual parent Slurm state and
refuses physical failure or user cancellation. Only TIMEOUT, PREEMPTED,
NODE_FAIL, or BOOT_FAIL may continue. It checks native snapshot completeness,
IDs, exact fixed masses, W=0, finite fields, positive volume, and positive
internal energy before copying the latest accepted checkpoint. It never
uses `khdiag.*.failed.dat` as a restart state and never overwrites its
parent output. If the parent completed normally, the `afternotok` job does
not run. The [selection code](prepare_kepler_polar_timeout_resume_20261026.py)
and [batch script](kepler_polar_timeout_resume_20261026.slurm) are saved
with the run.

The [time-step and throughput audit](KEPLER_TIME_STEP_SCALING_2026-10-26.md)
quantifies how one central cell limits the entire 384² RK4 update and records
an independent queued q1 GPU local-periodic-sync A/B job. Its active-cell
work ratios are not measured wall-time accelerations.

The [cell-spin status and conservation gate](KEPLER_CELL_SPIN_STATUS_2026-10-26.md)
separates the implemented W=0 face-velocity correction from the still-absent
intrinsic cell spin, and specifies the angular/energy audit needed if a spin
branch is introduced after the polar physical comparison.

## Exact-orbit geometry comparison at native t=0.001

The new [arbitrary-checkpoint auditor](audit_kepler_polar_checkpoint_20261026.py)
compares the actual native density with independent SciPy Voronoi areas after
placing the same generators on exact softened circular orbits. It was run on
the completed q1 and q2 native t=0.001 preflights. All state gates passed
(same IDs and fixed masses, W=0, finite fields, positive volumes and internal
energies). For the initial disk 0.5≤r<2:

| Layout | Native mass-weighted absolute density change | Exact-orbit geometry prediction | Native versus geometry residual |
|---|---:|---:|---:|
| q1 | 4.13267e-6 | 4.13304e-6 | 1.50683e-8 |
| q2 | 2.03556e-6 | 2.03625e-6 | 1.51133e-8 |

This checks a very early interval only. At t=0.001, disk generator positions
are within ~5.5e-15 mass-weighted of their exact orbits, so the measured
density signal is essentially the fixed-mass Voronoi geometry signal. The
spoke/ring structure visible in the error panels is present in both native
and orbit-only geometry, consistent with ring-count transitions rather than
an early hydro instability. See the [q1 receipt](KEPLER_POLAR_Q1_T001_GEOMETRY_AUDIT_2026-10-26.json),
[q2 receipt](KEPLER_POLAR_Q2_T001_GEOMETRY_AUDIT_2026-10-26.json),
[q1 native density figure](KEPLER_POLAR_Q1_NATIVE_DENSITY_T001_2026-10-26.png),
and [q2 native density figure](KEPLER_POLAR_Q2_NATIVE_DENSITY_T001_2026-10-26.png).
The same auditor and [density-field renderer](plot_kepler_polar_native_checkpoint_20261026.py)
are ready for the first full checkpoint at t=0.20625 and later snapshots.

## Native q=2 first continuation checkpoint (2026-10-27)

The complete q2 native step-2700 snapshot is at `t=0.20626460064652694`
(job 436890). In the initial disk `0.5<=r<2`, mass-weighted absolute density
change is `1.368102377e-4`; independent exact-orbit unweighted Voronoi
geometry predicts `1.368122657e-4`. The native-versus-ballistic absolute
density residual is `1.617825283e-8`. At the same time q1's native disk
density change is `4.9132123e-4`, about 3.59 times larger. This comparison
measures initial-layout geometry and very early native evolution; it does not
measure long-time shear heating.

All snapshot states are finite, IDs retain their exact initial masses,
Laguerre weights remain zero, and volumes and internal energies are positive.
The q2 disk's mass-weighted absolute log-entropy change is `3.74917e-5`;
its mass-weighted position error against exact softened circular orbits is
`1.22608e-10`. The native disk face graph is identical to the exact-orbit
face graph at this checkpoint, even though 22,894 initial disk faces have
disappeared and 22,934 new ones have appeared. Face turnover alone therefore
does not establish spurious heating. The outer `r>=2.75` region already has
mass-weighted absolute radial velocity about `0.0181`, so the periodic seam
and ambient remain separate long-time risks.

See the [checkpoint audit](KEPLER_POLAR_Q2_T020625_AUDIT_2026-10-27.json),
[radial audit](KEPLER_Q2_RADIAL_T020625_2026-10-27.json),
[face audit](KEPLER_POLAR_Q2_FACES_T020625_2026-10-27.json), and
[native density field](KEPLER_POLAR_Q2_NATIVE_DENSITY_T020625_2026-10-27.png).

## IC admissibility and benchmark interpretation (Astra review, 2026-10-27)

The q=2 polar ring layout is a mathematically valid finite-anisotropy
sampling of the same continuum Kepler problem: it starts with the prescribed
density, pressure, potential, and circular speed; every site thereafter
moves under the native equations, without a ring constraint or steering.
Fixed initial masses and W=0 remain valid. Its radial-to-azimuthal spacing
ratio is about 2.03 in the disk. This tests whether a shear-aligned initial
Voronoi geometry reduces error, but it must be named as an aligned-IC result.
Hopkins' published effective-resolution construction begins with evenly
spaced particles and removes low-density-region particles; it does not
prescribe this q=2 polar layout. The paper explicitly distinguishes
geometry-tuned disk meshes from general solver behavior
([Hopkins 2015, §4.2.4](https://arxiv.org/html/1409.7395#S4.SS2.SSS4)).

The hypothesis that *closer adjacent shells alone* cause spurious heating is
not established: q=1 already has cross-shell spacing slightly greater than
its same-shell chord, and even q=0.5 polar geometry beats the Cartesian
control in the earlier geometry-only experiment. The q=2 contact length is
0.67 times q=1, but its contact-length-weighted shear-speed-squared integral
is 1.337 times q=1. Neither proxy measures irreversible heating. Preserve
q1, q2, and the Cartesian control; compare native entropy, face work, and
torque at equal physical times through at least t=10. A small ring-phase or
position-perturbation q=2 control will test dependence on exact alignment.
Only a same-IC solver comparison can support a claim of solver improvement.
The subsequent [q2 phase-jitter geometry control](KEPLER_Q2_PHASE_JITTER_INTERPRETATION_2026-10-27.md)
shows that 5% independent site offsets raise its exact-orbit t=600 mean
density drift from `1.97e-4` to `1.47e-3`, whereas coherent ring offsets of
the same amplitude raise it only to `2.11e-4`. Along-ring regularity is part
of this IC's geometric advantage; no native heating conclusion follows.

### Astra follow-up after the native t=.825 comparison (2026-10-27)

Astra reaffirmed that q2 is an admissible finite-resolution sampling of the
constant-pressure softened circular equilibrium, with fixed initial masses
and no constraint on subsequent generator motion. It is a **shear-aligned IC
control**, not a replay of the paper's particle IC or a general solver
validation. The native [q1/q2 comparison](KEPLER_Q1_Q2_NATIVE_T0825_COMPARISON_2026-10-27.md)
now extends the observed result to `t=.825`: disk mass-weighted absolute
density change is `1.16521e-4` for q2 versus `5.23520e-4` for q1. The
independent exact-orbit q2 Voronoi prediction is `1.16545e-4`, with a
pointwise native-versus-ballistic disk residual of `8.76967e-8`.

This strengthens the **early geometry** interpretation. It does not show
that q2 prevents spurious shear heating at `t>6`. Neither q1 nor q2 has a
cross-shell generator nearer than its same-shell neighbor, and the face
motion audit finds a larger integrated material-versus-Voronoi motion gap
on same-shell disk faces. Distance ordering alone cannot establish a heat
source. The smallest discriminating native control remains unperturbed q2
versus 5% independent site-jitter q2 under the same solver, boundary, and
timestep, with entropy, torque, orbital error, and density relative to each
IC's exact-orbit geometry checked through `t>6`. The staged jitter jobs
currently stop at `t=.001` and `.20626460064652694`; they cannot settle
this late-time question by themselves.

## Native q=1 third continuation checkpoint (2026-10-27)

The complete q1 step-8119 checkpoint is at `t=0.6187644470884642`
(job 436889). All 147,456 states are finite; initial masses and W=0 are
exactly preserved; volumes and internal energies remain positive. In the
initial `0.5<=r<2` disk, mass-weighted absolute density change is
`5.631402531e-4` versus `5.632820318e-4` for independently tessellated
exact-orbit generators. The native-versus-exact-orbit density residual is
`2.628681493e-7`, and mass-weighted disk position error is `5.79746e-9`.
The native and exact-orbit disk face graphs are identical. Absolute disk
log-entropy change is `1.252197e-4` (signed mean `-4.85972e-6`); this does
not by itself measure irreversible heating. Total energy relative change is
`-2.00210e-11`.

In the outer `r>=2.75` region, mass-weighted absolute radial speed has
grown to `0.05426` and absolute density change to `0.05605`. These outer
errors remain much larger than disk errors and warrant monitoring for
propagation across the disk boundary. This snapshot is still far short of
the historical `t>6` failure window and the paper's `t=600` image.

See the [geometry audit](KEPLER_POLAR_Q1_T061875_AUDIT_2026-10-27.json),
[radial audit](KEPLER_Q1_RADIAL_T061875_2026-10-27.json),
[face audit](KEPLER_POLAR_Q1_FACES_T061875_2026-10-27.json), and
[native density field](KEPLER_POLAR_Q1_NATIVE_DENSITY_T061876_2026-10-27.png).
In the paper's `0.8<R<1.2` annulus, the [L1 audit](KEPLER_POLAR_Q1_GIZMO_METRICS_T061876_2026-10-27.json)
gives velocity `5.23054e-9` and plateau-corrected density `3.98664e-4`.
The literal paper density reference gives `0.01000021` because the actual
prescribed plateau is `1.01`, not `1.0`; these very early values do not
predict either threshold-crossing time.

## Full q2 native continuation toward the historical t>6 window (2026-10-27)

The accepted q2 binary and original full 147,456-site IC are retained for
the next physical-time comparison. The original 8-cell H200 continuation
jobs **436985** (`afterok:436890`, to `t=2.221`) and **436986**
(`afterok:436985`, to `t=6.5`) were staged with independent complete-snapshot
and physics gates. Both were canceled before execution after the later
ballistic halo-width screen showed that 8 cells are insufficient. The
[first original runner](kepler_polar_q2_to2p221_20261027.slurm) and
[second original runner](kepler_polar_q2_to6p5_20261027.slurm) are retained
as provenance. The replacement 28-cell dependency branches and quantitative
screen are documented below.

The [restart preparer](prepare_kepler_polar_long_restart_20261027.py)
requires a passing parent audit, the exact expected snapshot time, a complete
native snapshot, zero cumulative floor/dual-energy injection in the parent
RK4 ledger, and the accepted executable SHA-256 before copying into a
new directory. The [long-time audit](audit_kepler_polar_long_gate_20261027.py)
checks all cell IDs, exact masses, W=0, positive and finite state, zero
floor/dual-energy injection, disk density and entropy drift, orbit and
angular-momentum error, total energy, and the paper's annular L1 gates.
It was exercised on the complete q1 `t=.001` preflight; all checks passed.
The restart preparer was separately exercised with a real q1 native
checkpoint and temporary receipt; source time, binary hash, and copied
inputs were verified; the additional zero-injection parser passed on the
completed q1 native log. Both original queued directories passed
`scripts.sha256` before cancellation. Earlier pending jobs 436983 and 436984 were canceled
before execution (both elapsed `00:00:00`) so this stronger parent gate
could be staged in new directories without modifying queued inputs.
No `t=2.221` or `t=6.5` native result is claimed yet.

## q1 t=0.691 GPU halo abort and independent recovery test (2026-10-27)

Original q1 continuation 436889 ended `FAILED` after 3:20:04 at RK4
step 9071, accepted step-start `t=0.691008245204794`. The k4 tessellation
reported `[GPU_TESS_INVALID] volume=0 velocity_fit=0 geometry=0 halo=4`:
four cells' conservative polygon search radii crossed the 8-cell padded
grid. The code aborts rather than using potentially incomplete periodic
competitors. This is a **GPU halo sufficiency failure**, not evidence of
negative internal energy or physical disk breakup. The preceding RK4 ledger
had `n_ie_le0=0` and zero floor, dual-energy, and entropy-floor injections.
The `khdiag.step009071.start_geometry.dat` file is a diagnostic accepted
step-start state; it is **not** used as a restart.

The [read-only radial audit](KEPLER_Q1_FAILURE_START_RADIAL_T069101_2026-10-27.json)
of that diagnostic has disk radial speeds of order `10^-8`, while the outer
`r>=2.75` region has mass-weighted absolute radial speed `0.06021` and
density change `0.06767`. The [periodic short-pair screen](KEPLER_Q1_FAILURE_START_PAIRS_T069101_2026-10-27.json)
finds 140 short generator pairs, all with endpoints whose initial radii
are at least `2.9187`; none are in the `r<2.35` disk/guard region. This
spatial evidence points to an outer-boundary computational failure, but
does not identify the four clipped cell IDs or prove that all effects are
confined there.

The invalid-state fallback 436895 was canceled before execution because
its safety policy permits only timeout/preemption/node failure, not a
geometry `FAILED` parent. Independent H200 job **436989** instead starts
from q1's last complete, separately audited `t=0.618764447` step-8119
snapshot, retains the exact same binary and physics, and changes only
`GFS_KEPLER_HALO_CELLS=8` to `12`. Its
[recovery preparer](prepare_kepler_q1_halo_recovery_20261027.py) verifies
the parent failure signature, initial state, accepted checkpoint, zero
injection, and binary SHA before copying. The preparer passed a real-data
dry run; the [runner](kepler_polar_q1_halo12_recovery_20261027.slurm) and
auditors are byte-verified in a separate run directory. The original q1
failure output is preserved.

**Intermediate native observation:** the halo-12 recovery accepted restart
step 960 at `t=0.691620`, beyond the original halo-8 abort near
`t=0.691008`. The four preceding tessellations reported zero halo-clipped
cells; the RK4 ledger still had positive internal energies, zero floor and
dual-energy injection, and `dEtot/|Etot0|=-9.492e-12` relative to this
restart's reference. This verifies passage through the former geometric
barrier, not completion of the `t=.825` gate or correction of the periodic
outer-state physics.

**Observed q1 `t=.825` result:** job **436989** completed the native
evolution and saved complete `khout.010856.dat` at exact `t=.825`. Its
original early and long audits passed, but its final radial auditor could
not import a module omitted from that work directory, so Slurm marked the
batch `FAILED 1:0` after the physics output was complete. The saved state
was independently re-audited without rerunning the simulation: the
[strict 15-check receipt](/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/halo12_from061876_to0825/audit_long_t0825_strict.json)
and [radial receipt](/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/halo12_from061876_to0825/audit_radial_t0825.json)
both pass. Initial and final cell areas sum to `35.999999998` and
`35.9999999996`; maximum `rho/(m/V)-1` is `5.92e-8`. The initial-disk
mass-weighted absolute density change is `5.2352e-4`, absolute entropy-log
error `1.5332e-4` (signed mean `-8.45e-6`), position error `1.11e-8`,
and whole-box relative total-energy
change `-3.12e-11`. Floor/dual-energy injections remain zero. The outer
edge-within-0.2 density change is already `0.1697`; halo sufficiency has
been repaired through `.825`, but periodic outer-state physics has not.
The increasing absolute entropy error is a numerical-noise measure here;
the signed result does not show net irreversible heating at this time.

The radial auditor is now self-contained and was rerun in 3.96 s on the
saved checkpoint. The first dependent q1 pair **437208/437209** could never
start because 436989's final packaging error made `afterok` false; the
replacement pending pair **437210/437211** was canceled before GPU use so
its first runner could require the strict 15-check receipt directly. Current
q1 stages **437212** and **437213** target exact `t=2.221` and `t=6.5` in
separate directories. The first has no Slurm parent dependency because
436989 is terminal `FAILED`; its restart preparer requires the passed
strict receipt, complete exact-time checkpoint, accepted binary SHA, and
zero-injection parent ledger before running. The second requires
`afterok:437212`. Both use fixed mass, W=0, variational pressure, RK4,
and a 28-cell periodic halo. Each runner checks staged script hashes and
writes a fresh long-gate and radial audit before dependent work starts. The
[first](kepler_polar_q1_halo28_to2p221_20261027.slurm) and
[second](kepler_polar_q1_halo28_to6p5_20261027.slurm) scripts passed
shell syntax, Python syntax and staged SHA-256 checks. These q1 and q2
branches are independent. Neither pending q1 stage is yet a `t>6` result.

After the q1 postprocessing import failure, the exact staged Python imports
were executed successfully in both pending q1 long-stage directories, both
pending q2 long-stage directories, and the running q2 halo-12 recovery
directory. The pending q2 long-stage `scripts.sha256` checks also passed.
This preflights the previously missing module path without altering any
active physics binary or native checkpoint.

An independent periodic SciPy Voronoi reconstruction of the saved q1
accepted states now sharpens the cause. For cells within 0.5 of a 6x6 box
edge, the largest ratio `2*circumradius/(0.5+physical_edge_distance)` rises
from `0.31615` at t=.2063 to `0.59950` at t=.4125, `0.88180` at t=.6188,
and `0.999876` at the failed step-start t=.6910. The last four near-limit
cells lie symmetrically along the periodic seam near the corners; their
independent polygon area is `0.003708395`, matching the saved native volume
to ~`1e-9`, about 15.2 times the original Cartesian outer-cell area
`0.000244141`. The true periodic cell radius has grown to `0.250105`;
its required W=0 search reach `2R` is almost exactly the 8-cell halo width
of `0.5`. A small k4-stage geometry change can therefore cross the bound
without any negative-volume or energy event. Generator coordinates remain
inside the 6x6 physical box in all saved states. Wider padding permits
this evolving outer geometry to be represented, but does not remove the
periodic seam disturbance.

At t=.2063 the [q1](KEPLER_Q1_PERIODIC_HALO_GEOMETRY_T020625_2026-10-27.json)
and [q2](KEPLER_Q2_PERIODIC_HALO_GEOMETRY_T020625_2026-10-27.json)
outer halo-risk ratios, radius p99, and radius maxima are **identical**.
The q2 outer layout is shared with q1, so q2 is at high risk of the same
8-cell abort near t=.69; this is a prediction, not yet a q2 failure.
The [q1 t=.412](KEPLER_Q1_PERIODIC_HALO_GEOMETRY_T041254_2026-10-27.json),
[t=.618](KEPLER_Q1_PERIODIC_HALO_GEOMETRY_T061876_2026-10-27.json), and
[failed step-start](KEPLER_Q1_PERIODIC_HALO_GEOMETRY_T069101_2026-10-27.json)
receipts preserve the independent measurements. The boundary-taper
control remains necessary to isolate the outer forcing itself.

The [wrapped exact-orbit control](KEPLER_Q1_BALLISTIC_PERIODIC_HALO_T069101_2026-10-27.json)
uses the *initial* q1 generators, moves them on exact softened circular
orbits to `t=0.691008245`, wraps positions into the 6x6 box, and builds an
independent periodic Voronoi diagram with no gas force or heating. It
predicts **four** 8-cell halo-bound violations at the **same four IDs**
`5760, 141695, 147087, 368` that have the largest near-limit radii in
the native accepted step-start reconstruction. Its maximum
`2R/(0.5+edge)` is `1.000416`; the native step-start is `0.999876`, and
the subsequent native k4 tessellation reports four clipped cells.
This strongly attributes the first clipping event to the initial outer
orbital geometry combined with periodic wrapping, rather than accumulated
disk hydro heating. The wrapped exact orbits are a *geometric control*,
not a physical exterior solution across the nonsmooth periodic seam.

## Conditional q2 halo12 continuation (2026-10-27)

Because q1 and q2 have identical early outer-seam geometry, a separate
q2 contingency is queued without changing or canceling the running q2
job 436890. H200 job **436995** requires **both**
`afterok:436989` (q1 halo12 recovery passed) and `afternotok:436890`
(original q2 did not complete). Its [preparer](prepare_kepler_q2_halo_failure_recovery_20261027.py)
then additionally requires the *same isolated GPU halo abort signature*,
zero cumulative floor/dual-energy injection, the accepted binary SHA,
and a complete positive fixed-mass W=0 q2 accepted state. The revised
selector accepts only the specifically verified `start_geometry` diagnostic
that copied the accepted RK4 `k1_input`, never a failed-stage payload; it
starts from that state rather than t=0. If either prerequisite is absent,
this contingency consumes no evolution time. The earlier conservative
selector selected the ordinary step-8119 checkpoint in a q1 regression test;
the q2 accepted-step update is audited below.

The first q2 halo12 recovery remains limited to `t=.825`. Its original
12-cell long continuations **436996** and **436997** were canceled before
execution after the [future-halo width audit](KEPLER_Q2_BALLISTIC_HALO12_T2221_WIDTH_AUDIT_2026-10-27.json)
below. The original 8-cell long continuations **436985** and **436986**
were also canceled before execution. Their submitted scripts are retained
as historical receipts; no result is claimed for them.

Two mutually exclusive 28-cell continuation branches are now staged and
submitted: **437074** (`afterok:436890`) advances the original q2 `t=.825`
checkpoint to `t=2.221`, followed by **437075** to `t=6.5`; **437076**
(`afterok:436995`) and **437077** do the same from a successful halo12
recovery checkpoint if the original q2 job fails. The
[first](kepler_polar_q2_halo28_to2p221_20261027.slurm) and
[second](kepler_polar_q2_halo28_to6p5_20261027.slurm) runners preserve the
accepted q2 binary, fixed mass, W=0, original IC, and RK4 with the
variational W0 pressure operator. AV5 HLLC intermediate values are computed
but do not drive pressure force or work in this branch.
They change only the periodic halo from 8/12 to 28 cells and give the first
leg a 24-hour limit. All four new stage directories passed script SHA and
shell-syntax checks, with H200:1 and 32 GB requested. The parent checkpoint
preparer will require positive full-state checks, zero floor/dual-energy
injection, exact accepted binary SHA, and a complete `t=.825` snapshot.
None of the new long jobs has run yet.

**Observed q2 outcome, 2026-10-27:** Original full-q2 job **436890** ended
`FAILED` at the same physical time as q1,
`t_start=0.69100824520479398`, step 9071. Its final accepted RK4 ledger at
step 9057 has `dEtot/|Etot0|=-5.616e-11`, positive internal energy, and zero
cumulative floor, dual-energy, and entropy-floor injection. Three geometry
stages report `halo_clipped_cells=0`; the fourth reports exactly **4** and
`[GPU_TESS_INVALID] volume=0 velocity_fit=0 geometry=0 halo=4`. The saved
`khdiag.step009071.start_geometry.dat` contains the accepted `k1_input`
state at the failed step start; the failed `k4` geometry is not serialized.
Thus the q1-based outer-halo prediction is now confirmed in
the independent q2 native run; this is a geometric abort, not a measured
late-time disk-heating event.

The exact staged q2 [recovery selector](prepare_kepler_q2_halo_failure_recovery_20261027.py)
was dry-run on this actual failed q2 parent in an isolated temporary
directory. It passed the binary, failure-signature, zero-injection,
completeness, fixed-mass, W=0, and positivity checks, selecting the latest
complete **ordinary** checkpoint `khout.008119.dat` at
`t=0.6187644470884642`. This conservative dry-run preceded the direct
serialization audit below.
Job **436995** initially waited for the independent q1 halo12 proof
**436989**; the first q2 recovery changes only the halo from 8 to 12, and
its preparer repeats these checks in the actual job. The dependency was
subsequently reassessed below.

**Accepted-step restart update (2026-10-27).** Native source review shows
the optional `start_geometry` file is copied from the RK4 `k1_input` payload
before any stage mutation, with the same binary header, particle payload,
and 52-byte trailer as `khout`. `kh_readdata` reads those same fields. The
q2 file has the complete 147,456 IDs, exact original cell masses, W=0,
positive state, area sum `36.000000012`, `rho=m/V` to `5.92e-8`, and disk
mean orbit displacement `1.73e-9`; its time and trailer are
`0.691008245204794` and `dtold=7.6420175901e-5`. Its SHA-256 is
`b393d585c49c7f43a72be1cef4b99a37c331ea8df6d3126dd2bb4c9c825a8dd8`.
The [revised selector](prepare_kepler_q2_halo_failure_recovery_20261027.py)
now requires the exact accepted-step log signature, matching complete
binary trailer, original mass, W=0, positive state, periodic area and
density consistency, then copies this accepted file to the standard name
`khout.009071.dat` in the recovery directory. An isolated dry-run on the
actual failed run selected this file and step 9071. Staged script hashes
passed after replacement; job 436995 was held during the update and then
released. This saves redoing `0.072244` time units, about 945 global steps.
It is still a recovery **attempt**: its native continuation and physical
audit must pass before any later promotion.

**Parallel-release update (2026-10-27).** q1 halo12 had accepted
`t=0.703576`, beyond the original `t=0.691008` halo-8 failure. The staged
q2 scripts passed SHA-256 checks, and the exact staged selector was rerun
on the actual failed q2 parent: it selected the accepted k1-start step
9071 with source SHA-256
`b393d585c49c7f43a72be1cef4b99a37c331ea8df6d3126dd2bb4c9c825a8dd8`.
The q2 physical comparison is independent of q1's final `.825` result,
so the unnecessary `afterok:436989` Slurm dependency was removed from
job 436995. `scontrol show job` then reported `Dependency=(null)` and
`Reason=QOSGrpGRES`: q2 is now independently waiting for its H200 slot.
The q1 `.825` audit remains a separate gate; neither run has yet proved
long-time disk stability.

**Full-layout global-area gate (2026-10-27).** The long-checkpoint
[auditor](audit_kepler_polar_long_gate_20261027.py) now additionally requires
the initial and endpoint Voronoi areas each to sum to the periodic 6×6 box
area within `1e-5`, and maximum `|rho/(m/V)-1|<1e-6` over all IDs. This
closes the global initialization/geometry blind spot found in the earlier
sparse IC. End-to-end replay on the complete full-q2 native `t=.001`
preflight passed all 15 checks: area sums `36.00000000116415` and
`35.99999999787542`, maximum density relation error `5.90153e-8`.
The revised auditor SHA-256 is
`8d882ce39780992f687c4440bd0cd1e90f1ec781d95f4edbaa2c2b1903109346`.
Pending q2 recovery **436995** and its conditional `t=2.221` / `t=6.5`
jobs **437076/437077** were held during deployment; their staged
`scripts.sha256` checks passed and all three jobs were released. Running q1
**436989** retains its startup-verified auditor; apply the strengthened
auditor post hoc to its accepted `.825` endpoint before using the result
as global-area evidence.

The original-q2 `afterok:436890` long branch **437074–437075** can never
run after the verified 436890 failure and was canceled without execution.
The timeout-only fallback **436896** was also canceled: its documented
selector permits TIMEOUT/PREEMPTED/NODE_FAIL/BOOT_FAIL but refuses this
physical geometry `FAILED` state. The separate halo-recovery branch began
as **436995 → 437076 → 437077**. Job 437076 is now running; before 437077
started, it was replaced by the audited
[t=4.5 and t=6.5 segmented chain](KEPLER_Q2_HALO28_SEGMENTED_T65_2026-10-27.md)
**437250 → 437251** to preserve a physical checkpoint and avoid a possible
48-hour timeout. The q2c continuation **436894** was an independent
full-count control and is assessed below.

Because q2c retains the same outer generators, its halo-8
continuation was exposed to the same outer clipping. Its independently
checked initial native snapshot has 147,456 sites, total periodic area
`36.0000000000873`, and total mass `13.1217219`; it does not have the
sparse-mass defect. A separate conditional H200:1 recovery **437179** was
submitted with `afternotok:436894`, using the same accepted binary and
physical settings but halo 12. Its byte-verified
[runner](kepler_polar_q2c_halo12_recovery_20261027.slurm) applies the actual
failure-signature, complete-checkpoint, fixed-mass, W=0, positivity, and
zero-injection selector before any evolution. It refuses a non-halo failure
and never restarts from a failed RK4 stage.

**Observed q2c outcome, 2026-10-27:** 436894 ended at accepted step-start
`t=0.69101477440242343`, step 5681. Its first RK4 geometry stage had zero
halo clips; the next had exactly four, followed by
`[GPU_TESS_INVALID] volume=0 velocity_fit=0 geometry=0 halo=4` and `EXIT:1`.
The native source copies `k1_input` before any stage mutation and writes
`khdiag.step005681.start_geometry.dat` with the standard snapshot header,
payload, and trailer. The failed stage geometry is not serialized. The
staged selector (SHA-256 `53a93cb3a18c01bec926b44f781be22861ad1b020d88f70476e5342caae6ffce`)
checks that accepted-step signature, original ID/mass completeness, W=0,
positive state, periodic area, and density consistency. An isolated dry-run
on the actual failed q2c parent selected the accepted file with SHA-256
`ad89c84b42994043f4fef79d5e4c399da004480bad1fd7ce789399b21bf67bfb`,
step 5681, under the standard `khout.005681.dat` name. It avoids repeating
the interval from the last ordinary output at `t=0.6187947519066646`.
The following endpoint assessment resolves this conditional recovery.

**Resolved q2c endpoint, 2026-10-27:** recovery 437179 advanced the accepted
step-5681 state to exact `t=.825`, step 6821, and the native executable wrote
`EXIT:0`. The batch job reported `FAILED 1:0` only because its final radial
analysis copy imported an unstaged module. The already completed native
[strict receipt](KEPLER_Q2C_T0825_STRICT_2026-10-27.json) passes all 15
finite/area/fixed-mass/W=0/positivity/no-injection/accuracy gates. CPU job
437248 recomputed the missing [radial receipt](KEPLER_Q2C_T0825_RADIAL_2026-10-27.json)
from the original initial and exact-time endpoint snapshots in one second;
all five radial-state checks pass. No GPU restart was necessary.

At the common `.825` physical time, q2c's original disk has 49,007 cells
versus q2's 49,038. Its mass-weighted absolute density change is
`1.57610e-4` versus q2 `1.16521e-4` (1.353 times larger). The independent
[q2c exact-orbit Voronoi control](KEPLER_Q2C_BALLISTIC_2026-10-26.json)
predicts `1.57660e-4`; the native aggregate differs by only `4.99e-8`.
Q2c absolute log-entropy error is `1.03994e-4` versus q2 `1.29736e-4`
(0.802 times), while its signed mean is `-7.01e-6` versus `-1.90e-5`.
Q2c took 6,821 steps versus q2 10,856, so this entropy difference does
not isolate center geometry from timestep effects. Removing the central
generator lowers the early absolute entropy error but worsens the early
disk density geometry; it is not a demonstrated cure for `t>6` heating.

A separate CPU-only Slurm job **437249** completed an ID-matched
[q2c native-versus-exact-orbit Voronoi audit](KEPLER_Q2C_NATIVE_BALLISTIC_T0825_2026-10-27.json)
in three seconds, using the two original q2c snapshots. All five native
state checks pass. In the initial disk, the independent geometry-only
absolute density change is `1.576602e-4`, versus native `1.576103e-4`.
The **cellwise** mass-weighted absolute native/geometry density residual
is `1.05676e-7` (95th percentile `3.40073e-7`), about `0.0670%` of the
native density-change norm. The same q2 residual is `8.76967e-8`.
Thus the q2c density result is geometry dominated in a cellwise residual
norm, not merely an accidental agreement of aggregate means. This early unweighted
Voronoi control does not test gas stability after `t>6`. An independent
[q2c periodic-image certificate](KEPLER_Q2C_PERIODIC_SAFETY_T0825_2026-10-27.json)
from CPU job **437252** checks the complete native initial state: initial
area sums to `36.00000000008731`, `rho=m/V` holds within `5.78e-8`, and
the minimum omitted-image safety ratios for all 49,007 original disk
polygons are `40.58` at `t=0` and `40.56` at exact-orbit `t=.825`.
This validates the early disk polygon comparison against the omitted
periodic images; it does not make the outer periodic flow a steady IC.
The [endpoint report](KEPLER_Q2C_NATIVE_T0825_RESULT_2026-10-27.md) also
records the post-run packaging failure, the independent EOS check, and the
[common-scale q2/q2c native density fields](KEPLER_Q2_Q2C_NATIVE_DENSITY_T0825_2026-10-27.png).

The [q2c accepted-state radial audit](KEPLER_Q2C_FAILURE_START_RADIAL_T069101_2026-10-27.json)
passes complete-ID, exact fixed-mass, W=0, positive-state, and density
consistency gates at `t=0.6910147744`; its static total-energy change is
`-2.78018e-10` relative. Compared with q2 at `t=0.6910082452`, the
mass-weighted absolute entropy-log change in the initial `r<0.1` bin is
`5.62231e-5` versus `9.53822e-6`, and mean orbit position error is
`1.61714e-5` versus `1.10435e-5`. The q2c bin has 158 cells versus 129
for q2 and different enclosed mass, so this short comparison does not
establish a converged or generally worse center treatment. It does show
that removing the central generator has not eliminated early core error.

The [accepted q2 step-start radial audit](KEPLER_Q2_FAILURE_START_RADIAL_T069101_2026-10-27.json)
adds spatial evidence at `t=0.691008`: all 147,456 IDs are complete,
fixed mass and W=0 hold exactly, volumes and internal energies are positive,
and the static total-energy change is `-5.59e-11` relative. No initial disk
cell crossed either disk edge. At initial `0.5<=r<0.75`, the mass-weighted
absolute radial speed is `7.87e-9`; for initial `r>=2.75`, it is `0.06021`
with `0.06767` absolute density change. Cells currently within 0.2 of a
periodic box edge have `0.15339` absolute density change and `0.14812`
absolute radial speed. This supports the outer-boundary interpretation; it
does not establish what happens after recovery or at `t>6`.

## Ballistic future-halo width screen (2026-10-27)

To avoid starting an under-padded full native continuation, the independent
periodic Voronoi auditor wrapped the original q2 generators after exact
softened circular motion to `t=.825`, `2.221`, and `6.5`. It used a
2.5-unit periodic image margin for the two later times: the margin exceeds
twice the largest reconstructed edge-cell circumradius by factors 1.53 and
1.21 respectively. The relevant inequality is
`2*circumradius < distance_to_physical_edge + halo_width`. Among audited
edge cells:

| Exact-orbit time | Minimum integer halo cells | Violations at 8 | At 12 | At 24 | At 28 |
|---:|---:|---:|---:|---:|---:|
| 0.825 | 10 | 24 | 0 | 0 | 0 |
| 2.221 | 25 | 256 | 200 | 24 | 0 |
| 6.5 | 24 | 200 | 140 | 0 | 0 |

The [early](KEPLER_Q2_BALLISTIC_HALO12_T0825_WIDTH_AUDIT_2026-10-27.json),
[t=2.221](KEPLER_Q2_BALLISTIC_HALO12_T2221_WIDTH_AUDIT_2026-10-27.json),
and [t=6.5](KEPLER_Q2_BALLISTIC_HALO12_T65_WIDTH_AUDIT_2026-10-27.json)
receipts include width scans, image-margin checks, and risk-cell IDs. This
is a geometry-only risk screen, not a valid exterior gas solution after
crossing the nonsmooth periodic seam. A 28-cell halo is adequate at these
sampled ballistic times; it is **not** a proof of sufficiency at all
intermediate times or in the native gas run. The continuation gates remain
necessary, and a wide global halo is a diagnostic cost, not a satisfactory
long-term performance strategy.

## Same-time native q2 checkpoint at t=0.412548 (2026-10-27)

Original q2 job **436890** wrote complete fixed-mass W=0 snapshot
`khout.005403.dat` (147,456 cells, 84,934,724 bytes). The independent
[checkpoint audit](KEPLER_POLAR_Q2_T041254_AUDIT_2026-10-27.json) finds
finite state, exact fixed masses and W=0, positive volumes and internal
energies. For initial `0.5<=r<2`, its mass-weighted mean absolute native
density change is `1.296369e-4`, versus `1.296449e-4` from the exact-orbit
Voronoi control. Their cellwise absolute residual averages `2.81e-8`.
The same-time [q1 audit](KEPLER_POLAR_Q1_T041254_AUDIT_2026-10-26.json)
has native density change `5.709482e-4`: q2 is smaller by a factor of
`4.40` at this early time. The q2 absolute entropy-log error is
`6.58034e-5` versus q1 `9.58481e-5`; q2 signed entropy-log change is
`-5.15230e-6`. These entropy numbers are error scales, not measured
irreversible heat production.

The [q2 face audit](KEPLER_POLAR_Q2_FACES_T041254_2026-10-27.json) finds
its native and exact-orbit disk edge graphs identical (145,671 faces),
despite 43,667 lost and 43,707 new faces relative to t=0. Its
[paper-annulus audit](KEPLER_POLAR_Q2_GIZMO_METRICS_T041254_2026-10-27.json)
at `0.8<R<1.2` gives `L1_v=1.0654e-9` and density deviation from the
actual initial 1.01 plateau `1.2356e-4`. The
[radial audit](KEPLER_Q2_RADIAL_T041254_2026-10-27.json) reports
outer `r>=2.75` mean absolute radial speed `0.03645`, while the disk's
radial speed remains at the `1e-9` level. The
[native density field](KEPLER_POLAR_Q2_NATIVE_DENSITY_T041254_2026-10-27.png)
shows the ring-scale residual with the same normalization as the earlier
q1 image.

The [independent periodic-halo reconstruction](KEPLER_Q2_PERIODIC_HALO_GEOMETRY_T041254_2026-10-27.json)
has no 8-cell violation yet, but its maximum risk ratio is `0.59949567`,
matching q1 at the same time. Thus q2's disk geometry is measurably better
while the shared outer periodic-seam risk persists. This checkpoint proves
neither reduced irreversible heating nor survival to `t>6`.

## Exact native q2 t=0.825 gate (2026-10-27)

The halo-12 recovery **436995** completed at exact `t=.825` with all 15
full-state checks passing and zero floor/dual/entropy injection. Its
[same-time q1/q2 result](KEPLER_Q1_Q2_NATIVE_T0825_RESULT_2026-10-27.md)
records q2 disk absolute density drift `1.16521e-4`, versus q1
`5.23520e-4`. Independent exact-orbit Voronoi reconstruction predicts the
q2 drift `1.16545e-4`, with cellwise native residual `8.77e-8`; its central
disk polygons pass a periodic-image safety factor of at least `41.08`.
q2's signed mean log-entropy change is `-1.90e-5`, so the early density
advantage is an aligned-geometry result, not evidence of suppressed late
irreversible heating. The dependent halo-28 path to `t=2.221`, then `t=4.5`
and `t=6.5` still requires its own native gates; the
[segmented-chain record](KEPLER_Q2_HALO28_SEGMENTED_T65_2026-10-27.md)
has the verified dependencies and runtime headroom calculation.

## Current q1 native continuation observation (2026-10-27)

The q1 H200 continuation **437212** remains RUNNING. Its latest log sample
reached `t=0.959996`, global step 12,657, with `dt=7.64202e-5`. The latest
printed RK4 ledger reports `dEtot/|Etot0|=-2.360e-11`, `dEhyd/Ehyd0=2.255e-5`,
`n_ie_le0=0`, and zero floor, dual-energy, or source-floor injection. Its
exact `t=2.221` gate has not run; dependent q1 stage 437213 remains pending.

The subsequent log read reached `t=0.985041` / global step 13,001. Its latest
ledger reports `dEtot/|Etot0|=-2.756e-11`, `dEhyd/Ehyd0=2.716e-5`, no
nonpositive internal-energy cells, and zero floor/dual/source-floor
injection. The t=2.221 q1 audit remains outstanding.

At 21:17 KST, job **437212** remained RUNNING and had advanced to `t=1.1244`,
global step 14,861. Its latest ledger reports `dEtot/|Etot0|=-4.578e-11`,
`dEhyd/Ehyd0=5.542e-5`, no nonpositive internal-energy cells, and zero
floor/dual/source-floor injection. The exact `t=2.221` q1 endpoint audit is
still pending.

At the 21:34:51 live log check, q1 job **437212** had reached `t=1.15404`,
global step 15,261. The latest RK4 ledger reports
`dEtot/|Etot0|=-4.225e-11`, `dEhyd/Ehyd0=6.227e-5`, no nonpositive internal
energy, and zero floor/dual/source-floor injection. The exact `t=2.221`
endpoint audit and dependent `t=6.5` run are still pending.
