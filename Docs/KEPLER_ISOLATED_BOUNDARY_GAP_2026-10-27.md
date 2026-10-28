# Kepler benchmark boundary gap, verified 2026-10-27

## What the published comparison actually specifies

[Hopkins (2015), §4.2.4 and its resolution footnote](https://arxiv.org/html/1409.7395#S4.SS2.SSS4)
specify an effective 256² particle resolution over `[-2,2]²`, with the
inner and outer particle populations reduced to follow the prescribed
surface-density profile. The footnote explicitly distinguishes particle
methods' ability to use a vacuum edge from the softened low-density edge
introduced for mesh-method comparison. Its Figure 7 shows MFM at `t=600`.
The text also reports evolution to approximately 1,000 inner-disk orbits
without showing a corresponding image or an unambiguous reference radius.
At `r=.5`, the **unsoftened** nominal 1,000-orbit time is `2221.441469`;
the present Plummer force with `epsilon=.01175` gives `2222.361499` for
1,000 circular periods. Thus `t=600` is the plotted comparison and the
longer orbit count is a separate final duration claim to interpret and
verify. The earlier wording here incorrectly attributed `2221.44` to
the softened frequency.

There is a radius/count ambiguity within the paper's own prose: the Figure 7
caption says the disk at `r~0.5` has completed `>250` orbits by `t=600`,
whereas the paragraph below says `t=600` is about `160` inner-disk orbits.
From the stated softened circular speed, `P(0.5)=2*pi*(0.5^2+epsilon^2)^(3/4)`
is about `2.22236`, so `t=600` is about `270` orbits at `r=.5`, consistent
with the figure caption. Keep the plotted endpoint exact at `t=600`; report
orbit counts at explicit radii. Interpret the prose's approximate 1,000
orbits as an additional duration target, with `t=2222.361499` at `r=.5` as
the conservative radius-explicit realization, rather than silently equating
it to the Figure 7 snapshot.
This passage does **not** specify the 6×6 periodic torus used by the current
native q1/q2 runs. It also does not establish the exact outer boundary
algorithm for our Voronoi code, so an isolated boundary must be validated
on its own terms before calling the native result a matched replay.

The current full 384² q1/q2 run fills the 6×6 box with a circulating
low-density exterior and wraps both coordinates. Its outer velocity and
gravity derivatives are discontinuous across the identified seams. The
[seam audit](KEPLER_PERIODIC_SEAM_CONTROL_2026-10-26.md) finds a median
cross-seam velocity jump about 77 times the local sound speed. Enlarging
the GPU halo corrects missing periodic competitors; it does not make this
outer state a stationary Kepler solution on a torus. At native q1 `t=.825`,
the mass-weighted absolute density change is `0.1697` for cells within 0.2
of a box edge, versus `5.24e-4` in the initial disk. This is spatial
evidence of an exterior problem, not proof of causal penetration into the
disk at late time.

The corrected sparse q2 IC is closer to the paper's particle-count intent:
it preserves 49,038 disk sites, compared with about 48,255 uniform sites
expected from the paper's effective 256² density over the annulus. This
1.6% count comparison is an inference from area, **not** an exact IC match.
It still uses the same periodic outer rule and has its separate sparse-core
gravity-resolution limitation.

## Native implementation audit

- The periodic-baseline `Exam/KH/kh.c` Kepler branch selects
  `periodic_postStage_blend`, which wraps x and y after every stage. The
  fixed-polygon run is a separate opt-in branch.
- The available `mkLinkedList2D_sod` path reflects at x walls and remains
  periodic in y; `walls_xy_postStage_blend` reflects at all four walls for
  the Sedov case. Neither is a vacuum/outflow Kepler boundary.
- The periodic-base CUDA `GFS_VARPRESS_W0` kernel rejects unsupported
  ghost-neighbor faces. This restriction applies to its periodic face path;
  it does not describe the later fixed-polygon overlay.
- The separate
  [fixed-polygon source overlay](KEPLER_FIXEDPOLY_SOURCE_OVERLAY_2026-10-27/README.md)
  clips each cell to a 1,024-sided polygon and represents polygon edges with
  explicit boundary-face records. With `GFS_VARPRESS_W0=1`, it adds the
  static-wall relative-sweep pressure work and uses the bounded-area pressure
  gradient, without a ghost-pair force. Its manifest records the build mode,
  and the `t=30` runtime log prints the variational-W0 activation marker.
- This is a closed, fixed wall. A generator crossing still aborts; the
  `30→60` continuation stopped at `t=48.3753` on an outer-cell crossing.
  The branch has no vacuum or outflow boundary and no validated post-impact
  wall impulse/torque ledger.

The periodic-base details come from the preserved
[native source overlay](KEPLER_NATIVE_SOURCE_OVERLAY_2026-10-27/README.md);
the bounded-wall details come from the separate fixed-polygon overlay and
the recorded runtime receipts. The running and queued H200 jobs retain their
current periodic physics and remain controls for that boundary condition.

## Required isolated-boundary branch

1. Keep the already accepted q1/q2/sparse periodic states as controls and
   continue live chains only while their geometry and physics gates pass.
   Their disk, core and edge measurements determine whether the torus
   artifact reaches the disk and whether sparse sampling is acceptable.
   The sparse 32-cell-halo chain ended at `t=2.75514`; its update below
   supersedes the earlier expectation of a `t=6.5` endpoint.
2. Implement a separate opt-in finite-domain/free-surface or outflow path in
   the user's CUDA+MPI+OpenMP solver. Define what happens at every boundary
   face, including its geometric speed, pressure force, thermal work, mass
   policy and escaping energy/angular-momentum flux. The W=0 variational
   operator must be adjoint to the actual bounded-cell volume derivative.
   An omitted-neighbor face or reflected wall alone is not sufficient.
3. Preserve identical disk generators, initial density and pressure when
   comparing periodic and isolated exteriors. Validate initial cell areas
   against independent bounded Voronoi geometry, fixed masses where claimed,
   positive states, and a boundary ledger that closes total energy and
   angular momentum including allowed outward flux. Exercise static uniform
   pressure, free expansion, a circular-orbit control, and a short native
   Kepler evolution before promoting a long run.
4. Measure wall time and GPU/MPI scaling of the new branch. Boundary work
   should scale with boundary faces, not require a whole-domain subcell
   quadrature or repeated global solve. Only after matched `t>6`, `t=120`
   and `t=600` density/orbit/entropy results should it be compared with
   the paper's longest image. Continue toward the paper's reported
   approximately 1,000 inner-orbit duration only after that plotted gate.

An opt-in one-rank native fixed-polygon result now exists through `t=.001`;
the later section records its limits. No accepted long-time isolated-boundary
result exists yet.

## Updated direction after sparse halo-32 failure (2026-10-27)

The patched sparse A10 continuation **437279** stopped at accepted
`t=2.755142729` because four outer cells exceeded the 32-cell periodic halo
in the following RK4 stage. Its [accepted-state audit](KEPLER_SPARSE_T251_GEOMETRY_ABORT_2026-10-27.md#halo-32-failure-at-accepted-t2755142729-2026-10-27)
still passes disk physics and exact-orbit geometry. Independently, 36
initially seam-near cells have entered the `r<2.35` guard, and an outer cell
sets the global hydro step (`0.000979`, half the inner-core orbital bound).
Thus the earlier step 1 expectation that this *particular* sparse periodic
chain would supply a `t=6.5` control is superseded. Retain its complete
accepted state; do not treat a larger halo as a cure for the non-equilibrium
outer torus. The separate full-q2 H200 chain remains an active periodic
control.

### Fixed convex boundary: exact W=0 operator contract

For a fixed convex domain `Omega`, let each cell be
`Omega intersection Voronoi_i`. A static boundary segment has laboratory
normal speed zero. Interior face `ij` retains the current W=0 geometric
speed and pair pressure force. Since boundary segments do not move, the
geometric volume derivative is the sum over **interior** faces only. The
existing CUDA force term `(P_i-P_j)*b_ij` remains the negative gradient of
the cell internal-energy potential. Its summed force equals minus the
pressure traction on the fixed wall, so fluid momentum need not be conserved
by itself; wall impulse must enter the momentum ledger.

The periodic-only CUDA **energy work** does require an extra term. It uses
the interior relative sweep `sum A*(w_n-v_i.n)`, relying on the periodic
cell's interior face vectors closing to zero. At a bounded cell they do not:
their sum is minus the boundary face vectors. For every static boundary
segment with outward normal `n` and length `A`, add the relative sweep
`-A*(v_i.n)` and therefore add `+P_i*A*(v_i.n)` to the current total-energy
RHS (and the equivalent contribution to any authoritative thermal64 path).
No boundary mass flux or ghost force is needed for this *closed, fixed-wall*
model. A specular site/velocity wall rule and wall impulse/torque ledger
remain necessary if a generator reaches the wall; a regular polygon only
approximates a circular slip boundary and must be converged in side count.

The independent [16-cell, 64-sided mathematical audit](KEPLER_BOUNDED_VARPRESS_ADJOINT_SMALL_MESH_2026-10-27.json)
and [reproducer](audit_kepler_bounded_varpress_adjoint_20261027.py) check
two random meshes with uniform and nonuniform pressure. Bounded areas tile
the polygon to `3.6e-15`, reciprocal face geometry to `1.3e-15`, geometric
area rates to `4.1e-9`, pressure force against `-grad(U)` to `5.8e-8`, and
wall momentum traction to `3.2e-15`. The unchanged periodic energy-work
formula has residual power `0.44-0.76` in these tests; the static-boundary
correction reduces it below `8e-16`. This establishes a bounded-domain
**discrete algebraic contract**, not a native CUDA or long-time solution.

### Matched bounded q2 initial generators and areas

The [bounded-IC builder](make_kepler_bounded_polygon_ic_20261027.py) now
preserves every accepted sparse-wide generator with initial `r<2.65`, in
the original order, and drops only the 580 sites farther out. The result
has **55,351** sites, including the identical **49,038** original disk sites
and **1,476** core sites. Its fixed domain is a 1,024-sided polygon centered
at `(3,3)` with circumradius 2.75; every generator begins at least 0.103
inside that radius. No velocity taper, steering, weight change or viscosity
has been introduced. This is a new finite-domain physical problem, with
an unmodified disk IC, for isolating the periodic seam's effects.

The [glass](KEPLER_Q2_BOUNDED_R275_N1024_2026-10-27.glass.bin),
[one-time W=0 area sidecar](KEPLER_Q2_BOUNDED_R275_N1024_2026-10-27.areas.bin)
and [manifest](KEPLER_Q2_BOUNDED_R275_N1024_2026-10-27.json) are complete.
Independent enclosing guard rings at radii 4 and 5 give a maximum area
disagreement `1.18e-10` relative; all 55,351 clipped areas sum to the
polygon area `23.758145361526` within `8.3e-13`. Only 72 cells actually
touch the boundary. An independent byte comparison confirmed all retained
coordinates equal the accepted sparse IC, and the disk initial areas differ
from the accepted periodic sidecar by at most `1.18e-10` relative. This
validates the **initial IC geometry**. The native first-step checks below
have since exercised CUDA geometry and energy work. Generator impact handling,
MPI padding and long evolution remain open before a physical stability claim.

### Native bounded preflight, wall faces, and first continuation

The one-rank A10 native CUDA W=0 preflight `437295` evolved the matched
55,351-site IC to exact `t=.001` with the fixed 1024-gon, fixed mass,
`P=10^-6`, no steering or viscosity. The executable itself exited zero;
its Slurm status was failed only because the first postrun auditor tried to
serialize a NumPy boolean. Rerunning that corrected auditor on the unchanged
snapshots gives [14 passed checks](KEPLER_Q2_BOUNDED_NATIVE_T001_V2_2026-10-27.json).
The additional [independent evolution audit](KEPLER_Q2_BOUNDED_NATIVE_EVOLUTION_T001_2026-10-27.json)
passes 17 checks: total relative energy change `3.62e-14`, mass-weighted
disk absolute density change `2.04e-6`, exact-orbit geometry prediction
`2.04e-6`, native/ballistic cellwise mean residual `1.51e-8`, and no
injection or wall crossing. These are **one-step** results only.

A separate accepted native CUDA wall-face capture `437304` completed and
produced 1,096 physical faces on all 1,024 polygon sides, owned by 72 cells.
The [direct wall-face audit](KEPLER_Q2_BOUNDED_NATIVE_WALLFACES_T0_2026-10-27.json)
checks the support plane (`8.88e-16` maximum), outward area-vector direction
(`2.90e-14`), relative sweep (`8.67e-19`) and pressure-work identity (zero
recorded residual). The uniform-pressure whole-wall traction and torque
cancel to roundoff. This validates the initial static-wall CUDA work contract;
it **does not test a site-wall collision**, reflection, or the wall impulse
and torque ledger after impact.

The first continuation `437302` from the accepted `t=.001` state reached an
accepted step start `t=.045616533734626756`, then the next RK4 stage hit
`[GPU_CLIP_CAPACITY]` with 48 intermediate polygon corners. The complete
accepted state is
`/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/to0825_v2/khdiag.step000024.start_geometry.dat`
(SHA-256 `08a7500e5d72c965c771ad51fd90d87a3cd7872a439edfb146df7a566d18c76a`).
The failed RK stage is not a restart state. The geometric failure is not a
demonstration of disk heating or wall crossing. We moved polygon clipping to
immediately after the eight nearest real halfplanes, before the remaining
competitor stream; this retains a conservative superset for the search while
avoiding irrelevant intermediate supports near the boundary. Native candidate
binary SHA-256 `2c4ac265d9bb9e181bc4af2bb5412af65eefb9a45cc0a3ad6c38dfe5bef13327`
was built from the preserved
[source overlay](KEPLER_FIXEDPOLY_SOURCE_OVERLAY_2026-10-27/MANIFEST.json).
Continuation `437310` did not start physics: the restart argument `1` requires
the accepted file to be named `khout.000001.dat`; the first staging called it
`khout.000024.dat`. The corrected-name `437311` reached the initial tessellation
but a 48-corner wall clip overflowed. Raising the host/device polygon buffer
to 64 (`437313`) moved the overflow to a subsequent real-neighbor clip at the
same accepted start. The [independent accepted-state reconstruction](KEPLER_Q2_BOUNDED_CORNER_CAPACITY_T0045616_2026-10-27.json)
shows the **final** clipped cells need at most 30 corners; two guard radii agree
to `1.18e-10`, and native areas agree within `5.84e-8` relative. These aborts
therefore identify avoidable *intermediate clipping complexity*, not a genuine
64-sided physical cell.

The next candidate keeps the 64-corner abort but preclips 24 more nearest real
bisectors, in addition to eight directional anchors, before the fixed wall.
The complete seed/adaptive competitor stream still follows. Native binary
SHA-256 `3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e`
completed as `437315` from the exact same accepted state at `t=.825`
(Slurm `COMPLETED 0:0`, 11m16s A10). The [17-check native audit](KEPLER_Q2_BOUNDED_NATIVE_T0825_2026-10-27.json)
passes fixed mass, W=0, positive EOS, polygon area, zero injection, no wall
crossing, and total relative energy error `6.67e-10`. Disk mass-weighted
absolute density change is `1.165206e-4`, entropy-log change `2.51848e-6`,
and native versus exact-orbit density residual `8.47e-8`.

The [same-generator periodic comparison](KEPLER_Q2_BOUNDED_VS_PERIODIC_DISK_T0825_2026-10-27.json)
matches all 49,038 initial disk sites exactly. At identical `t=.825`, bounded
versus sparse-wide periodic disk mean absolute differences are `1.20e-11` in
relative density change, `9.65e-12` in entropy-log change and `3.25e-14` in
position. The [density-field image](KEPLER_Q2_PERIODIC_BOUNDED_NATIVE_DENSITY_T0825_2026-10-27.png)
uses a common paper-style `[-2,2]^2` view and color scales. The result shows
that the boundary change has not measurably affected this initial disk at
`.825`; it does **not** prove boundary independence at `t>6`.
The [matched radial-bin audit](KEPLER_Q2_BOUNDED_VS_PERIODIC_RADIAL_T0825_2026-10-27.json)
locates the difference: the original 1,476 core sites have identical density
and entropy metrics; the 4,580 cells in `2<r<2.35` differ in relative density
change by a mass-weighted `2.38e-6`; the `2.35<r<2.65` exterior differs by
order `0.02-0.03`. Those exterior cells also have different initial masses
because the polygon truncates their Voronoi areas, so the comparison does not
rank one exterior as physically superior.

The exact accepted `t=.825` endpoint is now being continued by A10 job
`437317` toward `t=6.5`, with intermediate `.825`-spaced native dumps and
an endpoint physics gate. No earlier accepted state was reset or discarded.

### Independent IC interpretation (Astra consultation, 2026-10-27)

Astra judged the polar q2 and bounded derivative **conditionally admissible**
for a controlled, closed-wall Kepler experiment. `dr/chord≈2.03` is a valid
resolution choice; fixed one-time `m_i=Sigma(r_i)V_i(0)` and autonomous W=0
evolution introduce no ring constraint. But the IC is highly ordered: an
independent 5% angular site jitter raises the exact-orbit Voronoi geometry
error about 7.4 times. The q2 success cannot alone establish general solver
accuracy. The [original paper](https://arxiv.org/html/1409.7395#S4.SS2.SSS4)
also distinguishes special cylindrical alignment from its general Kepler
benchmark. The `r=2.65` site cutoff leaves no vacuum gap because the
remaining exterior cells fill the polygon to `R=2.75`.

At the polygon wall, a circular Kepler flow is not exactly tangent to each
straight side. The maximum continuum normal component is approximately
`v_K(2.75) sin(pi/1024)≈1.85e-3`; this is a geometric estimate, not a measured
native boundary velocity. A wall-side-count test and exterior-sampling test
are therefore required. The fixed wall's relative-sweep correction is a
consistent discrete energy-work term, but native wall-face checks do not
validate collision reflection or the integrated wall impulse and torque.
Comparison with GIZMO Figure 7 may be made only with the distinct closed
boundary and special polar layout explicitly identified. A direct replay
claim requires a matching initial layout and boundary contract.

The [5% independent site-jitter bounded IC](KEPLER_Q2_BOUNDED_JITTER005_R275_N1024_2026-10-27.json)
now tests that alignment dependence with the *same* 1024-gon and exterior
sites. It changes angles of 50,481 polar-ring sites, preserves their radii,
retains 49,038 disk sites, and recomputes one-time W0 masses from independently
clipped areas. The two guard constructions agree to `1.15e-10` relative and
tile the polygon. Native A10 preflight `437319` completed at `t=.001` with
all 14 geometry/EOS/fixed-mass/energy checks; its endpoint SHA-256 is
`6c54048b0dd9effd7a548a479de0616d5fb126697c8d43d5f74a57c9fbe04fa1`.
The same-binary `t=.825` continuation `437320` completed (`0:0`, 11m39s A10).
Its [17-check native audit](KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T0825_2026-10-27.json)
passes fixed mass, W=0, positive EOS, polygon tiling, zero injection, and no
wall crossing. The jitter control changes the generator configuration and initial
cell masses as a finite-volume discretization; it is not a same-site
trajectory comparison or a change in the target continuum Kepler profile.
The [independent exact-orbit Voronoi audit](KEPLER_Q2_BOUNDED_JITTER005_BALLISTIC_T0825_T65_2026-10-27.json)
predicts disk mean absolute density drift `1.6210e-3` for jitter versus
`1.1654e-4` for aligned at `.825` (**13.91 times**); at `t=6.5` it predicts
`1.4788e-3` versus `1.8573e-4` (**7.96 times**). These are geometric
predictions for sites on exact circular trajectories, not a native result or
evidence of entropy production. They show that the aligned q2 IC strongly
reduces partition noise even before considering the hydrodynamic operator.

At the same native `t=.825`, the jittered disk mean absolute density change
is `1.61996e-3` versus `1.16521e-4` aligned (**13.90 times**). The respective
exact-orbit geometry predictions are `1.62104e-3` and `1.16545e-4`;
native-minus-geometry absolute residuals are `2.15725e-6` and `8.47182e-8`.
Disk absolute log-entropy changes are `4.17001e-6` and `2.51848e-6`, both
small at this early time. Thus the large density contrast is already present
in the geometry-only trajectories; it cannot be attributed solely to native
heating. The matched jitter follow-on `437327` is submitted with
`afterok:437320` and an independent parent-audit gate toward `t=6.5`.

### Accepted intermediate restart and hardware timing control

The bounded A10 chain `437317` saved a complete accepted step at
`t=1.651374332206546` (`khout.000427.dat`, SHA-256
`b3993e3569307fceab64a5c53cb3a46fe69edac6abaf979fd73010a5337d649f`).
The [17-check direct native audit](KEPLER_Q2_BOUNDED_NATIVE_T165137_2026-10-27.json)
passes fixed mass, W=0, area tiling, EOS, zero floor injection and no wall
crossing. Its disk mass-weighted absolute density change is `1.24654e-4`,
independent exact-orbit geometry predicts `1.24705e-4`, and their cellwise
mean absolute residual is `3.17532e-7`. Disk entropy-log absolute change is
`4.07983e-6`; total energy relative error is `1.33936e-9`. The A10 chain
continues toward `t=6.5` from the accepted `t=.825` state without resetting.

H200 job `437323` independently restarts the byte-identical `t=.825` state
with the same executable SHA-256
`3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e`
and otherwise identical physics and two CPU threads. Its exact stop is the
A10 accepted time `1.651374332206546`. The [runner](kepler_q2_bounded_h200_speed_to165137_20261027.slurm)
and [immutable-input stager](stage_kepler_bounded_h200_speed_20261027.py)
make this a controlled wall-time and physics comparison; no H200 outcome is
claimed until the native endpoint audit passes.

The aligned `t=6.5` job `437317` also has conditional follow-on `437328` to
`t=10`. Both follow-ons use the [parent-audit restart gate](prepare_kepler_bounded_followon_20261027.py):
they hash the accepted endpoint, check all parent physics gates, preserve the
same CUDA binary, and copy that checkpoint into a fresh run directory only
after the prerequisite job succeeds. A completed `t=.825` parent was used
to preflight the helper and its full input-hash check on a small temporary
stage before submission.

### Independent accepted intermediate physics checks

The jittered native chain `437327` saved the **same exact accepted time**
`t=1.651374332206546` as the aligned chain. Its
[17-check endpoint audit](KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T165137_2026-10-27.json)
passes. At that time the disk mass-weighted absolute density change is
`1.46650e-3` jittered versus `1.24654e-4` aligned (11.76 times), while the
independent exact-orbit Voronoi predictions differ by 11.78 times. The
native-minus-geometry density residuals are `9.10975e-6` and `3.17532e-7`;
the jittered residual is larger but remains below 1% of its own absolute
density drift. Absolute disk log-entropy changes are `6.22665e-6` and
`4.07983e-6`. Their **signed** mass-weighted means are respectively
`-2.61e-8` and `-1.93e-6`, so the measured early contrast is not evidence
of positive net shear heating.

The aligned chain also saved exact `t=3.300240722632501`. Its
[17-check endpoint audit](KEPLER_Q2_BOUNDED_NATIVE_T330024_2026-10-27.json)
passes fixed mass, W=0, area tiling, EOS, energy, disk orbit, zero injection
and no wall crossing. Disk absolute density drift is `1.51927e-4`, versus
`1.51996e-4` from exact-orbit geometry; their cellwise mean residual is
`1.04511e-6`. Absolute disk log-entropy change is `5.92553e-6`, and its
signed mass-weighted mean is `-3.86317e-6` (the initial `0.5<=r<0.75`
band is `-1.30634e-5`). This indicates a small numerical entropy bias
toward cooling in the aligned layout at this time; the cause and its
behavior past `t=6` remain unresolved.

Both runs also saved equal native time `t≈2.4758079598036` (difference
`3.51e-14`). Their separate [aligned](KEPLER_Q2_BOUNDED_NATIVE_T247581_2026-10-27.json)
and [jittered](KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T247581_2026-10-27.json)
17-check audits pass. Disk absolute density change is `1.43447e-4`
aligned versus `1.47095e-3` jittered (10.25 times); the exact-orbit
geometry predictions differ by 10.24 times. Native-minus-geometry
cellwise mean residuals are `6.50583e-7` and `2.25052e-5`, so the
jittered residual is growing but still about 1.5% of its total mean
density change. Disk absolute log-entropy changes are `5.04334e-6` and
`7.95049e-6`. Their signed means are `-3.12642e-6` aligned and
`-1.00479e-7` jittered; no positive net heating is established at this
intermediate time. The aligned cooling bias and jittered residual need
the planned `t=6.5` comparison before a stability conclusion.

### Long-duration wall risk visible in the accepted exterior

No generator crosses the 1024-gon through the audited aligned
`t=3.300240722632501` snapshot. The closest site still has `0.102998`
radial-normal distance to the polygon side; the entire fluid's angular
momentum differs from `t=0` by only `-3.05e-11` relative. However the
exterior is evolving: for initial `r>=2.35`, mass-weighted absolute density
change rises from `0.0542` at `t=.825` to `0.1323` at `t=3.30`; among
sites currently at `r>2.35`, mass-weighted absolute radial speed rises
from `2.38e-5` to `1.71e-4`.
The nearest site's instantaneous radial speed is `1.94e-4` outward; its
distance divided by that speed is about `530` time units. This quotient
is **not** a collision-time prediction because the speed can reverse over
an orbit. It does show that the present no-impact abort cannot be assumed
safe through the paper's `t=600` image or longer endpoint without tracking
the wall margin, wall impulse/torque, and exterior transients. The current
closed-wall IC remains a controlled boundary experiment, not a direct
isolated-vacuum replay.

### Radial tail check on accepted native states

The [radial-tail reproducer](audit_kepler_bounded_radial_tails_20261027.py)
reads only SHA-verified snapshots whose 17-check native gate already passed.
It resolves signed/absolute entropy, native-minus-exact-orbit Voronoi density
residual, and orbital velocity error in six initial-radius bands. Its
whole-disk entropy and density means must reproduce the parent audit to
`1e-12` before it writes a receipt.

At matched `t=2.4758079598`, the inner `0.5<=r<0.75` band has mean absolute
density residual `2.047e-6` in the aligned run versus `6.946e-5` with 5%
site jitter (33.9 times). The corresponding circular-velocity-vector RMS
errors are `1.902e-7` and `2.026e-6` in units of local circular speed.
The jittered inner-band signed log-entropy change is `+2.50e-8`, near zero;
the aligned value is `-9.85e-6`. Thus the geometrically noisy run already
has larger density and orbit residuals without a measured positive mean
heating signal at this time. The accepted equal-time
[aligned](KEPLER_Q2_BOUNDED_RADIAL_TAILS_T247581_2026-10-27.json) and
[jittered](KEPLER_Q2_BOUNDED_JITTER005_RADIAL_TAILS_T247581_2026-10-27.json)
receipts preserve the annular and tail values.

The later aligned [t=3.30024 tail receipt](KEPLER_Q2_BOUNDED_RADIAL_TAILS_T330024_2026-10-27.json)
finds an inner-band mean absolute density residual `3.375e-6`, 99th
percentile `3.363e-5`, and maximum `4.447e-5`; its velocity-vector error
is `3.522e-7` RMS in circular-speed units. These values are small but rising
from `t=2.4758`; they rule out an already large hidden inner-cell orbit
failure at `t=3.3`. They do not establish late-time stability. Apply the
same diagnostic to both `t=6.5` endpoints before attributing any change to
shear heating or adding a spin/stress operator.

The jittered chain subsequently saved the same `t=3.3002407226326` as
the aligned chain. Its [17-check native audit](KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T330024_2026-10-27.json)
passes fixed mass, W=0, polygon area, positive state, zero injection, energy,
orbit and no-wall-crossing checks. Disk mean absolute density drift is
`1.48560e-3`; exact circular-site geometry predicts `1.47219e-3`.
The mean absolute native-minus-ballistic cellwise residual is `4.51163e-5`,
up from `2.25052e-5` at `t=2.4758`, versus `1.04511e-6` in the aligned
run at `t=3.30024`. Its [radial-tail audit](KEPLER_Q2_BOUNDED_JITTER005_RADIAL_TAILS_T330024_2026-10-27.json)
places the largest mean residual in the `0.5<=r<0.75` band:
`1.34349e-4` jittered versus `3.37538e-6` aligned (39.8 times).
The jittered band's 99th percentile is `4.86330e-4`, and its maximum
is `7.77121e-4`. Its signed mass-weighted log-entropy change is
`+1.96e-8`, effectively zero at this scale; the aligned band has
`-1.30634e-5`. The growing orbit/geometry residual currently precedes
any measured positive *mean* entropy production. This is evidence to
measure the force/mesh coupling at later times; it does not yet isolate a
specific force term or establish a long-time instability mechanism.

The [common-scale native density fields at `t=3.30024`](KEPLER_Q2_BOUNDED_ALIGNED_JITTER_DENSITY_T330024_2026-10-27.png)
show the accepted aligned and 5% site-jitter cases in the same paper-style
`[-2,2]^2` view. The upper row uses one absolute surface-density scale;
the lower row shows disk-relative density change on a shared signed scale.
The corresponding [plot reproducer](plot_kepler_bounded_aligned_jitter_density_20261027.py)
requires both SHA-verified 17-check receipts and equal native times. The
visible fine texture in the jittered lower panel is consistent with its
larger measured density drift; this picture alone does not identify a
heating mechanism or predict `t=600` morphology.
The [equal-time comparison receipt](KEPLER_Q2_BOUNDED_ALIGNED_JITTER_COMPARISON_T330024_2026-10-27.json)
quantifies that distinction: at `t=3.30024`, jitter raises total mean
absolute disk density drift 9.78 times, while exact-orbit geometry alone
predicts a 9.69-fold ratio. The native-minus-ballistic residual is 43.2
times larger in the jittered case, but accounts for only 3.04% of that
case's total mean absolute drift (0.688% in the aligned case). The
[comparison reproducer](compare_kepler_bounded_aligned_jitter_20261027.py)
requires accepted physics/tail receipts, exact snapshot hashes and common
physical time; it can be reused at `t=6.5`.

The aligned chain's next saved accepted state is exact
`t=4.126612443667704`. Its [native 17-check receipt](KEPLER_Q2_BOUNDED_NATIVE_T412661_2026-10-27.json)
passes, including fixed mass, W=0, positive EOS, area tiling, zero
injection and no wall crossing. Whole-disk mean absolute density drift is
`1.60583e-4`, almost the exact-orbit geometry value `1.60546e-4`;
native-minus-ballistic residual is `1.47152e-6` and relative total-energy
error is `3.35e-9`. The [radial-tail receipt](KEPLER_Q2_BOUNDED_RADIAL_TAILS_T412661_2026-10-27.json)
finds inner `0.5<=r<0.75` mean absolute density residual `4.78350e-6`
(up from `3.37538e-6` at `t=3.30024`), 99th percentile `5.15453e-5`,
and maximum `7.00418e-5`. The same band's signed log-entropy change is
`-1.62598e-5` (earlier `-1.30634e-5`). The trend remains small and toward
cooling in this aligned branch; a `t>6` checkpoint is still required.

The [detailed wall-margin audit](audit_kepler_bounded_wall_margin_20261027.py)
now records ten nearest sites, wall-normal and radial speeds, fixed mass,
W=0, and polygon containment from native snapshots. Its
[series receipt through `t=6.5`](KEPLER_Q2_BOUNDED_WALL_MARGIN_SERIES_TO65_2026-10-27.json)
includes the accepted `.825` state. For the two earlier SHA-verified
aligned states, this calculation gives minimum generator-to-wall clearance
`0.1029979874` at `t=3.30024` and `0.1027645983` at `t=4.12661`.
The same exterior generator (ID 55293) is closest at both times; all 72
near-wall sites remain more than 0.1 inside. Its instantaneous velocity
projected on the nearest polygon-side normal changes sign from
`-1.685e-3` to `+1.188e-3` across the two checkpoints. This includes
the tangential Kepler motion across polygon side normals. Neither an
instantaneous clearance/normal-speed quotient nor a two-point linear
extrapolation is a reliable wall-impact time. The slowly decreasing
clearance nevertheless keeps the missing reflection and wall impulse/
torque ledger relevant for a `t=600` claim.

### Exact-orbit geometry at the paper's longest durations

The [bounded IC long-orbit reproducer](audit_kepler_bounded_ballistic_long_20261027.py)
starts from the SHA-verified accepted `t=0` native states for aligned and
5% site-jitter layouts, fixes each cell's initial mass, and rotates every
generator on the exact Plummer circular orbit. It rebuilds ordinary W=0
Voronoi areas without evolving gas. All evaluated disk cells pass the
omitted-periodic-image safety test by a factor greater than 41. The
[machine-readable receipt](KEPLER_Q2_BOUNDED_BALLISTIC_T600_T1000ORB_2026-10-27.json)
reports disk mass-weighted mean `|V_i(0)/V_i(t)-1|`:

| Exact-orbit time | Aligned | 5% site jitter |
|---:|---:|---:|
| `6.5` | `1.85727e-4` | `1.47879e-3` |
| `600` | `1.96990e-4` | `1.46597e-3` |
| `2221.441469`, nominal unsoftened 1,000 periods at `r=.5` | `2.12966e-4` | `1.48904e-3` |
| `2222.361499`, actual softened 1,000 periods at `r=.5` | `2.10264e-4` | `1.47702e-3` |

The exact-orbit fixed-mass Voronoi density error does not grow to a large
mean value even at these long comparison times. If native gas evolution
later breaks, this static-IC geometry control by itself cannot explain a
large disruption. It cannot prove that the native pressure operator,
thermal work, boundary response or long-time integrator will remain
stable. The 1,000-orbit reference radius in the paper remains ambiguous.

### First accepted native bounded `t>6` result (2026-10-27)

A10 **437317** reached exact `t=6.5` and completed `0:0` in 1h20m48s.
The native endpoint `khout.002927.dat` has SHA-256
`6639d430782aa8d198991ecfac597a0caaae085155cfa452e16f14b2c6aab857`.
The [copied and rehashed 17-check physics receipt](KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json)
passes fixed mass, W=0, positive EOS, polygon area, no injection,
energy, disk density/orbit and no wall crossing. Total-energy relative
error is `5.263e-9`; disk mass-weighted mean absolute density change is
`1.86493e-4`, versus exact-orbit Voronoi `1.85727e-4`. Their cellwise
mean absolute residual is `2.88254e-6`, 1.55% of the native absolute
density change. Disk mean absolute log-entropy change is `7.47903e-6`.

The [radial-tail receipt](KEPLER_Q2_BOUNDED_RADIAL_TAILS_T65_2026-10-27.json)
gives a disk signed mass-weighted log-entropy change `-4.93911e-6`,
so this state still does not show positive *mean* irreversible heating.
For the inner `0.5<=r<0.75` band, the mean absolute native-minus-ballistic
density residual rose from `4.78350e-6` at `t=4.12661` to `1.15159e-5`;
its 99th percentile is `1.39557e-4`, and its signed log-entropy change
is `-2.53381e-5`. These small but growing residuals require later-time
control before claiming long-term stability.

The [wall-margin series](KEPLER_Q2_BOUNDED_WALL_MARGIN_SERIES_TO65_2026-10-27.json)
finds minimum generator clearance `0.1013591` at `t=6.5`, down from
`0.1027646` at `t=4.12661`; the same outer generator ID 55293 is
nearest. Its instantaneous radial speed is `+7.805e-4`, but neither
that speed nor its wall-normal component predicts an impact time.
Because the present fixed-polygon code still aborts on a crossing and
has no validated wall impulse/torque ledger, this accepted short-duration
result cannot be extrapolated to the paper's `t=600` image.

The already submitted unchanged-binary A10 follow-on **437328** started
automatically after the `t=6.5` native audit passed and targets exact
`t=10`. The separately guarded faster-binary continuation **437343**
remains dependent on 437328 and targets `t=20` only if the `t=10` gate
passes. Neither later endpoint is a result yet.
