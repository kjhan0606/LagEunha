# Corrected sparse Kepler t=10 A10 continuations (2026-10-27)

The accepted H200 wide IC has 55,931 sites, a complete fixed-mass W=0
checkpoint at exact `t=.825`, and passing early physics/state audits. Its
H200 `t=10` continuation **437175** remained pending because of the H200
group GPU quota while full-q2 **437076** continued to `t=2.221`.

The accepted wide native executable contains an `sm_86` CUDA image, confirmed
with `cuobjdump -lelf`. A10 nodes had unallocated GPU cards. An isolated
[A10 runner](kepler_q2_sparse_area_wide_to10_a10_20261027.slurm) and
[stager](stage_kepler_sparse_wide_a10_to10_20261027.py) were created under
`/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_profile_wide_area_t001/to0p825/to10_a10`.
The A10 runner differs from the accepted H200 wide runner only in job name,
partition, `--gres=gpu:A10:1`, time limit, and a device log. It uses the
same 55,931-site accepted binary, physics parameters, timestep setting,
halo 32, RK4, one MPI rank, two OpenMP threads, fixed mass, W=0, and
the exact `.825` checkpoint. Its staged script SHA manifest, shell syntax,
and Python imports passed. The restart selector rechecked the parent audit,
source checkpoint, binary and parameters before copying native inputs.

A10 job **437253** started on `syn02`; `device.log` identifies an NVIDIA
A10 and the native log reached `t=.869616`, step 449, with zero halo clips,
positive internal energy, zero floor/dual/entropy injection, and a total
energy drift of `3.825e-9` relative at that point. This is an observed
early interval, not a `t=10` result. The output and required exact-time
audit remain pending.

To prevent duplicate wide evolution, H200 wide job **437175** is held
(`JobHeldUser`) while A10 **437253** runs. If A10 fails, release 437175
after diagnosing its final state; if A10 passes the exact `t=10` physics
and state audits, the pending H200 duplicate can be canceled without
executing. No running H200 job was changed.

The independent 54,312-site narrow IC now also has an isolated
[A10 runner](kepler_q2_sparse_area_narrow_to10_a10_20261027.slurm) and
[stager](stage_kepler_sparse_narrow_a10_to10_20261027.py) in
`/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_profile_buffered_area_t001/to0p825/to10_a10`.
Its runner differs from the accepted H200 narrow runner only in the same
resource/time/device-log fields as wide; its physics and exact `.825`
restart are unchanged. Staged script hashes, shell syntax, and Python
imports passed. A10 job **437254** started on `syn02` concurrently with
wide **437253**. The two `device.log` files identify distinct NVIDIA A10
GPU UUIDs. Pending H200 narrow **437174** is held (`JobHeldUser`) while
437254 runs, with the same release-on-failure and cancel-on-validated-success
policy. Neither A10 run has a `t=10` result yet.

## First common saved checkpoint

Both native jobs wrote complete step-852 checkpoints at actual times
`1.6502423888` (narrow) and `1.6502423914` (wide). Independent CPU
Slurm jobs **437257** and **437256** completed the
[narrow radial audit](KEPLER_SPARSE_NARROW_A10_T1650_RADIAL_2026-10-27.json)
and [wide radial audit](KEPLER_SPARSE_WIDE_A10_T1650_RADIAL_2026-10-27.json).
All five finite/fixed-mass/W0/positive-state/`rho=m/V` checks pass.

| Mass-weighted initial `0.5<=r<2` disk statistic | Narrow | Wide |
|---|---:|---:|
| Absolute relative density change | `1.24963e-4` | `1.24537e-4` |
| Absolute log-entropy change | `3.90231e-6` | `4.06547e-6` |
| Signed log-entropy change | `-1.69342e-6` | `-1.93383e-6` |
| Absolute orbit-position error | `1.76853e-8` | `1.05130e-8` |

No generator crossed the inner `r=.5` disk edge in either case. Their
total-energy changes are `2.32735e-8` and `2.32878e-8` relative.
Current-edge cells within 0.2 of a box face already have about `0.2653`
absolute density change in both. These results show an early outer
periodic-seam disturbance while the original disk remains close to its
own IC. They do not establish that the seam has or has not caused the
later `t>6` disk behavior, nor do they prove net irreversible heating.

The [matched time-series auditor](audit_kepler_sparse_a10_series_20261027.py)
was checked on the two saved times in CPU job **437258** (`COMPLETED 0:0`,
one second). Its [partial receipt](KEPLER_SPARSE_A10_SERIES_PARTIAL_T1650_2026-10-27.json)
has `target_reached=false` as required before `t=10`. From `.825` to
`1.650242`, the disk mean signed log-entropy increments are
`-1.063e-6` (narrow) and `-1.273e-6` (wide); mean absolute increments
are `2.474e-6` and `2.567e-6`. The 49,038 exactly matched original disk
sites differ in relative density change by `1.50e-7` at `.825` and
`6.35e-7` at `1.650242` in mass-weighted absolute norm. Their entropy-log
differences are `8.10e-8` and `3.23e-7`. Current-edge mean absolute
density change grows from `0.1303` to `0.2653`, while disk density remains
near `1.25e-4`. This is a measured spatial trend, not proof of a causal
edge-to-disk signal. The same auditor will inspect all saved points through
`t=10` after both native jobs finish.

## Audits queued after the native endpoints

The [CPU consistency runner](kepler_sparse_a10_t10_consistency_20261027.slurm)
has SHA-256 `2b047c35327ff9974454bc0f41f2932464cec9c0bc0f5485a15bf3860dd7d7ff`.
It requires the exact `t=10` physics receipt, all of its checks, the expected
site count, and an endpoint inside the corresponding run directory before
running the independent state/EOS audit. The wide state job **437259** has
`afterok:437253`; the narrow state job **437260** has `afterok:437254`.
Each writes `audit_consistency_t10.json` in its own run directory, refusing
to overwrite an existing receipt.

CPU time-series job **437261** has `afterok:437259:437260` and will compare
all complete saved times through exact `t=10`. Its intended receipt is
[the t=10 matched series](KEPLER_SPARSE_A10_SERIES_T10_2026-10-27.json).
These jobs are dependencies, not completed `t=10` results. A physics or
state-audit failure keeps the comparison from running and must be diagnosed
before any long continuation is accepted.

## Geometry separation at the historical failure window

The [t≈6.6 posthoc runner](kepler_sparse_a10_t6p6_posthoc_20261027.slurm)
has SHA-256 `5f2ed9f7464232baafb1cbd5cd513548866e914dfb9a306261104ab38bc2eb1a`.
Jobs **437262** (wide, `afterany:437253`) and **437263** (narrow,
`afterany:437254`) select the unique complete saved checkpoint within
`|t-6.6|<.02`. They independently audit radial state and compare each
native disk density field with the exact-orbit Voronoi prediction. An
`afterany` dependency preserves this diagnostic if the native run fails
*after* saving the historical failure-window snapshot; a failure before
the snapshot produces an explicit missing-checkpoint error.

The [exact t=10 geometry runner](kepler_sparse_a10_t10_geometry_20261027.slurm)
has SHA-256 `ad193f989d6083bbe2a28a72408ed0525823c53a71334330ed73c7e5806efcfd`.
Jobs **437264** (wide, `afterok:437259`) and **437265** (narrow,
`afterok:437260`) require matching passing physical and state receipts
before writing `audit_t10_native_ballistic.json` in their respective run
directories. These two time points will show whether early native density
agreement with exact-orbit tessellation survives the former `t>6` wall.
Neither geometric comparison by itself measures irreversible heating.

The [t=10 seam-provenance runner](kepler_sparse_a10_t10_seam_20261027.slurm)
has SHA-256 `20ae1a387cdd5ed7a3e4799e493c447d04dd8ded4cb92d9e7e7defc645be8ed4`.
Jobs **437269** (wide, `afterok:437259`) and **437270** (narrow,
`afterok:437260`) require the same physics/state lineage and write
`audit_t10_seam.json`. This reports where initially exterior fixed-mass
generators have moved and how their density and entropy evolved. A crossing
is a generator-centre proxy; even zero crossings cannot exclude an inward
pressure or torque signal through faces. Inspect these receipts alongside
the radial and geometric audits before starting an expensive continuation.

## Both A10 runs stopped at a GPU clipping capacity (2026-10-27)

Wide **437253** ended `FAILED 1:0` after the accepted step-1329 start
state `t=2.5176961142503815`; the next RK4 geometry stage reported
`[GPU_CLIP_CAPACITY] nc=48 capacity=48`, then
`[GPU_TESS_INVALID] volume=1 geometry=1 halo=0`. Its complete accepted
regular checkpoint is `khout.001302.dat` at `t=2.4761277758878566`.
The extra `khdiag.step001329.start_geometry.dat` is structurally complete
at the accepted step start, but has not been promoted as a restart.
Narrow **437254** ended `FAILED 1:0` at the end-of-step final-volume
geometry after the complete regular `khout.001318.dat` at
`t=2.4753977144699952`; its log has the same capacity-48 diagnostic,
`[GPU_FINAL_VOLUME] invalid=1 geometry=1 halo_clipped=0`. No exact
`t=10` result exists. The original H200 t=10 duplicates **437175** and
**437174** were subsequently canceled; their unchanged binary would likely
repeat the capacity failure. The now impossible `afterok` audits 437259–437261,
437264–437265 and 437269–437270 were canceled. The `afterany` t≈6.6
audits 437262–437263 ran and failed because the parent jobs stopped before
that checkpoint.

Independent CPU radial jobs **437273–437274** passed every state check on
the two last complete regular checkpoints, with fixed mass, W=0, positive
volume/internal energy and `rho=m/V`. The independent periodic Qhull
corner-count jobs **437275–437276** used 2-unit image halos and reproduced
the native areas to maximum relative `5.86e-8` (wide) and `5.73e-8`
(narrow), with global area 36. Their maximum *final* cell corner counts
are only 10 and 11; neither has a final polygon with 32 or more corners.
Thus the 48-corner overflow is an **intermediate GPU clipping construction
problem**, not evidence that a physical final cell needs 48 faces at these
accepted snapshots. The failing RK4 stage itself still needs diagnosis.
Receipts: [wide corners](KEPLER_SPARSE_WIDE_A10_T2476_CORNERS_2026-10-27.json),
[narrow corners](KEPLER_SPARSE_NARROW_A10_T2475_CORNERS_2026-10-27.json),
[wide radial](KEPLER_SPARSE_WIDE_A10_T2476_RADIAL_2026-10-27.json),
[narrow radial](KEPLER_SPARSE_NARROW_A10_T2475_RADIAL_2026-10-27.json).

The separate five-step A10 profiling job **437271** completed at exact
`t=.835`, with all native physical gates passing. Its first eight device
profiles put Voronoi kernel A at about 58–70 ms per call, versus <1 ms
each for pressure, face scan/roundtrip and face scatter; host packing and
CSR cost about 16–22 ms. The later running jobs reached 440–460 ms
tessellation calls near `t=2.45`, and then the clip-capacity stop.
Do not interpret the early profiling duration as a late-time speedup.

An independently staged [support-coalescence branch](KEPLER_SPARSE_T251_GEOMETRY_ABORT_2026-10-27.md)
has since crossed the wide failure point in A10 job **437277**, reaching
exact `t=2.55` with its native state, radial, exact-orbit geometry and seam
receipts passing. It changes the clipping construction while retaining
the 48-corner capacity and physical settings. This short success requires
a longer resumed trajectory before the original `t>6` question is answered.
Independent CPU job **437278** then tessellated all 55,931 corrected wide
native generator positions at exact `t=2.55` with 55,691 periodic images.
Its [Qhull receipt](KEPLER_SPARSE_WIDE_CLIPCOMPACT_T255_CORNERS_2026-10-27.json)
finds maximum final polygon size 9 and maximum relative native-vs-independent
cell-area difference `5.84e-8`, with total area 36. This checks the new
clipping construction against the full native state after the former abort;
the `t=10` and `t>6` gates remain open.
