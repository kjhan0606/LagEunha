# Kepler periodic-image fallback: evidence and native contract (2026-10-27)

## Why this is being considered

The one-rank H200 q2 continuations use the accepted CUDA+MPI+OpenMP W=0
RK4 solver. Their outer periodic tessellation can run out of candidate
images even when the disk is geometrically quiet. A wider global image halo
protects the current runs, but it adds many sites to every RK stage. This
document defines a possible **opt-in** local repair; it does not promote the
repair to the accepted solver or change the running jobs.

## Evidence at exact-orbit t=2.221

The [offline auditor](audit_kepler_local_periodic_fallback_20261027.py)
starts with an eight-cell image halo, identifies polygons whose candidate
search reaches the grid edge, and takes the union of otherwise missing
periodic images within twice each incomplete polygon's circumradius. It
then builds a new Voronoi diagram using this reduced image set and compares
the affected cell areas with an independently built 2.5-unit-margin diagram.
The margin exceeds twice the full circumradius for every audited risk cell.

| IC | Clipped real cells | Base images | Reduced images | Wide images | Total reduced sites / total wide sites | Maximum affected-cell area difference |
|---|---:|---:|---:|---:|---:|---:|
| Full q2, 147,456 real, [receipt](KEPLER_Q2_LOCAL_HALO_FALLBACK_T2221_RECONSTRUCTION_2026-10-27.json) | 256 | 44,968 | 118,535 | 205,191 at 28 cells | 265,991 / 352,647 | 5.22e-12 relative |
| Sparse narrow q2, 54,312 real, [receipt](KEPLER_Q2_SPARSE_LOCAL_HALO_FALLBACK_T2221_RECONSTRUCTION_2026-10-27.json) | 32 | 460 | 9,091 | 53,074 at 32 cells | 63,403 / 107,386 | 1.97e-13 relative |

Only 36 full-q2 and eight sparse new periodic images actually support
affected Voronoi faces, but an online algorithm does not know those supports
before reconstruction. The larger conservative candidate unions are the
relevant implementation counts. At this sampled time they reduce the total
site input by about 25% and 41%, respectively, relative to the specified
global halos. They are **not measured wall-time improvements**.

This is a W=0, exact-orbit, wrapped-position geometry result. It neither
describes valid periodic exterior gas evolution nor proves native stability
between samples or at t=120/600. The outer velocity/potential seam remains a
separate physical boundary problem.

## Measured cost of a modest global halo increase

The [matched native timing receipt](KEPLER_Q1_HALO8_12_MATCHED_COST_2026-10-27.json)
compares the same 147,456-site q1 CUDA binary and physical interval
`0.630053<=t<=0.649927`: 263 accepted steps in the original eight-cell
halo run and 263 in its 12-cell recovery. Mean padded sites rise from
`51,555` to `80,930` (1.570 times). The printed mean CPU step time rises
from `0.8` to `1.0` seconds (about 25%), the GPU wrapper's `SoA`-labeled
interval from `37.18` to `50.68` ms (36%), and logged GPU step total from `75.46` to `101.39` ms
(34%). Tessellation time itself rises only from `48.27` to `51.61` ms (7%).
The CPU timings are rounded to 0.1 s per line, and this narrow window does
not predict the cost of a 28-cell halo or a local fallback. It does show that
global padding already costs substantial host refresh and per-step time
before the longer Kepler run begins. Preserving full periodic correctness
without globally padding every stage is therefore a measured performance
need, not just an image-count preference.

## One-rank host path to profile before changing it

Source review of the active `LagEunha` worktree on `syntax` identifies a
specific host path. `RkS/mpirks.c::ppadding` marks candidate images with an
OpenMP loop, then fills their particle payloads serially. Its OpenMP branch
first allocates `100+0.1*N` particle slots and resizes to the counted exact
image number before filling. On one rank,
`adv_justppadding` performs no MPI exchanges, but still calls
`DeletePtlOutOfLocalBox`; `ppadding` then calls
`DeleteInsideLocalBox`, giving two further serial image scans.
`Exam/exam.c::mkLinkedList2D_oExam` bins real and padded particles serially,
while `Exam/exam_gpu_extract.c::fillParticleSoA` is OpenMP parallel. This
explains why linking with OpenMP does not mean every host stage uses all
allocated cores.

These are source-level candidates, **not measured bottleneck shares**. The
`[GPU step] SoA` interval is timed in `Exam/exam_gpu_extract.c` from the
start of its `tess_on_device` branch through the **whole host periodic-padding
rebuild, alias/order checks, SoA fill, and GPU re-upload**. It is not a
measurement of `fillParticleSoA` alone. The earlier 14-step
[local device-refresh A/B](KEPLER_TIME_STEP_SCALING_2026-10-26.md#independent-polar-local-sync-ab-queued)
made this interval fall from `36.25` to `0.1` ms with bitwise-equal saved
fields, but does not establish a sustained whole-step speedup or behavior
after the seam deforms. Independent evolved-state jobs 437180/437181 are
already queued for that test. The current logs still do not separately time
padding, the two scans, or linked-grid insertion. Profile those sections on
identical accepted restart states before changing their algorithms. An
exact-size allocation after image counting is a possible byte-preserving
micro-optimization, but its step-time benefit is unmeasured and likely much
smaller than eliminating unnecessary global images or increasing the
physically valid sparse-core timestep.

The active full-q2 and corrected sparse-narrow `params.dat` files are
byte-identical (SHA-256 `2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43`).
They specify `Kappa=0`, `entropy_mode=0`, `av_mode=5`,
`cd_amax=cd_amin=0`, and zero physical
viscosity, matching the *static* eligibility checks of the existing
`GFS_GPU_LOCAL_PERIODIC_SYNC=1` branch for one MPI rank. Its device kernel
still verifies owner IDs, periodic shifts, equal velocity/mass/energy and
zero weights at runtime; a failed alias check falls back to the ordinary
host refresh. The option is specific to this W=0 one-rank case and is not
yet promoted for the evolved halo-12 or sparse halo-32 continuations. The
low-priority matched jobs 437180/437181 will test the evolved branch before
any timing-based promotion.

## Live halo-28 cost sample (2026-10-27)

The unchanged full-q2 H200 continuation **437076** began from its audited
`.825` checkpoint with a 28-cell global periodic halo. At accepted
`t≈.835`, the log reports about `218,940` padded images and rounded
CPU-stage cost `1.9 s/step`; the GPU wrapper reports about `0.21 s` and
its tessellation calls about `0.107 s` each. The preceding halo-12 q2
endpoint had about `80,000` padded images, rounded CPU-stage cost
`1.1 s/step`, and GPU wrapper about `0.096 s`. Across a short live interval,
21 halo-28 accepted steps consumed 53 seconds of Slurm elapsed time, or
`2.52 s/step` including all runtime overhead. The halo-12 recovery's
entire `.691→.825` path averaged `3060.390/1785=1.714 s/step`.
These are different time intervals, so they quantify a practical cost
increase, not a controlled same-stage microbenchmark.

If both the observed `7.64202e-5` timestep and `2.52 s/step` were constant,
the remaining `.835→2.221` interval would take about **12.7 hours**, and
the subsequent `2.221→6.5` interval about **39.2 hours**. The 48-hour
second-job limit therefore has finite headroom; timestep shrinkage or more
expensive geometry may require a checkpointed shorter segment. These are
capacity projections only. The active job's physical gates and actual wall
time decide any promotion. The live data reinforce that global halo width
has a material host cost and that sparse IC and local-image work remain
relevant to an AREPO-scale evolution rate.

## Native implementation contract

1. Keep the accepted source and binary immutable. Build an isolated opt-in
   variant. Keep the accepted stage timestep, EOS, mass, and W=0 physics.
2. In `Exam/exam_gpu.cu`, distinguish **halo clip only** from invalid volume,
   velocity fit, or polygon geometry. The current
   `gpu_tessellate_and_build_faces()` aborts (or marks an RK4 geometry trial
   failed) after reading these counters. An opt-in recoverable status must
   return the clipped owner IDs and incomplete polygon circumradii without
   converting unrelated failures to halo retries.
3. On a one-rank W=0 first implementation, use the complete real-particle
   stage catalogue to generate unique `(owner ID, x shift, y shift)` images
   within `2 R_partial` of the clipped owners, plus a numerical margin.
   Preserve the full particle payload for each image: position, velocity,
   mass, thermal state, density, pressure, weight, and owner identity.
   Rebuild the linked grid, CellCSR, SoA and GPU buffers so every candidate
   has a valid index. A coordinate-only wrap inside the tessellation kernel
   is insufficient because gradients and face/flux code dereference its
   neighbor index and fields.
4. Re-upload the **original host RK4-stage state** and rerun the complete
   tessellation and face build. Kernel A overwrites device density and
   pressure; retrying it against its previous output would change the
   physical stage. Do not write any failed attempt back to host particles.
   Keep `nbp` as the number of real owners and count new images only in
   `n_total`. Verify face reciprocity, unique image mapping and conservative
   force/flux bookkeeping after retry.
5. Require an MPI-wide decision before retry or abort so all ranks take the
   same path. The initial one-rank variant must reject multi-rank use until
   remote owner payloads and image exchange are implemented. If a clipped
   cell cannot be certified from available candidates, use a controlled
   global-halo rebuild or fail explicitly; never accept a clipped polygon.
6. This `2 R_partial` selection is specific to W=0 Voronoi. For a weighted
   Laguerre extension, use a weight-aware search bound and re-audit it; the
   existing kernel includes `global_w2_max - w2i` in its search radius.

The host insertion point is `getAccVoro2DBlend_GPU()` in `Exam/exam.c`,
between its GPU tessellation call and result download/writeback. It already
owns the padded AoS, linked grid and packed SoA for a stage. The CSR builder
in `Exam/exam_gpu_extract.c` maps linked-list entries to real/pad indices;
any appended image must participate in that mapping. The detection point is
the existing `PolygonOut::halo_clipped` and device counter in
`Exam/exam_gpu.cu`.

## Promotion gates

- At a fixed accepted checkpoint, compare opt-in local fallback and a wide
  global halo at the same RK4 stage: affected-cell polygons, areas, face
  neighbors, reciprocal flux/force pairs and unchanged real-particle input.
- Continue both at the same timestep and compare ID-wise fixed mass,
  positivity/floor injection, density, signed and absolute entropy change,
  angular momentum, energy and the GIZMO disk annulus norms. Include an
  outer-seam diagnostic so boundary failure is not called disk heating.
- Measure wall time per RK4 stage, tessellation, padding/CSR, H2D/D2H and
  fallback frequency on H200. A site-count reduction alone is insufficient.
- Only after these gates, test multi-rank halo exchange and any W2 != 0
  generalization. The existing wide-halo H200 jobs remain the reference.
