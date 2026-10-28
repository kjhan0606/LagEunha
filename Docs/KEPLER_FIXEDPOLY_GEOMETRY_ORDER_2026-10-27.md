# Fixed-polygon CUDA clipping order and accepted-state replay

The native 1,024-gon W=0 Kepler branch passed its first `t=.001` step and
the [wall-face audit](KEPLER_FIXEDPOLY_WALL_DUMP_2026-10-27.md), but the
first attempted continuation to `.825` stopped at accepted step start
`t=.045616533734626756`. The original `GPU_CLIP_CAPACITY` report had a
**real generator** support (`neighbor_soa=30966`) at 48 corners; its
accepted-state snapshot is
`/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/to0825_v2/khdiag.step000024.start_geometry.dat`,
SHA-256 `08a7500e5d72c965c771ad51fd90d87a3cd7872a439edfb146df7a566d18c76a`.
It is a complete accepted restart, not the failed RK4 stage.

Moving all 1,024 wall cuts before the real candidate stream moved the
failure rather than removing it. A replay of the same checkpoint using
binary SHA-256
`2c4ac265d9bb9e181bc4af2bb5412af65eefb9a45cc0a3ad6c38dfe5bef13327`
failed at the initial geometry calculation: several **wall** supports
`-1000-k` filled the 48-corner buffer. The preserved A10 run is
`/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/capacity_replay005_clipfirst_wallprobe_v1`.

The [two-stage patch](KEPLER_FIXEDPOLY_COARSE_FINAL_CLIP_2026-10-27.patch)
uses eight *actual* polygon halfplanes, every 128th side, immediately
after the first nearest-real-generator cuts. Their intersection is a
superset of the exact 1,024-gon because these eight constraints are a
subset of the 1,024 constraints. Real candidates then clip this smaller
superset. Finally **all** 1,024 wall halfplanes are applied. Thus this
changes clipping order and intermediate polygon size, while retaining
the intended exact final fixed polygon and conservative candidate search.

Source SHA-256 `c40fe698f5f0be4323e62aff7312177486440f4f1136d39416276536bde30e56`
and binary SHA-256
`8a038be1ab7b6b13bc55b3560fa2a88c607e195fed0b8d6a006b6a5f3b745872`
were used in A10 replay **437314**. The physics executable exited `0`
after three accepted RK4 steps to exact `t=.05`. The Slurm script then
exited `1` because the first version of its Python auditor attempted to
serialize a NumPy boolean. The auditor was fixed and rerun on the saved
output, without repeating the simulation. Its
[passed receipt](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/capacity_replay005_coarsefinal_v2/audit_t005_post.json)
finds fixed mass, W=0, positive state, bounded area tiling, zero injection,
energy change `3.14e-12` relative to the accepted restart, and maximum
incremental disk-orbit position error `1.26e-11`.

A concurrent refinement added 24 more nearest-real preclips before the
coarse wall cut. The current source SHA-256 at that stage was
`2194c4aa81ea0eee030c5b07955794d0b7dc8d25dbc32f5049d14c71add67865`
and binary SHA-256
`3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e`.
The exact seven-file [native source snapshot](KEPLER_FIXEDPOLY_NEAR32_SOURCE_2026-10-27/README.md)
is preserved outside `syntax`'s temporary filesystem.
Its independently staged A10 continuation **437315** had reached
`t=.226022` without a geometry abort by 07:11 KST; the accepted endpoint
and physical audit are recorded below.

The extra 24-nearest search currently runs for every site when fixed-poly
mode is enabled. If it becomes a significant cost, restrict it to cells
whose initial square can touch the fixed polygon and verify a paired
one-step state/geometry result. Do not trade away the exact wall or the
geometry-abort guard merely to reduce runtime. The long `t>6` stability
question remains open.

The same near32 binary passed a **fresh** A10 `t=.001` native and first-wall
preflight, job **437316**, both auditors exit `0`. The
[native receipt](/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_bounded_r275_wallpreflight_near32_a10_20261027/audit_t001.json)
retains fixed mass, W=0, positive state and energy to `3.62e-14` relative;
the [wall receipt](/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_bounded_r275_wallpreflight_near32_a10_20261027/audit_wall_call1.json)
again finds 1,096 wall fragments covering all 1,024 sides. A guarded
[t=3 continuation stager](stage_kepler_bounded_to3_20261027.py) and
[A10 runner](kepler_q2_bounded_poly1024_to3_a10_20261027.slurm) are prepared.
The stager refuses to run until the near32 `.825` receipt exists and every
physics check passes; no t=3 job has been submitted at this update.
The near32 native `t=0` snapshot and the earlier `t001_v2` `t=0`
snapshot have different raw hashes, but an ID-ordered comparison found
exact equality in all ten parsed physical fields checked: positions,
velocities, mass, volume, internal energy, density, pressure, and W=0.
This supports using the fresh near32 initial state as the reference for
the continuing `.825` chain; it does not assert raw byte equality.

### Native `.825` gate and next run

A10 job **437315** completed `0:0` at exact `t=.825` in 11m16s.
Its [physics receipt](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/to0825_near32_v1/audit_t0825.json)
passes all 17 checks: disk density absolute mass-weighted change
`1.1652e-4` versus exact-orbit Voronoi `1.1654e-4`, their pointwise
residual `8.47e-8`, disk entropy-log absolute change `2.52e-6`,
total-energy relative change `6.67e-10`, and zero floor/dual/entropy
injection. It uses the `3fc66a...` binary and endpoint SHA-256
`134e4e76a475c6923f1106ffb944fd0d51a0df82e7d4c86ea2814057ce9a12f6`.

A separate worker had already submitted A10 **437317** to exact `t=6.5`
from this same binary and endpoint, with checkpoints every `.825`. The
prepared `to3_near32_v1` directory and runner were **not submitted** to
avoid duplicate evolution. Job 437317 and its later disk/outer physical
audit, rather than the initial and `.825` gates, decide whether the
historical `t>6` wall has been crossed.

The [A10 runtime receipt](KEPLER_BOUNDED_A10_RUNTIME_T0825_2026-10-27.json)
measures 402 accepted steps and 1,608 GPU tessellation calls in
671.391 solver seconds. Tessellation took 264.034 seconds (39.3% of
wall time; median call 166.229 ms). The `CPU step` timers covered
404.7 seconds (60.3%); **their update interval includes the GPU
tessellation calls**, so these times must not be added. The remaining
266.691 seconds lie outside those CPU step timers and need finer
instrumentation. Optimizing tessellation alone has a limited maximum
speedup. Median accepted `dt` was `.00194`; extrapolating both cost and
`dt` unchanged to `t=600` gives about 143 A10 hours, **only a cost
indicator**, not a prediction of long-time runtime. The inner step/update
and unmeasured integration paths need profiling alongside GPU tessellation before
promoting a practical long-time method.
