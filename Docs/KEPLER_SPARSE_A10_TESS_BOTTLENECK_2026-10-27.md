# Native sparse Kepler A10 tessellation bottleneck (2026-10-27)

Jobs **437253** (wide, 55,931 sites) and **437254** (narrow, 54,312 sites)
were accepted fixed-mass W=0, RK4, CUDA+MPI+OpenMP evolutions before both
stopped near `t=2.51` at the 48-corner GPU clipping capacity. No
running binary or physical setting was changed for this diagnosis.

The native `[GPU tess]` and `[GPU step]` records near `t≈2.43–2.45` give
the following median of the latest 32 tessellation calls and latest eight
force calls, respectively. Each RK4 step makes four of each call. The logs
were read from only their last 256 KiB, without scanning the GPFS tree.

| Run | Real sites | Median periodic padding | Tessellation/call | Force/call |
|---|---:|---:|---:|---:|
| Wide | 55,931 | 55,680 | 455.571 ms | 54.1 ms |
| Narrow | 54,312 | 53,016 | 439.652 ms | 49.9 ms |

These medians place four tessellation calls near 1.82 s (wide) or 1.76 s
(narrow) of a contemporary roughly 2.2 s accepted step. The force call
includes SoA handling and upload as reported by `[GPU step]`; the tess
call includes CSR, packing, upload, device tessellation, download, and
writeback. Thus **the measured target for acceleration is tessellation**.
The present logs do not distinguish its GPU kernel from host preparation;
the ratio is not a prediction of speedup from any particular change.

The [preserved CUDA source](KEPLER_NATIVE_SOURCE_OVERLAY_2026-10-27/Exam/exam_gpu.cu)
has opt-in `GFS_GPU_TESS_PROFILE=1` device event timing and the
[host path](KEPLER_NATIVE_SOURCE_OVERLAY_2026-10-27/Exam/exam.c) prints
CSR/pack/upload/GPU/download/writeback timing for the first eight calls
and every 1,024th call. An isolated, low-priority
[five-step A10 profiler](kepler_sparse_a10_tess_profile_20261027.slurm)
was staged at
`/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_profile_wide_area_t001/to0p825/profile_t0835_a10`.
Its script and dependency hashes passed, and its diff from the accepted
wide A10 runner changes only resource/time fields, exact stop `.835`,
`GFS_GPU_TESS_PROFILE=1`, and the corresponding endpoint audit.
Slurm job **437271** completed (`0:0`) on another A10 node and reached
exact `t=.835`; its [physical receipt](KEPLER_SPARSE_A10_TESS_PROFILE_T0835_2026-10-27.json)
passes every gate. In its first eight device calls, Voronoi kernel A plus
the frozen-density copy takes 57.7–70.1 ms. Pressure, prefix scan, host
roundtrip and face scatter are each below 0.7 ms; host CSR and packing
together take about 16–22 ms. The later `t≈2.45` tessellation cost is much
larger at a similar periodic-padding count. Evolved geometry and candidate
search are stronger leads than a static padding-count explanation, but
these observations come from different physical times and nodes.
The prospective A10 and H200 long runners also enable this sparse
profiling output. They have not been staged or submitted with the updated
scripts, so the current native jobs remain untouched.

The [offline local-halo reconstruction](KEPLER_Q2_SPARSE_LOCAL_HALO_FALLBACK_T2221_RECONSTRUCTION_2026-10-27.json)
uses 9,091 rather than 53,074 periodic images at exact-orbit `t=2.221`,
but has no native GPU force/work validation. Keep global halo 32 until a
new implementation proves face/volume equivalence and accepted-step
wall-time improvement on the same physical state.

At the accepted `t≈2.475` states, independent periodic Voronoi polygons
have only 10–11 corners at most. The capacity-48 abort is therefore an
intermediate clipping-construction problem. A separate support-coalescence
branch has crossed the wide failure point; see the [paired abort
note](KEPLER_SPARSE_T251_GEOMETRY_ABORT_2026-10-27.md).

## Evolved-state GPU profile from the clipcompact continuation

The independently staged A10 continuation **437279** restarted the accepted
wide `t=2.55` checkpoint with the support-coalescence fix and
`GFS_GPU_TESS_PROFILE=1`. Its first eight device profiles report
`density_copy_A=461.471–487.412 ms` (median `465.523 ms`) per tessellation.
The corresponding complete host/device tessellation calls take
`484.140–508.401 ms`, while host CSR and packing together take approximately
`15–18 ms` per call. Kernel B, scan, host roundtrip and scatter are each
below `0.7 ms`. These are actual evolved-state measurements with the
corrected native executable, not the earlier `t=.835` control extrapolated
to a later geometry. The raw evidence is the first eight
`[GPU tess device profile]` and `[GPU tess profile]` lines in
`/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_profile_wide_area_t001/to0p825/clipcompact_to6p5_a10/log`.

At the later live `t≈2.717` checkpoint, each of four RK4 tessellation calls
took about `0.58–0.64 s`, and the printed full step about `3.3–3.4 s`.
The log reports roughly `1,980` cells requiring adaptive stencil expansion
per call, versus `55,931` real sites; it does not yet identify which cells
or how much of kernel A time those expansions consume. The source-level
next diagnostic is per-cell candidate/clip/box-retry work distribution or a
matched evolved-state A/B optimization. It would be premature to attribute
the 8-fold kernel growth relative to `t=.835` to periodic padding count,
clipcompact itself, or a single exterior region. Job 437279 was still
RUNNING when this profile was recorded; none of these timings is a `t>6`
physics result.

## Same-checkpoint GPU work distribution at t=2.55

The accepted clipcompact checkpoint at `t=2.55` was restarted for exactly
one step to `t=2.551` in two separate A10 jobs: baseline **437289** and
instrumented **437290**. Both completed `0:0` in eight seconds. Their staged
checkpoint SHA-256 was the same
`503dad7e489fdbf86b8869d3297ddb343675d6d200e73423a59d95e8420fea46`;
the runner and physical settings were identical. The accepted baseline binary
was `886638100a4ae057ccd6c3624c1b712141dea0656ccaa791934d0f00d44834fc`.
The [isolated work-counter patch](KEPLER_SPARSE_TESS_WORK_PROFILE_2026-10-27.patch)
reconstructs the probe's CUDA source SHA-256
`ff69bfdc19adfe0f49a5a72b04e2bd1af8c9a2474b288139afe49137e3e3964f`
from the accepted clipcompact source. Its probe binary SHA-256 was
`6bf4af2b5fb09d524aca6e2a5b90625e1c7ac278d58a3a92b99e0b889e65aebc`.

Both [baseline](/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_profile_wide_area_t001/to0p825/tessbase_t2551_a10/audit_physics_t2551.json)
and [instrumented](/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_profile_wide_area_t001/to0p825/tesswork_t2551_a10/audit_physics_t2551.json)
one-step endpoints pass all 15 physical and eight state/EOS checks. The
[ID-aligned paired audit](KEPLER_SPARSE_TESSWORK_PAIR_T2551_2026-10-27.json)
finds exact equality in every one of its 16 parsed particle fields and in
the snapshot metadata. The raw files differ in five *unparsed* bytes at
record offsets 1–5 for 51,074 cells; their hashes therefore differ. This
is a parsed-state equality result, not a claim of whole-record byte equality.

Each of the four instrumented RK4 tessellation calls traversed about
`21.1 million` neighboring entries in the eight-direction nearest pass and
another `21.1 million` in the seed pass, then about `1.79 million` in the
adaptive pass. About `1.34 million` half-plane clipping calls were made.
At the first call, `2,000` of `55,931` real cells expanded past the seed
stencil; `146` of them scanned more than 1,000 adaptive candidates.
One initially exterior cell, **ID 55899**, scanned `91,538` adaptive
candidates and invoked `52,825` clips, about `5.1%` and `3.9%` of the
respective whole-call totals. At `t=2.55` it is at `r=3.4576`, only
`0.168` from the periodic seam, and has a `0.1123` area, versus a typical
disk area near `0.00024`. The worst nearest/seed pass count was `420` at
disk ID 46912. The actual full GPU work lines are in the instrumented
job's `log`; the [pair reproducer](audit_kepler_tesswork_pair_20261027.py)
records the extreme site positions.

This locates a severe candidate/clip straggler in the periodic exterior and
also reveals two full seed-neighborhood scans for every cell. The counters
do not assign wall time to individual cells, so they do not prove the
straggler alone caused the `0.46–0.50 s` kernel time. The new bounded IC in
[the boundary audit](KEPLER_ISOLATED_BOUNDARY_GAP_2026-10-27.md) omits
ID 55899 because its initial radius is greater than `2.65`; that branch
must be profiled directly before claiming a speedup. The periodic sparse
continuation subsequently failed at `t=2.75514` from a distinct halo-width
limit, as recorded in [the abort note](KEPLER_SPARSE_T251_GEOMETRY_ABORT_2026-10-27.md).
