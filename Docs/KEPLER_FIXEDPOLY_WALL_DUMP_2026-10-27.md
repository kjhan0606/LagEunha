# Native fixed-polygon wall-face diagnostic, 2026-10-27

The bounded Kepler CUDA force kernel now represents polygon edges with
negative supports `-1000-k` and adds their relative-sweep energy correction.
The existing `GFS_GPU_FACE_DUMP_CALL` diagnostic used to discard every
negative neighbor before writing `face_dump_r*_call*.tsv`. It therefore
could not show the wall faces whose pressure traction, torque, and energy
correction need checking. This was an observation gap, not evidence of a
force-kernel failure.

The [diagnostic-only patch](KEPLER_FIXEDPOLY_WALL_DUMP_2026-10-27.patch)
retains physical wall faces in that dump, labels them `wall=2`, and writes
the negative support as `id_j`. Set `GFS_GPU_FACE_DUMP_WALL_ONLY=1` for a
compact one-call dump rather than all interior faces. The patch changed
only the host-side diagnostic output in `Exam/exam_gpu.cu`; it did not change
the CUDA geometry, force, or energy arithmetic. It was applied to the
isolated source `/tmp/Kepler_sparse_fixedpoly_20261027` on `syntax` after
checking its original source SHA-256
`75c9377267e150b468c7745d5cd9a8d9e367081735d2e19d3ca473f6d41a3a9b`.
The patched source SHA-256 is
`7d6299b1ceadd3c1b1a6adce81ac91b02449ad9931341c151cf0c1a9dd1de31b`.
The patch also passed a dry-run against the captured source snapshot.

For the native **initial** force call, use
`GFS_GPU_FACE_DUMP_CALL=1 GFS_GPU_FACE_DUMP_WALL_ONLY=1 GFS_GPU_TESS_DUMP=0`.
The [wall-face auditor](audit_kepler_fixedpoly_wall_faces_20261027.py)
checks all 1,024 polygon sides, their total and individual lengths,
support-plane position, outward normal orientation, uniform initial
pressure, cancellation of total wall force and central torque, and
`work=-P*udotds` for every wall face. It reports the pressure traction
and the relative-sweep correction separately: a stationary wall does no
physical pressure work. The auditor compiled and passed a synthetic
1,024-side polygon capture, with relative side-length error below
`1.6e-13`; it has **not yet been run on native CUDA output**.

The diagnostic is only an initial-wall gate. Later calls can have
nonuniform pressure and nonzero wall force/torque; their physical
interpretation needs a time-integrated wall ledger and the disk/outer
radial audit. A passed initial face dump does not prove long-time Kepler
stability or a match to the GIZMO paper.

For a later force call, run the same auditor with `--allow-nonuniform`.
It still gates exact wall coverage, normals, and the per-face energy
formula, but reports the nonzero pressure force and central torque
without treating them as an initial-uniform-pressure failure. Both modes
passed separate synthetic 1,024-side captures; the nonuniform capture
produced a deliberately nonzero wall force. A single later force call is
an instantaneous diagnostic, not a time-integrated wall ledger.

## First native A10 result

The immutable [A10 retry run](/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_bounded_r275_wallpreflight_a10_retry1_20261027/staging.json)
used binary SHA-256
`bff5384b0a57e90beb8b736ef0285c7fc6622eb674cb43de0c4f824e76cef369`.
Slurm **437309** completed with exit `0:0` in seven seconds. Both the
[native `t=.001` receipt](/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_bounded_r275_wallpreflight_a10_retry1_20261027/audit_t001.json)
and [first-force wall receipt](/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_bounded_r275_wallpreflight_a10_retry1_20261027/audit_wall_call1.json)
passed every stated check. The earlier A10 job **437305** exited before
physics because the stager copied the binary without execute permission;
its log is preserved. The stager now preserves source modes. The duplicate
pending H200 job **437303** was canceled after the A10 pass.

The native first call reported 72 boundary cells and 1,096 wall-face
fragments covering every one of 1,024 polygon sides. Their total length
`17.278732489024375` matches the exact perimeter
`17.278732489024378`; the worst individual side-length relative error
is `1.02e-11`. The maximum face-to-wall plane residual is `1.34e-15`.
Uniform initial pressure gives fluid wall force
`(3.98e-14, 7.80e-21)` and central torque `1.71e-21`. The summed
relative-sweep energy correction nearly cancels, while its absolute sum
is `3.11e-7`, so omitting it on individual boundary cells would be wrong.
The native first step retained fixed mass and W=0, positive state,
area tiling, zero floor injection, disk circular positions to
`2.52e-15`, and relative total-energy change `3.62e-14`.

This **does not promote the branch to a long run**. The independently
staged A10 evolution **437302** stopped at accepted step start
`t=0.0456165337` after a `GPU_CLIP_CAPACITY` error in a subsequent RK4
stage. Its failure is geometric and precedes the disk's historical
`t>6` stability test. A concurrent source change moved fixed-polygon
clipping before the remaining candidate stream to address the transient
48-corner overflow; that change needs its own native regression and
performance gate before any longer physical inference.

## Instantaneous wall state at accepted `t=.825`

The near32 branch subsequently passed its `.825` physical gate. An
independent one-step A10 probe **437318**, starting from that exact
checkpoint, dumped its first force call at `t=.825` and stopped at
`t=.826`. The [wall receipt](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wall0825_probe_near32_a10_v1/audit_wall_call1.json)
passes all geometry and per-face energy-formula checks: the same 72 cells
and 1,096 fragments cover every side with maximum length error
`1.01e-11` relative. Wall pressure now varies by up to `8.07%` relative
to the first boundary face, so the appropriate nonuniform-pressure audit
reports rather than rejects the force. The instantaneous fluid wall
traction is `(1.70e-10, 2.27e-10)` and central torque `5.29e-13`.
These are small at this early time, but a single force call does not
bound their time integral or later influence on the disk. The signed
cell-energy correction is `-2.01e-9`, compared with absolute correction
`3.12e-7`; neither is physical work by the stationary wall.

The paired [one-step receipt](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wall0825_probe_near32_a10_v1/audit_t0826.json)
passes fixed mass, W=0, area, EOS, positivity, zero injection, and an
energy change of `4.26e-14` relative to its `.825` restart. It does not
replace the ongoing `.825→6.5` physical evolution.
