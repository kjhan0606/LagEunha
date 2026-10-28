# Selective fixed-wall preclipping: native CUDA paired result

The accepted near32 fixed-polygon kernel adds 24 nearest-real-generator
halfplane clips before the coarse polygon wall for **every** cell. This
prevents transient corner overflow near the wall, but the great majority of
local seed squares cannot intersect the wall. The [candidate patch](KEPLER_FIXEDPOLY_WALL_SELECT_PRECLIP_2026-10-27.patch)
recomputes, on every local-square retry,

`hypot(x_i-3,y_i-3) + sqrt(2)*halfbox > polygon_apothem`.

When this conservative test is false, the square lies inside the polygon's
incircle, so the extra 24 preclips and both wall-cut loops are unnecessary.
The original eight directional preclips, full real-generator candidate
stream, complete 1,024-plane wall cut whenever the test is true, and the
corner-capacity abort remain. This changes the order and number of
intermediate clips, not the intended final cell geometry.

The comparison used two isolated syntax-local source trees. Both contained
identical fast CPU `Exam/exam.c` (SHA-256
`6b1d54c9de7ed7ce710293d92daac7e055434acf09a3674dd4026b9cde0c58f9`)
and were rebuilt with the same `make -j2 USE_CUDA=1 GFS_VARPRESS_W0=1`
flags and module stack. The control GPU source SHA-256 was
`2194c4aa81ea0eee030c5b07955794d0b7dc8d25dbc32f5049d14c71add67865`,
GPU object `552c4bdd57c360e2617fa6b439edf9b02cae583582cd1da7b41dc5414056c61f`,
binary `56d96440005d82b424f00041f0631540ab22628932b238df285a54216ba4bc7b`.
The candidate GPU source SHA-256 was
`9011dd6de77c2e938400e0fd7bce7dbbbedbd671eee028dabf066f643b498a6e`,
GPU object `1aaa579e232db023c00eb343bfdb0d83ee22f2060cbe8f51cd9ee22a2d4b34b2`,
binary `87d60c1cb059e1710dd2749c84af5c98ba17fe077e6fe080b9c8b42d5fd40968`.
These recompiled binaries are different from the live near32 binary; the
paired comparison holds their CPU source and build flags fixed.

The [stager](stage_kepler_bounded_wallselect_pair_20261027.py) and
[runner](kepler_q2_bounded_wallselect_pair_t085_a10_20261027.slurm)
restarted both from the **same accepted** `t=.825` native state (SHA-256
`134e4e76a475c6923f1106ffb944fd0d51a0df82e7d4c86ea2814057ce9a12f6`)
and ran to exact `t=.85` on one A10 allocation, first in control→candidate
order (job **437338**) and then candidate→control order (job **437340**).
Both jobs completed `0:0`. Each run passed its short native physics gate;
all ten ID-ordered parsed physical fields, including position, velocity,
mass, volume, internal energy, density, pressure and `W=0`, were **bitwise
equal** between control and candidate endpoints. The
[first pair receipt](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_pair_t085_a10_v1/pair_audit.json)
and [reverse pair receipt](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_pair_t085_a10_reverse_v1/pair_audit.json)
contain the complete field and timing audits.

| Run order | GPU tess median, control→candidate, ms | Whole RK4 median, control→candidate, s |
|---|---:|---:|
| Control then candidate | 177.4565 → 166.0840 | 1.188104 → 1.129439 |
| Candidate then control | 176.4865 → 165.4650 | 1.191964 → 1.138266 |

The 52 GPU tessellation calls in each run give candidate/control summed-time
ratios 0.9378 and 0.9382 in the two orders (first four calls excluded from
those sums). The 12 post-warmup RK4 steps give summed-time ratios 0.9487
and 0.9563. Thus the short paired result supports about **6% less GPU
tessellation time and 4–5% less RK4 time** on A10. It is not a sustained
`t=600` speed prediction.

The candidate emits about 30–50 fewer raw face rows per roughly 333,100
rows per tessellation. A separate one-step native capture **437342** now
resolves that difference at the accepted `t=.825` state. Both variants
passed their physical replay gates to `t=.826`, with all ten parsed physical
endpoint fields bitwise equal. The
[full-face equivalence audit](KEPLER_WALLSELECT_FACE_EQUIVALENCE_2026-10-27.json)
compares all `face_dump_r0_call1.tsv` rows by directed owner/neighbor/wall
support. The control has 333,100 rows and the candidate 333,056, but **both
have the same 333,052 support keys**. The control has 48 segments of length
at most `1e-14`; the candidate has four. Forty-eight support keys have a
different segment count. Once segments with the same support are combined,
the maximum face-length difference is `1.36e-16` and the maximum face area
vector difference is `1.36e-16`. Thus the 44-row difference in this first
call is extra negligible segments, not a lost finite face. The
[auditor](audit_kepler_wallselect_face_equivalence_20261027.py) and immutable
[native capture](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_face_call1_a10_v1/staging.json)
preserve the calculation. This is one accepted state; it does not prove
arbitrary late-time geometry.

The [separate corner-capacity stager](stage_kepler_bounded_wallselect_capacity_replay_20261027.py)
replayed the previously failed accepted `t=.045616533734626756` state to
exact `.05` as A10 **437339**, completed `0:0`. Its
[native receipt](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_capacity_replay005_a10_v1/audit_t005.json)
passes the geometry abort, bounded area, fixed-mass, positive-EOS, orbit,
zero-injection and energy gates; relative energy change from restart was
`3.14e-12`. This checks that the optimization did not reintroduce the
specific early 48-corner failure.

The live `t=.825→6.5` job **437317** and its dependent `t=10` job **437328**
retain their original accepted near32 executable. The candidate remains
isolated until a later accepted-state comparison and the main `t>6`
physics gate justify promotion.

## Evolved exterior checkpoint at t=5.775473

The original near32 job saved a complete
[17-check accepted state](KEPLER_Q2_BOUNDED_NATIVE_T577547_2026-10-27.json)
at `t=5.775473491592131` (snapshot SHA-256
`95ef285cc048bcf733876792d3d29d2917563d98b5915f9eb55cb44bf5044bdd`).
The same-compiler control and candidate binaries used above were restarted
from **that identical state** in a new A10 pair **437347** to exact `t=5.8`.
Both short native gates pass; all ten ID-ordered physical endpoint fields
are again bitwise equal. The [late pair receipt](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_pair_t58_a10_v1/pair_audit.json)
and [stager](stage_kepler_bounded_wallselect_late_pair_20261027.py) preserve
the inputs and comparison.

Across 48 post-warmup GPU tessellation calls, summed time falls from
`7.637828` to `7.182815` seconds (5.96%); over 12 post-warmup RK4
steps, summed time falls from `13.119127` to `12.446384` seconds (5.13%).
This extends short-window physical equality and the modest speed gain to
an evolved exterior state. It does not establish sustained long-time
equivalence; the live original `t=6.5` and dependent original `t=10`
jobs remain untouched.

## Guarded sustained promotion after the existing t=20 chain

The [t=20→30 runner](kepler_q2_bounded_wallselect_to30_a10_20261027.slurm),
[native gate/stager](prepare_kepler_bounded_wallselect_followon_20261027.py),
and [prestager](stage_kepler_bounded_wallselect_to30_20261027.py) are
queued as A10 Slurm **437348** with verified `afterok:437343`.
Job 437343 itself depends on the unchanged-binary `t=10` gate. At job
start, 437348 refuses to copy a restart unless the parent `t=20`
17-check native receipt passes, its endpoint SHA/time match, and its
executable is the expected fast-CPU/old-GPU binary SHA-256
`d92b6cd3b8491678cdfafa7610608378db89ca9bd16585a52a51aaf26f3b3367`.
It then stages the immutable fast-CPU/selective-GPU binary SHA-256
`87d60c1cb059e1710dd2749c84af5c98ba17fe077e6fe080b9c8b42d5fd40968`,
runs to exact `t=30` with 2.5-unit checkpoints, and applies the same
native 17-check gate. The physical settings remain those of 437343.
The parent **437343** subsequently passed all 17 checks at exact `t=20`
and **437348** began from that checkpoint. The accepted receipt and metrics
are archived in the [fast RK4 record](KEPLER_BOUNDED_FAST_RK4_2026-10-27/README.md).
The `t=30` endpoint remains pending its live run and audit.

## Compatibility with the earlier fast executable

An additional [read-only three-way comparison](KEPLER_WALLSELECT_LEGACY_FAST_COMPAT_T085_2026-10-27.json)
uses the existing `t=.825→.85` outputs of the earlier fast CPU/original GPU
executable (`d92b6cd3…`), the same-compiler rebuilt control (`56d96440…`),
and the selective-wall candidate (`87d60c1c…`). The executable and endpoint
hashes were checked, all three runs passed their 11 short native gates, and
their accepted restart SHA-256 is identical. All ten ID-ordered physical
fields (`x`, `y`, `vx`, `vy`, `mass`, `volume`, `ie`, `den`, `pressure`, `w2`)
are bitwise equal pairwise at `t=.85`. This removes the short-state binary
compatibility uncertainty for the queued `t=20→30` candidate continuation;
it does not establish arbitrary-state or sustained bitwise equivalence.
