# Ada result: retain the halo20 failure and validate an independent halo40 optimization

The dedicated Ada seven-case experiment completed all numerical evolutions
to their specified endpoints. Both executable qualifiers and all full-case
15-check native physics gates passed. The final strict seven-case comparison
**failed**, because two force-stage dumps had different vertex coordinate
bits when the halo was reduced from 40 to 20. The original failed receipt and
exit status remain unchanged. No halo20 performance follow-up or production
promotion was performed.

| Comparison | Full trajectory/state | Four force/RHS calls | Four geometry calls |
| --- | --- | --- | --- |
| control40 to vertex40 | Exact | Exact | Exact |
| control40 to vertex20 | Exact | Exact | Calls 2 and 4 differ |
| vertex20 to vertex20_reuse | Exact | Exact | Exact |
| vertex20_reuse to vertex20_reuse_local_sync | Exact | Exact | Exact |

All full trajectories use the same accepted t=2.6549183918680868 checkpoint,
59 steps to exact t=2.66, fixed mass and W=0. The endpoint x/y/vx/vy/mass/
volume/U/density/pressure/weight fields are bitwise equal in each comparison,
and every observed nonpositive/injection ledger field is zero. This does not
turn a failed geometric-reproducibility gate into a pass.

## Failure localization and independent reference

Only cell ID262 differs, among 147456 owned cells. The maximum coordinate
difference is 1.3877787807814457e-15 in call 2 and 1.1934897514720433e-15 in
call 4. These are real finite coordinate differences, not signed zero.
CSR offsets, support IDs, periodic image displacements and the complete
ax/ay/energy-rate/dissipation/dt arrays agree exactly. The other two force
calls' geometry agrees exactly.

A bounded read-only diagnostic constructs adjacent W0 bisector intersections
using exact rational arithmetic from the binary input generator coordinates
and recorded periodic images. Cell262 has six faces. The adjacent-plane
condition number in the infinity norm reaches about 164.89; no adjacent
pair is singular. Compared with these ideal-position intersections, the
largest errors in halo40/halo20 are respectively 2.13e-15/3.51e-15 for call2,
and 4.48e-15/3.29e-15 for call4. The minimum slack to the other recorded
support planes is about 3.03e-5. Direction subtraction is exact for these
recorded supports.

These references diagnose a well-defined ideal polygon and quantify both
outputs' errors. The dump does not contain the actual rounded CUDA threshold
coefficients, intermediate interpolation values or candidate clipping order.
The exact intersection reference therefore does not establish the error bound
of the native arithmetic path or prove that a changed candidate order caused
the difference. Unrecorded-generator constraints are not independently checked
by this small audit. No tolerance-based replacement acceptance gate was added.

Astra reviewed the completed failure and advised retaining its FAIL while
separately investigating supporting-plane conditioning and actual clipping
order if halo20 is pursued. A future geometry-equivalence criterion must have
independently derived arithmetic/error bounds and reject wrong supports/images
and unacceptable perturbations; a tolerance fitted to this observed delta is
insufficient.

## Independent halo40 fallback

Astra recommended testing halo40 with vertex certificate OFF, padding reuse
ON and device alias ON. This selects implementation optimizations separately
from the smaller-halo coordinate issue. One new native case was run using
the same a22 executable and checkpoint, while the already completed control40
was reused through a reference link. The baseline trajectory was not replayed.

The new combination passes exact comparison of all four complete native
geometry/operator dumps and all 59 accepted steps/endpoint fields. The native
physics/injection gates pass, and banners confirm both optimizations executed
without host alias fallback. Actual reference/current GPU UUID and hostname
are checked, with honestly distinct recorded execution sessions. This is a
separate passing halo40 receipt, not a revision of the failed halo20 test.

The [halo40 ABBA stager/auditor](audit_kepler_halo40_net_wall_20261029.py)
requires that independent exact halo40 receipt and a normal terminal status
before preparing any performance evolution. It then runs canonical671 halo40
and combined a22 halo40 in A–B–B–A order, independently restarting the same
checkpoint to exact2.66. Geometry dumps are zero; native endpoint output and
physics audits remain enabled. Four native endpoints add about 0.34 GB.

Bash time measures mpirun startup, evolution and endpoint I/O, excluding the
post-run Python auditors. All four runs must retain exact fields/trace/59-step
target, native physics gates and the actual requested execution branches. The
predeclared clarity rule is fastest canonical wall / slowest combined wall
greater than1.05; a favorable median alone is insufficient when repeats vary.
The finite Ada ABBA is executing at this record. Its final numerical and wall
results have not yet been accepted.

## Existing trajectory and promotion limits

Main438957, pending439115/439029 and the exact-t3 existing-allocation hook
remain unchanged. Astra advised preserving them: Ada does not establish A100
throughput or MPI scaling. Production promotion additionally requires target
device/rank validation, canonical-output wall evidence and a short bridge at
the latest accepted transition checkpoint, followed by a frozen continuation
and intact audit/restart chain. Do not interrupt an already started segment.

At scheduled t3, repeat the existing read-only locality diagnostic to assess
whether the exterior bottleneck persists. The earlier eight-cell exterior
result prioritizes boundary/IC consistency in later physics work. Individual
stepping, spin and weights do not themselves remove the counterstream cause.
The t>6–10, GIZMO t120/t600 and explicitly interpreted long-orbit requirements
remain unverified. This optimization work does not replace those endpoints.

Artifacts: [original failed seven-case receipt](KEPLER_VERTEX_HALO_ADA_RESULT_2026-10-29.json),
[rational-reference diagnostic](audit_kepler_vertex_halo_roundoff_20261029.py),
[diagnostic output](KEPLER_VERTEX_HALO_ADA_ROUNDOFF_2026-10-29.json),
[halo40 combination stager/auditor](audit_kepler_halo40_combo_20261029.py),
[separate halo40 pass](KEPLER_HALO40_COMBO_ADA_RESULT_2026-10-29.json),
[halo40 case stage](KEPLER_HALO40_COMBO_ADA_STAGE_2026-10-29.json),
[ABBA stage](KEPLER_HALO40_NET_WALL_ADA_STAGE_2026-10-29.json).
