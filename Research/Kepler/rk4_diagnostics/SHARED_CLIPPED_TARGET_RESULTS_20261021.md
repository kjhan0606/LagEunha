# Same saved Kepler CONTROL: coverage and clipped integration connection

2026-10-21. Driver implementation/execution; one coherent Astra review in
`ASTRA_SHARED_CLIP_GEOMETRY_REVIEW_20261021.md`. This is progress toward the
original native CUDA/MPI/OpenMP Kepler objective, not its completion.

## Scientific status

The original saved epoch-1 finite kinematic CONTROL and its immutable heat
generator remain byte-identical, SHA256
`74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005`.
No new accepted RK4 evolution, restart/reset, pressure/heat alteration,
generator steering, viscosity, independent spin reservoir, commit or push
was performed. The production repository remains HEAD `8e25dfe0c2e41a7374c11e2e3102da2a8b028625`,
origin `git@github.com:kjhan0606/LagEunha.git`.

This diagnostic work is a CPU reference for geometry, physical moments and
error accounting. It does not certify production CUDA transfer or long
Kepler stability. The retained eight-ulp libm forward-error model is an
engineering assumption, not a formal CPU/CUDA library proof.

## 1. Boundary coverage completed, independently of integration fans

`shared_coverage_v2_mpi4` completed the winding/affine preflight:
384 interior queries, 1152 exact vertex queries, 384 half-open queries,
768 explicit endpoint events, 14880 changed-affine scalar comparisons,
zero failures. It retained 10 previous domains and completed another 24;
only the query families for sources 273/286 remained unproved.

The input coordinate intervals created root families wider than the fixed
root stopping threshold. Bisecting a family after its existence proof had
been discarded created phantom unresolved windows. The corrected canonical
root routine keeps a valid inclusion enclosure at the uncertainty floor.
The fixed256 work budget, original function, physical tolerances and
precision were unchanged.

`shared_boundary_witness_v5_mpi1` reuses these 34 domains/root packets and
completes ONLY the two outstanding source domains (14 membership queries).
Its merged `source_coverage.csv` is the authoritative ALL384 coverage file:

| Source category | Count |
| --- | ---: |
| Excluded | 344 |
| Whole source inside actual target68526 | 4 |
| Intersection candidate, certified empty | 10 |
| Actual nonempty clipped source | 26 |

Sources 118/246/254 all have complete boundary coverage. A complete boundary
is not a claim that an arbitrary physical star fan has a positive Jacobian.
Canonical roots and completed domains are reused by SHA. The original
817152-point cache/heat proof was not numerically replayed.

## 2. Four actual whole-source physical integrals completed

`shared_atlas_v4_mpi4` completes sources **240/244/248/252** (the actual CSV
IDs, correcting an earlier summary transcription). All 176 geometric
panels pass their same-function analytic certificate; all 45056 positive
common queries pass native-cache versus Q128 all12 checks.

The evaluator changes are exact algebraic reorganizations:

- Reference position is `A+s*((B-A)+t*(C-B))`.
- Stored P3 fields are affine part plus their exact nodal bends.
- Each source/edge residual uses a suitable exact Taylor anchor. Subtracting
  that cubic under exact P3 interpolation leaves the same `v0-I3(v0)`.
- The infinite analytic tail remains explicitly bounded at the original
  fixed256 cap. No polynomial refit or coefficient change is involved.
- The Gauss16 Cauchy bound removes a constant physical field times the
  exactly integrable linear reference Jacobian before bounding its variation.

These are FOUR-source integrals, not the complete target. They are retained
as one disjoint numerical group and counted once. Their potential error is
1.444114538894378e-21 from native/geometry/summation plus
4.319557566267092e-23 analytic remainder, against ORIGINAL whole-target tau
2.215823356869993e-21. The remaining sources must fit the same whole budget.
This partial bound does not establish either whole-target success or
impossibility.

## 3. All actual clipped-domain plans completed

`shared_clipped_target_v1_mpi4` first rejected 19/26 actual clipped sources:
one common original s/t axis was not certified transverse for all coarse
active faces. It did not attempt new analytic certificates or quadrature.

The following changes use the actual stored geometry before magnitudes:

- A signed restricted P3 Bernstein enclosure plus bounded analytic residual
  certifies inactive/excluding faces on each original chart.
- A gradient cone is centred on the actual stored local affine field,
  separating resolved rotation from its small nonlinear residual.
- The radial/angular axis is selected by its certified normalized
  transverse margin. This is a geometric choice, not an accuracy scan.
- Only the selected lower/upper boundary branches determine a slab metric;
  every original plane remains in the final physical ownership check.
- Real inverse Newton steps are damped when their full step would leave the
  source. The earlier unconditional half-slack restriction slowed already
  interior steps. The fixed16 iteration cap and final contraction enclosure
  remain unchanged.
- A noncontracting first interval-Newton graph inclusion is sharpened by
  certified monotone signs. A wide derivative cone is not immediately
  mistaken for an input-family uncertainty floor. The fixed256 cap remains.

One remaining chart, source298/chart5, is geometrically confined by original
face16 to original radial coordinate
`s <= 0.00146770090302326712042429588837485956`.
The certified derivative lower bound is
`3.44350531070710195769096843112610576e-05`, with exact apex face value
`-5.05403585409022956456202673234012469e-08`.
The radial tail is excluded by that original target inequality. A scaled
original-apex Duffy chart covers the entire remaining portion. Its basis,
stored coefficients and reference function are unchanged. No Cartesian
complex inverse is used at the apex.

All 26 source arrangements/plans complete in native v5 and later:
871 new panels, 222976 new positive queries planned, plus 45056 retained
whole-source queries = **268032**, within the original **817152** budget.
All 21528 independent exact-composed face value/AD2 comparisons pass.
The chart coordinate conversion handles actual scaled apex charts, with
positive orientation and corresponding determinant/Jacobian.

The progression v1/v2/v3/v4/v5 has 19/5/4/2/0 unproved source plans.
The partial 7297-panel v2 plan was too expensive and incomplete; it was
correctly rejected, not rescued by a larger cap. Branch-specific metrics
and margin-based axis selection reduce the complete plan to 871 panels.
`cost_budget_pass` is false until ALL plans exist, even if a partial count
would fit.

## 4. Physical arithmetic and publication contract

The original CPU reference explicitly used R80 fields and moments. The new
double graph is not called that same arithmetic contract. New clipped
queries are connected through `material_shared_reference80_ball.h` and
`material_shared_reference80_integral.h`: binary80 coordinates, field graph,
profile/all12 operations and per-source sums. Binary64 stored coefficients,
gradient/qnode/edge caches and their canonical uncertainty remain inputs.
No cache uncertainty is erased by conversion to R80.

`shared_reference80_v1_mpi4` checks 156 DIAGNOSTIC points on all 26 actual
clipped sources: 2652 field scalar and 1872 physical-moment comparisons pass
against the Q128 same-function evaluator. Binary record size is 1008 bytes,
including separate input/cache/field/moment/sum diagnostics. These fields
are cumulative diagnostics with different units, not an additive partition
of the total moment error. These are not target quadrature points.

The retained four-source group remains double. Thus a future assembled
packet is explicitly mixed arithmetic with Q128 geometry/global assembly,
not a homogeneous R80 replay and not a CUDA backend proof. Its entire
retained error is included once. There is no automatic precision scan or
reintegration of that group.

`material_shared_publication.h` accounts for final R80 value conversion and
outward error serialization BEFORE comparing all twelve ORIGINAL tau values.
The manifest is frozen from
`curved_whole_target_v3_mpi4/whole_target_ledger.csv`, not resampled scales.
Unit fixtures reject a conversion-only budget violation, reject incomplete
coverage and preserve the caller, and accept the deliberately loose fixture.
The thermal binding requires matching original coefficient/triangle/epoch
fingerprints; the generating recipe cannot reconstruct general evolved heat.

## 5. Remaining analytic gate

Native v5 completed all geometry/cost plans, but 22 clipped sources still
lacked complete analytic certificates. No new positive quadrature was run.

Native v6 introduces a uniform one-dimensional implicit correction about an
affine boundary predictor. Its residual is bounded using the directional
Hessian after first-order cancellation, and its correction radius follows
from the derivative gap. It keeps the existing rule, Cauchy radius,
function/tail/precision limits and physical tolerances.

In v6 **842/871** panels pass; 29 remain unproved:

- 27 binomial-domain enclosures in 18 sources exceed the unchanged 0.75
  Taylor-domain limit. For example source121/panel3 has argument centre
  0.01135350431319346 and radius 1.825187935147496. This is a failed complex
  enclosure, not a measured physical pole, negative J or solution failure.
- Two source301 panels have correction radius 0.00157600824172949 larger
  than their prior Hessian neighborhood. Their contraction 0.4428816764 and
  invariance 0.7214415800 are below one; the missing condition is neighborhood
  containment. It must not be described as a measured contraction failure.

Native v7 adds explicit linkage to the retained real root enclosure and
records the prior/real-root/correction radii. It retains **842/871** passing
panels: the same 27 binomial-domain and 2 prior-neighborhood rejections,
across 19 sources. Seven clipped sources have complete analytic bounds,
but global gating correctly prevents all new quadrature while any remain
unproved. Native runtime is 58.93 s (MPI4/OpenMP4). Its final actual result and
all integrity checks are recorded in `SHARED_CLIPPED_TARGET_LEDGER_20261021.json`.
Unobserved argument fields are tagged unobserved; zero placeholders are not
measurements. Per-panel pass/fail proof records are retained independently
of the per-source summary (which shows only the last failure).

Whole-target publication remains rejected until ALL actual clipped analytic
certificates, actual positive all12 integrals, whole twelve budgets and final
serialization pass. Caller epoch/target remains **9/1234** and both saved
state/heat checkpoints remain unchanged. Missing whole errors are `unproved`,
not zero. The complete target has not been numerically integrated yet.

## Next coherent implementation

Keep the completed geometry, roots, 176 whole panels/45056 actual integral
queries and 871 actual clipped plan. Inspect ONLY recorded failed analytic
neighborhoods first. Preserve affine dependence through the implicit graph;
the angular boundary can be represented by `q=s*t`, so its physical
Jacobian uses `q_upper-q_lower` and avoids artificial `1/s` derivative
cancellation. Exact apex factors require actual zero identities; a nonzero
apex offset and the genuine complex s=0 pole must remain. Select a finite
atlas from actual boundary/singularity/Hessian-neighborhood distances,
not progressively smaller guessed accuracy settings. Keep cap/order/precision,
tau, IC and guards fixed while correcting the analytic coordinates and bounds.
For the radial implicit branches, derive a valid Hessian neighborhood from
the actual transverse denominator and boundary singularity distances; a
correction cannot be accepted outside the neighborhood used to bound it.

After the complete analytic/cost gate, run the new clipped quadrature once,
combine the retained group once, and apply all twelve original budgets with
final serialization. A success then goes to full-target conservative
assembly/weight-law connection and accepted native RK4 evolution. It still
needs saved-state long-time t6/8/10 and the longest explicit GIZMO benchmark.

All builds/native computations in this bundle use lageunha, four MPI ranks
and four OpenMP threads per rank unless an explicitly single-rank witness.
No new Slurm simulation or monitoring loop was left running. No user input
or further per-test approval is pending. The long goal remains ACTIVE.
