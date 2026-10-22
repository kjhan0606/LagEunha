# Explicit curved topology events and positive source-complement assembly

2026-10-21. Native C++/OpenMP/MPI, authorized manual CPU execution on lageunha.
Same preserved finite h=2e-4 geometry CONTROL. No time evolution, restart,
new IC, artificial viscosity, thermal correction, weight solve or compression.

## Event isolation implementation

`material_curved_events.h` evaluates current compiled position, gradient and
Hessian with long-double intervals. It enumerates source vertices/reference
radial extrema, source-edge/power-plane crossings, active power vertices,
radial tangencies and angular-chart crossings. Two-dimensional preconditioned
interval Newton and one-dimensional interval Newton retain all surviving
boxes. Candidate/active pruning uses whole-box inequality bounds, not samples.
Fixed caps:4096 box evaluations per2D system,1024 per edge, depth48. These
are new event-isolation budgets, not an increase of the failed radial guard's
256-panel cap. Budget failures retain witnesses and prevent acceptance.

Native position/gradient envelopes are retained, including their explicit
engineering assumptions. The compiled polynomial is an exact algebraic
reexpression of stored coefficients; no future orbit or field reset is used.

Every possible branch-order event is placed in an explicit radial band.
Only in those bands the old interval branch guard encloses unresolved angular
windows. A normal simple branch's swept range is a resolved region, not
unknown topology. Radial measures use (high-low)*(high+low) to avoid subtracting
nearby squares. Nonzero uncertainty is retained throughout.

The reference domain is the triangle of stored vertices. New source predicates
use an interval inverse of those edges. Stored double geometry gradients
remain unaltered basis coefficients; they are not used to redefine the exact
reference-domain edges. The original data/source snapshots remain unchanged.

## Completed evidence

Same fixture throughout: source triangle(0,0),(1,0),(0,1), actual map
chi=(X+.05Y^2,Y), J=1, constant v=(.1,.2), original pressure/density, a closed
equilateral power cell of area3*sqrt(3)*1e-8 wholly inside the curved source.
Its original area tolerance is1e-11; event mass/K budgets remain unchanged.

- `curved_lens_v3`: old radial guard exhausted256 evaluations; area
  error1.5154e-11 and event mass bound1.2431e-6. Failed, preserved.
- `curved_events_v1`: explicit isolation finds15 event boxes in71 evaluations,
  area error1.9788e-21. Charging the entire event annuli overestimates the
  unknown region: mass bound3.7868e-11 exceeds5e-13. Failed, preserved.
- `curved_events_v2`: enclose only unresolved angular windows inside the
  isolated event bands. Same15 boxes/71 evaluations; area error1.9788e-21,
  event mass bound7.3451e-17 vs5e-13 allowance, K bound1.8363e-18
  vs6.9433e-13 allowance. Terminal exit0/AUDIT_COMPLETE,1.618s.
- `curved_events_v4`: exact-vertex source predicates; same fixture and gates
  pass, area error1.9729e-21. Terminal exit0/AUDIT_COMPLETE.
- `curved_events_v5`: certified pruning of inactive bisectors and source
  bands wholly outside the source; same fixture passes, area error1.9738e-21,
 1.781s. This is the fixture gate for subsequent actual controls.
- `curved_events_v3`: first actual W0 packet audit, all384 triangles pass;
 8578 evaluations,2930 event boxes, no budget failures. Total event mass
 bound5.079e-23, K bound4.711e-23; maximum allowance ratios7.04e-4
 and3.08e-4. MPI4/OpenMP2, .742s, terminal exit0/AUDIT_COMPLETE.
 It preceded the exact-vertex source-predicate correction; the actual audit
 must also be repeated with the corrected implementation before transfer.
- `curved_events_v6`: corrected-domain W0 audit passes all384 triangles,
 8578 evaluations/2930 boxes, no budget failures, .751s.
- `curved_events_v7`: existing nonzero weight pattern also passes all384;
 5332 evaluations/1684 boxes, no budget failures, .577s. It is a prescribed
 pattern, not a converged dynamical weight solve.

The first full outgoing/complement execution `curved_complement_v1` failed
on the positive-weight assertion. It is terminal (FAILED_EXIT_255); Intel
MPI terminated remaining ranks after one worker failed. Partial positive
packets, binary and logs remain preserved. The saved origin event boxes show
an exactly zero input vertex spuriously expanded to radius1.0457e-2475 by
unconditional nextafter widening of zero products/squares. Its quadrature
weights underflow. This is an interval arithmetic artifact, not a physical
overlap that may be deleted. Exact zero addition/multiplication/square/root
now remain exactly zero, and opposite exact point inputs cancel exactly.
No threshold clips a nonzero interval. The same fixture/actual event audit
and a single origin-triangle regression precede the next full execution.

These are EVENT bounds, not full cold-remap acceptance. For example v5 still
reports the separate angular coefficient/native uncertainty at radial nodes.
Neither its small area error nor global conservation certifies full per-target
cold precision, shared C0 gluing or positive complement mass matrices.

## New complement implementation, acceptance pending

`test_curved_complement_v1.cpp` independently integrates the actual whole
source at24x16 and computes only outgoing positive overlaps. Each shared
outgoing record is used once for subtraction from its source and addition
to its target. Original source plus target half-plane predicates define the
positive own-cell complement; no negative physical quadrature is stored.
The baseline is evaluated from the current packet, not read from old target
results. It does not restore the original IC or force individual target masses.

Positive8x12/16x16 rules estimate resolved-region quadrature error. Geometric
uncertainty contributions are kept separate. The fixed per-triangle budget
is initial event panels plus256 evaluations, depth24. Target m,p,L,K,U,
position/polar moments, potential and thermal-rate ledgers are saved. Resolved
I, intrinsic spin, Kint,Krot,Krest are diagnostics of the same retained field;
there is no independent rotational energy reservoir. Accepted estimate output
is explicitly labelled NOT_COLD_ACCEPTANCE. Mass-matrix positivity, shared
gluing and full cold error certification remain required before weight solve
or compression, followed by native CUDA/MPI t6/8/10 and longest-GIZMO evolution.

## Completed complement evidence and next restricted bundle

Exact-zero arithmetic passes the same tiny-cell fixture (`curved_events_v8`),
all384 W0 event audits (`v9`) and all384 existing weighted-pattern audits
(`v10`). `curved_bounds_v5` checks9600 native queries, excess0, unchanged
positive J/boundary separation results. Single origin triangle133 passes
the positive-weight regression in4.547s; there is no small-overlap cutoff.

The full complement's first completed data were recovered after an operational
runner-editing error caused an accidental duplicate and truncated native
stdout. See `curved_complement_v2_first_completed_recovery/RECOVERY.md`.
The first exit code is unknown; no AUDIT_COMPLETE is manufactured. All384
triangle rows and24 target/576 source-outgoing records are complete. All4
positive packets are independently intact,4,575,546 points, saved outgoing
mass roundtrip relative9.322e-20. Source vs independent finite control:
mass6.716e-13, K/coldU2.003e-11, U6.432e-13. Target m,U,I are positive;
minimum resolved Krest2.6732e-12. These do not prove complement-matrix positivity.

Numerical failure remains:344/384 triangles meet the original quadrature
estimate;40 exhaust the fixed budget. Worst ratio42.0018815736 is triangle225.
Its independent immutable-runner witness reproduces the same268 evaluations,
70,024 points and explicit exit2 in7.429s. Sum of per-panel normalized
estimates:potential42.00188, K4.17185, mass.32713, angular.17731, U.00327.
Its angular K upper at radial nodes is7.603e-9 of that source triangle's U.
This is not a full continuous radial/target-specific error certificate.
The worst potential-panel known-rule differences~1.2e-22 are covered by
the saved boundary-band bounds~1e-22--1e-20, making window-induced estimator
noise plausible but not established as the sole cause.

`curved_cache_equivalence_v1` verifies that caching and whole-source certified
target pruning preserve triangle133 moment CSV and point bytes exactly;
old4.5552s/new1.5227s, an observed2.992x gain for this control only. No full
solver/GPU or full-batch scaling claim is made.

The new bundled Astra decision is `ASTRA_CURVED_BRANCH_PLAN.md`. Implement
one common smooth representative for both rules, retain separate uncertainty,
add the missing THERMAL_RATE band bound, and attribute all error to targets.
Only225/137 plus the original tiny-cell regression next; no full rerun, weight
Newton, compression or evolution before that restricted bundle passes.
