# Original-function / operation-graph attribution: completed fixed bundle

2026-10-21. Driver implementation and execution on lageunha. This is the
single saved225/137/lens attribution bundle selected in
ASTRA_ERROR_ORACLE_PLAN.md and ../ASTRA_NATIVE_TRANSFER_PLAN_20261021.md.
No time evolution, new IC, new fixture, quadrature order, subdivision budget,
precision scan, coefficient refit, trace averaging or thermal repair.

## Decision

**Do not promote to full24-cell transfer, weight Newton, compression or
physical time stepping.** The independent point models pass, but the known
boundary error model exceeds the tiny target's entire physical tolerance.
Source-reference integration and target-local gluing certificates remain
absent. Missing terms are null/unproved, never zero. The original native
CUDA/MPI long Kepler goal remains unachieved.

## What was implemented

- `material_forward_ball.h`: operation-level forward errors with outward
  radius arithmetic; fixed binary128 support, explicit IEEE/libm assumptions.
- `material_original_oracle.h`: independent ORIGINAL stored-coefficient P3
  and residual basis oracle; stable local-cubic + Bv0 FMA Horner and an
  independent midpoint-representative check. No rounded compiled polynomial
  is substituted for the original physical reference.
- `material_native_forward.h`: actual original CPU basis/map operation graph,
  including the radial series stop/tail, stored-gradient compatibility and
  cast-coordinate radius, J and inverse, velocity gradient.
- `material_ray_forward.h`: actual degree-six scalar construction and Horner
  range forward model on each SAME accepted radial panel and root tube.
  No16384-operation blanket multiplier. Rounded plane scalars are exact
  inputs here; their construction error stays in the separate conversion term.
- `test_curved_error_ledger_v1.cpp`: geometry-only possible incidence fixed
  before reading errors, same saved-point native sum, outward serialization
  bridge, actual binary80 summation depth and nonzero whole-edge gaps and
  tangential derivatives. All possible off-diagonal contributors are kept,
  including saved zero flows. No observed-zero exclusion.

Normal IEEE nearest arithmetic and elementary functions within8 ulps are
explicit engineering assumptions. This is not a formal CPU/GPU libm proof.
The stable FMA function is an audited expression, not a change to the
executing native CUDA kernel or to any saved nominal moment/packet.

## Execution and invariants

| Case | Same saved positive points | Point graph / oracle excess | Degree-six graph checks | Degree-six mismatches |
| --- | ---: | ---: | ---: | ---: |
| 225 | 37567 | 0 / 0 | 1960 | 0 |
| 137 | 34498 | 0 / 0 | 2080 | 0 |
| lens | 32289 | 0 / 0 | 284 | 0 |

All104354 points checked, including all13 map/gradient fields and all12
moment expressions. The original13.72/12.63/10.54-second point audits and
the final1.42/1.80/.51-second attribution runs exited0. The final original
packet, saved positive quadrature and moment/geometry CSV SHA checks pass.
Every run used a frozen runner/environment and source/binary SHA ledger.
All processes are terminal; no monitor or background job remains.

Retained failures: original-oracle v1_225 and error-ledger v1_225 are BUILD
failures only. Error-ledger v2_225 is a retained preliminary attribution;
v3 separates the already charged plane-conversion path and preserves
possible zero-flow rows. These are code/model corrections, not resolution,
quadrature or precision parameter searches. The first independent point
audit is not repeated merely to rewrite output metadata.

The point audit's `field_forward.csv` has two groups named `stable`: the
first13 stable rows are representative-evaluation error, the next13 include
original-to-compiled conversion. The immutable file is preserved.

## Actual dependency paths

Regular labels/cuts use the compiled representative. Original conversion,
the actual ray arithmetic and the retained root guard move that boundary.
Charging the native integrand position envelope AGAIN as a regular-boundary
movement was a duplicate; it is removed here. Native position still enters
potential/angular/position moments. The original event bounds, including
unknown native owners, remain unchanged. Axis and source-fold certificates
retain their original separate arithmetic/coarea derivation.

Pre-save binary80 rho,w are enclosed within half a binary64 ulp of the
stored values; all moment operations are then outward enclosed. The old
sum error uses gamma_(target points + accepted panels +4) times the absolute
term sum, not a canceled global sum. Q128 canonical sums recover the
original-basis reference using the already recorded signed Q128 difference;
old/new Q128 accumulation, the binary80 difference cast, oracle radius and
subtractions are included. No new original-basis point solve was needed.

## Fixed possible-incidence budgets

The physical384-triangle state gives16-72 possible elementary contributions
per physical target. Count source baseline once, every possible own outgoing
subtraction once, and incoming contributions once. The source's represented
partial allocation counts every possible off-diagonal row, not only points
with observed flow. The lens has1 incoming contribution to101/103/104 and4
baseline/outgoing contributions to102. Gluing/source certificates remain
additional missing terms; this allocation is not a full acceptance proof.

Use the maximum ALREADY proved U lower bounds (never add overlapping proofs).
Target energy tolerance remains1e-10 U_lower and Udot tolerance1e-9 U_lower.

| Case | Worst incoming row | Known bound / partial allocation | Known bound / entire target tolerance | Largest measured same-packet difference / partial allocation |
| --- | --- | ---: | ---: | ---: |
| 225 | 68526 potential | 4.32180548 | .06752821 | .00044832 |
| 137 | 67942 potential | 2.81348978 | .06498082 | .00093096 |
| lens | 101 potential | 211.59890098 | 211.59890098 | .00105587 |

These upper bounds are NOT observed physical errors. Measured differences
compare the same discrete quadrature packet's native moments with the
original-basis reference, not the exact continuous transfer. Partial
allocation failure for225/137 alone does not prove their full-target
physical tolerance is exceeded. The lens fails even the full tolerance
under this error model, before adding the missing source/gluing terms.

Lens101: U_lower=6.9848668650464756e-13, tolerance=6.9848668650464756e-23.
Known regular-boundary bound=1.4718257046682706e-20; field=5.7698e-24,
serialization=4.8251e-24, sum=1.5213e-24, retained event=4.9528e-23.
Observed same-packet potential difference=7.37513109038022e-26.
Its quadrature sensitivity ESTIMATE=1.0442508999013603e-22 also exceeds the
physical tolerance; it is not a rigorous integration remainder proof.

## Dominant witness and required expression accuracy

`CURVED_BOUNDARY_DIAGNOSTIC_20261021.json` decomposes saved CSV bounds by
recovering a per-panel field envelope from the old bound / summed old
measure. This quotient is diagnostic, not a new outward certificate.

For lens101 potential the regular-boundary contribution is approximately:

- actual expanded-ray arithmetic:4.81533e-22;
- original->compiled / exact->rounded plane conversion:3.59328e-21;
- retained128-epsilon root guard:1.06434e-20.

The largest witness is accepted panel40, method2. Its numerator bound is
1.7978194765661250e-18. If EVERY boundary contribution scaled uniformly,
the boundary factor would have to be <=.000557416155, corresponding to
this numerator <=1.0021336201939550e-21. This optimistic design target
ignores the still missing source/gluing and the quadrature estimator; it is
not an assertion that replacing a guard by that number is valid. An
observed maximum cannot be substituted for a forward error model.

For225/137 the analogous optimistic uniform boundary factors are.08619855
and.34509716. Simply tightening native-field arithmetic cannot by itself
close these current root/conversion/axis/coarea bounds.

Touched whole-edge position/velocity gaps are nonzero:225 up to
2.89659e-18/2.07425e-15;137 up to2.57031e-18/9.11110e-16. Tangential F/V
gaps are recorded, as are deliberately broad whole-source reference V
jumps. Global exact C0 and target-area ambiguity are unproved. The
standalone lens has no paired material edge and cannot certify the
physical state's gluing by itself.

## Artifacts / next decision boundary

- `curved_original_oracle_v2_{225,137,lens}`: original point audits.
- `curved_error_ledger_v3_{225,137,lens}`: boundary, incidence, bridge, sum,
  touched-edge ledgers, contracts, frozen input/source/binary hashes.
- `CURVED_ORIGINAL_ERROR_LEDGER_20261021.json`: fixed target allocations,
  all known components, measured differences and explicit missing terms.
- `CURVED_BOUNDARY_DIAGNOSTIC_20261021.json`: dominant witness decomposition.

The selected bundle has reached its specified failure-report boundary.
No full24 transfer or physical timestep was released. One existing-Astra
evaluation may now select the NEXT coherent implementation bundle; do not
repeat precision/order/cap scans or replace the native CUDA/MPI long goal
with an endless sequence of isolated formal fixtures.
