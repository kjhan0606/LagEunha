# Common-branch integration and the radial parameterization witness

2026-10-21. Driver implementation/execution on authorized lageunha, native
C++/MPI, one rank and one OpenMP thread per independent control. No physical
state, PDE, IC, force, heat, spin or weight was changed. The h=2e-4 saved map
is a finite geometry control, not a new Kepler evolution or a restart.

## Completed bounded bundle: failed nominal/cold acceptance

The same scalar root representative is used by both positive rules; interval
root windows certify uncertainty rather than defining the integration cut.
Sources, current map coefficients and the native physical integrand are
unchanged. Actual branch IDs, root bands and per-panel moments are retained.

| Control | Panels | Old known ratio | Common representative ratio | Trace failures / fallbacks | Native exit |
| --- | ---: | ---: | ---: | ---: | ---: |
| triangle225 | 140 | 42.00188157 | 40.45587790 | 52 / 4 | 2 |
| triangle137 | 138 | 17.44412255 | 16.84904307 | 52 / 4 | 2 |
| same closed tiny cell | 44 | 13981709.21 | 13981708.86 | 77 / 10 | 2 |

Folders: `curved_branch_v1_225`, `curved_branch_v1_137`,
`curved_branch_v1_lens`. Every folder has immutable runner/environment,
source/binary/input SHA evidence, explicit native exit and failure status.
All regular branch-order comparisons pass; native velocity-reference-gradient
queries have zero enclosure excess. These facts do not imply a cold remap
pass. Whole-panel root tubes remain unresolved on120/134 branches.

The tiny cell is the original delta=1e-4 fixture with exact current area
3*sqrt(3)*1e-8 = 5.196152422706632e-8. Its new fine area error is2.641e-23,
so its small closed region was retained accurately. Its complete outgoing
moment estimator and trace gates failed; area alone does not certify them.
The 1e-11 figure in older narrative summaries is the area ERROR TOLERANCE,
not this fixture's exact area. No smaller/easier fixture was substituted.

The common representative did not explain away the previous failure.
Potential-energy integration dominates, and regular panels contribute
40.45582643/16.84904307, versus event-panel contributions5.15e-5/3.30e-9.

## One-bisection and arithmetic diagnosis: new actionable evidence

Only the two captured worst regular panels were investigated. No subdivision
cap, tolerance or rule order was raised. Original8x12/16x16 rules and the
native integrand were retained. Binary128 is a reference calculation only.
All variants integrate the same physical region/fields; the linear-radius
variant changes the quadrature parameterization with its positive Jacobian.

Triangle225 panel19: r=[.0234618747872393723954,
.0234653568997021398807]. Triangle137 panel49:
r=[.0117315880378319244546,.0117337277108039118632].

| Diagnostic | triangle225 panel19 | triangle137 panel49 |
| --- | ---: | ---: |
| Original sin-squared radius, coarse/fine ratio | 1.623706315 | .4104097851 |
| Radial-order change only,8x16 vs16x16 | 1.623704394 | .4104097091 |
| Angular-order change only,16x12 vs16x16 | 1.922e-6 | 3.707e-8 |
| Binary128 ray/root, original radius map | 1.623711570 | .4104064837 |
| Stable binary128 density too, original radius map | 1.623579495 | .4103821805 |
| Linear radius, ordinary long double | 7.227e-5 | 1.349e-5 |
| Linear radius, binary128 ray/root | 7.184e-5 | 1.079e-5 |
| Linear-radius fine vs one-bisection fine | 1.223e-5 | 1.906e-5 |

Ratios retain the original source tolerances, including1e-10*U_source for
K/potential/U and the original factor8 estimator. They are numerical
estimates, not certified true error or target-U-normalized full-transfer
uncertainty. Binary128 root shifts are only~2–4e-18 in chart coordinate.

This identifies the radius-rule error as dominant in these witnesses.
The universal sin-squared endpoint map can itself make an otherwise smooth
regular panel expensive for the low rule. Ordinary linear-radius quadrature
reduces the captured discrepancies by~22467/~30418 without more points.
This is not a proof that a linear map should replace endpoint regularization
at actual tangencies. Event/tangent intervals still need separate treatment.

Complete raw evidence: `curved_precision_v3_225` and
`curved_precision_v3_137`, each native exit0 and explicit
DIAGNOSIS_COMPLETE_NOT_PHYSICS_ACCEPTANCE. Timings .984s/2.047s exclude
compilation. `curved_precision_v2_225` retains the earlier four-mode witness.
`curved_precision_v1_225` is a preserved BUILD FAILURE: nested-main/ADL
errors, native program never started; no fabricated native exit/status.
The library guard, qualified call and explicit build-exit logging are fixed.

## Missing uncertainty terms remain explicit failures

The unknown-band THERMAL_RATE bound is now populated using the current
velocity-reference-gradient enclosure and div(v)=tr(V F^-1). No zero thermal
rate bound is inferred from a missing calculation.

The v1 whole-panel `analytic_delta_p` is a PROTOTYPE, NOT a complete certified
representation bound: it omits residual-B coefficient uncertainty and
source/power normal/constant conversion terms. The target contribution CSV
also lacks a complete field-rounding/gluing ledger and certified target-U
normalizers. Its numbers must not be used as complete error upper bounds.
The JSON explicitly has target_gluing_and_field_rounding_complete=false and
full_cold_transfer_proved=false. Positive m,U,I is not matrix positivity.

This bundle does not authorize full384 transfer, weight Newton, compression
or Kepler evolution on the strength of the two good panel estimates.
The full goal, including native CUDA/MPI t6/8/10 and the longest GIZMO
duration, remains ACTIVE and unproved. The next bundled review should select
a parameterization-aware bounded integration plus complete target error
ledger, rather than a precision/order/cap scan or a repeated IC refinement.

## Geometry-selected implementation v2: large improvement, still failed

The subsequent ONE bundled Astra review selected
`ASTRA_PARAMETERIZED_TRANSFER_PLAN.md`. That plan remains incomplete; the
geometry stage has explicit failures and the target field/gluing ledger
must still be finished before acceptance.

`material_curved_branch_bounds.h` now certifies a whole-radius root tube by
opposite endpoint signs and a derivative separated from0. The derivative
intersects two valid enclosures (direct Horner and Bernstein). It encloses
both original and fixed-representative coefficient families, with native
boundary uncertainty retained. At most12 bracket expansions replace the old
2048-node full-chart family isolation; no subdivision cap was increased.
Direct endpoints avoid interval dependency accumulation that had merged two
distinct simple branches into a single unresolved swept band.

Conversion bounds now include cubic AND residual B, exact-vertex source
inverse/constant, and power normal/constant/weight scalar arithmetic.
The n_exact*delta_chi term includes the dn*delta_chi cross term. Numerical
ray construction uses a separate pre-cancellation coefficient-norm envelope
(<16384 long-double operations, round-to-nearest and <=8ulp powl assumption),
plus root stopping uncertainty. These are engineering models, not a formal
libm proof. Angular strips use the conservative factor4 from a two-sided
chart-coordinate band and dtheta/du<=2.

Only event-free panels with these tube checks receive linear-radius
quadrature. Unresolved/event panels retain the old diagnostic rule and
remain explicit failures; no unconfirmed fold map was adopted.

| Entire captured triangle | Linear panels | Common estimator v1 -> v2 | Unresolved regular branches v1 -> v2 | Exit |
| --- | ---: | ---: | ---: | ---: |
| 225 | 129/140 | 40.45588 -> 1.0699255 | 120 -> 6 | 2 |
| 137 | 128/138 | 16.84904 -> .2098453 | 134 -> 4 | 2 |

Regular branch order failures remain0; native velocity-gradient enclosure
excess remains0. The original immutable partitions/rules/fields are used.
Nominal estimator<1 for137 does not certify geometry or the target thermal
ledger. Both complete controls are FAILED. No full384 remap was repeated.
`CURVED_PARAMETERIZATION_WITNESS_20261021.json` retains every remaining
regular root failure and the three dominant panels per triangle.

The same tiny-cell v2 regression retains its exact area to1.678e-23, with
24/44 linear panels, but its complete moment/trace checks still fail (exit2).
It is not relabelled a whole cold pass because its area is accurate.

## Certified folds: actual next geometric distinction

`test_curved_fold_v1.cpp` uses the SAME fixed current-map coefficients.
Point Newton supplies a candidate only. A single fixed1e-12 validation box
and strict Krawczyk inclusion prove a unique representative fold under the
recorded interval/libm model; f_r and f_theta_theta intervals exclude0.
There is no box-width scan, physical state change or new time evolution.

| Witness | Unique simple fold radius | Krawczyk contraction | Original source domain |
| --- | ---: | ---: | --- |
| 225, plane14/15 | .0234653568997044901 | 1.927e-12 | strictly INSIDE |
| 137, plane8/7 | .0117095453912223958 | 1.196e-12 | strictly OUTSIDE |

Both native fold diagnostics exit0, explicitly
FOLD_WITNESS_COMPLETE_NOT_COLD_ACCEPTANCE. `fold.csv` preserves all root,
derivative and inclusion intervals. The225 fold agrees with the preserved
event box. The137 fold lies outside the fluid triangle and was correctly
excluded by the existing source-filtered event guard; its absence from that
guard is not evidence of a missed physical event.

An independent exact-rational half-plane check uses the binary source
vertices and expands the printed fold boxes by1e-20 in each coordinate. It
proves the inside/outside statements above. Frozen source/binary SHA checks
for all11 completed branch/precision/fold controls match.
`CURVED_BRANCH_INTEGRITY_DOMAIN_20261021.json` retains those checks.

Next geometric actions follow the selected plan: a confirmed-fold
parameterization for225 with the critical-radius uncertainty retained;
source-clipped branch inclusion for137 instead of treating an exterior
fold as an interior event; exact handling of source-vertex radial extrema.
Then complete target-U lower bounds and field/gluing/summation budgets.
Do not extend this failed geometry to Newton/compression/runtime, raise the
cap, repeat the arithmetic sweep, or call a short control long-time success.
