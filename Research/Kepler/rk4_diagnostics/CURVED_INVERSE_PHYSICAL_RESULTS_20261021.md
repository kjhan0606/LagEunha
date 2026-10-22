# Physical225/137 inverse-target bundle: local terms closed, full target next decision

2026-10-21. Driver follows the SAME bounded
`ASTRA_TARGET_INVERSE_PLAN_20261021.md` after the tiny-target gate passed.
No additional Astra microreview during implementation. Read the JSON
ledger for exact values and preserved-run contracts.

## Decision and practical limits

The original nonzero-residual225/137 source functions now have native
physical-target quadrature, analytic integral remainder, source baseline,
and explicit shared-edge domain/field ambiguity envelopes. All selected
local terms are computed; source partition conservation is checked for
all12 moments. No original coefficient, state, weight, timestep or old
positive packet was changed. No ray/polar root, refit, heat correction,
steering or viscosity was introduced.

**Full-target acceptance is still unproved.** Worst local error / fixed
possible-incidence WORK allocation is2.38574 and2.17623, both dominated by
the shared material-cell border envelope. The same errors consume only
3.72772% and3.75212% of the unchanged WHOLE target tolerance. These are
upper-bound/work-allocation failures, not measured physical errors or
proof of a whole-target failure. Do not silently change the allocation or
normalizer to declare these partial rows passed.

One end-of-bundle Astra evaluation should choose between ONE canonical
whole-target assembly (every possible contributor, shared-edge charge once)
and terminating promotion of this broken representation. Do not start
another sequence of local guard/order/precision scans. No weight solve,
compression, transactional commit or physical timestep was released.

This is the same24-cell finite control, with384 material triangles, not
a long128^2 Kepler simulation. Production CUDA/MPI t6/8/10 and the GIZMO
longest duration remain unachieved.

## Actual native implementation

- `material_target_general.h`: algebraic expansion of the ORIGINAL P3+
  initial-v0 residual map into local cubic+Bv0, with Q128 coefficient
  intervals. Actual native binary80 FMA evaluation, F,V,J and velocity
  gradient; original-basis residual/compiler/native comparisons; residual/
  inverse-conditioning bounds and exact geometric source membership.
- `material_target_domain.h`: all source-edge/target-plane intersections,
  exact fixed-expansion determinant fallback for parallel power planes,
  feasible target vertices, oriented curved/straight boundary loops,
  strict interior anchor and whole-boundary positive star-fan Jacobian.
  Empty receivers require complete boundary exclusion and a target interior
  witness outside the source image. Saved zero flows are not exclusions.
- `material_target_analytic.h`: complex original-function map/gradient/
  Hessian, positive-J and radial/profile/potential branch majorants;
  anisotropic inverse contraction and a priori Gauss16 remainder bounds.
- `material_target_gluing.h`: original neighboring source traces restricted
  algebraically to the SAME reference edge before bounding their difference;
  unresolved ruled-strip area, endpoint caps, target exclusions and both
  source field envelopes. No arbitrary owner or manufactured glued state.

All kernels are native C++; Python only reads ledgers and uses exact decoded
binary arithmetic for bookkeeping. Fixed binary128 geometry/inverse fallback
and binary80 field/sum execution are CPU contracts on lageunha. No new CUDA
kernel performance or physical stability claim follows from these runs.

## Boundary, field and inverse checks

`curved_target_general_v1_{225,137}` both exit0:

| Source triangle | Boundary/plane intersections | Original/compiler/native checks | Model excess | Whole-source chi-I Lipschitz bound |
| --- | ---: | ---: | ---: | ---: |
| 225 | 6 | 468 | 0 | .01841450075 |
| 137 | 16 | 1066 | 0 | .02024268243 |

Whole-source-box J lower bounds are.99833085695 and.99891256290.
Maximum certified inverse radii at queried boundaries are7.68e-32 and
2.86e-31. The compiler retains original stored gradient compatibility and
the true mathematical partition-of-unity basis; it is not a fitted state.

## Positive physical integration and coverage

`curved_target_integral_v2_{225,137}` both exit0,8.43/8.27s:

| Source | Nonempty receiving target | Proved empty receivers | Closed loops | Boundary pieces | Positive points |
| --- | --- | ---: | ---: | ---: | ---: |
| 225 /68525 | 68526 | 1 | 1 | 3 | 3072 |
| 137 /67943 | 67942 | 2 | 1 | 3 | 3072 |

Each loop has two actual source curves and one straight target face. The
whole star Jacobian lower bounds are4.81769e-8 and1.06863e-8. Every physical
point is certified inside the original source. Four fixed sectors on each
boundary, each with a positive16x16 fan rule, give3072 points. These are
dimensionless local fan coordinates, not source polar charts. All12 moment
expressions, including thermal rate, use the SAME physical positive points.
Independent original nodal-basis/profile oracle excesses:0.

The original smooth reference density is evaluated with an equivalent
stable logistic expression in the native kernel. Its original tanh form
remains the independent oracle; no density profile or initial state changes.

## Actual analytic integration remainders

`curved_target_analytic_v1_{225,137}` both exit0. All12 sectors per source
pass the complex branch/J/inverse-box proofs. Fixed Cauchy radius3/4,
Gauss16, no precision/order/cap/tolerance scan. Separate complex-coordinate
discs are enclosed by an anisotropic inverse box; the domain also excludes
the density radial square-root/logistic poles and the physical potential
square-root branch. The curved-fan determinant is bounded using
J'=cross(gamma-O,gamma''), preserving the exact cancellation of the tangent
cross itself. No observed maximum replaces a bound.

As in the tiny-target report, the Gauss/Legendre remainder and Cauchy
derivative estimate give2*C16*R^-32*M, C16=(16!)^4/[33*(32!)^2]. See
[NIST DLMF3.5.19/21](https://dlmf.nist.gov/3.5#E19). No coarse/fine rule
difference is used as a rigorous remainder. Positive points are NOT
recomputed merely to add this analytic ledger.

| Target contribution | Native potential | Native field/geometry/node/sum error | Analytic potential remainder |
| --- | ---: | ---: | ---: |
| 225 ->68526 | -1.5402041928956028e-9 | 1.31464e-25 | 9.18168e-24 |
| 137 ->67942 | -4.7761432393706263e-10 | 4.11465e-26 | 2.92150e-24 |

These are enclosures under normal IEEE nearest arithmetic and the explicit
<=8ulp elementary-library engineering model, not formal CPU/GPU libm proofs.

## Source baseline and original shared traces

Final`curved_target_source_glue_v2_{225,137}` both exit0,1.17s each.
Source baselines use FOUR fixed positive red-subtriangles with16x16 rules
(1024 points), a fixed reference-Duffy Cauchy radius7/8, and their own
analytic majorants. No source baseline is inferred from target rule
agreement or silently set to zero. Original-basis oracle excesses:0.

Each source's THREE paired edges use16 fixed reference parts. Difference
of the original local polynomials is restricted along the shared reference
line BEFORE magnitude bounds; the common initial-v0 factor is retained
with the actual coefficient difference. This avoids charging unrelated
off-edge variations as trace mismatch.

| Source | Largest original position trace gap bound | Largest velocity trace gap bound |
| --- | ---: | ---: |
| 225 | 3.25274e-19 | 1.92784e-17 |
| 137 | 3.25529e-19 | 2.78583e-17 |

For each part, the ruled-strip image area is bounded by gap times average
mapped tangent length, plus conservative endpoint caps. Whole-strip target
half-plane exclusions identify possible receivers. Both source field
envelopes include density, J, v, velocity gradient/Udot and potential.
Source inverse extensions across the strip are justified by explicit local
Lipschitz checks against(1+L_source)/2; the chosen extension radius then
covers the physical gap with the inverse condition factor.

**This is a finite error envelope for the unassigned gap/overlap region.**
It does not construct an exactly continuous global material map, assign
an owner to a gap, or average/refit traces. Exact global C0 remains false/
unproved. Earlier tangent-derivative gap diagnostics remain preserved; the
new field envelopes also bound full source velocity gradients in Udot.
Whether this envelope suffices for native promotion must be judged with
the whole-target sum and representation contract, not by ignoring it.

## Full local error sum and fixed work allocation

The already proved U lower bounds, physical tolerances and original
geometry-only possible-incidence counts are unchanged. Source baseline is
counted once in own-source complement; outgoing field/quadrature errors
are added without cancellation. Incoming error has no source baseline.
Gluing is charged directly to the possible target, once per touched shared
edge in this local audit. A whole-target assembly must canonicalize edges
to prevent charging the SAME seam again from both neighboring sources.

| Incoming row | Field + geometry + node + sum | Integral remainder | Shared-edge envelope | Total | Total / whole-target tolerance | Total / fixed work allocation |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 225 ->68526 potential | 1.31464e-25 | 9.18168e-24 | 7.32865e-23 | 8.25996e-23 | .0372772 | 2.38574 |
| 137 ->67942 potential | 4.11465e-26 | 2.92150e-24 | 8.01356e-23 | 8.30982e-23 | .0375212 | 2.17623 |

Whole tolerances are2.21582335687e-21 and2.21469752757e-21. Work shares
are1/64 and1/58. The remaining source contributions are uncomputed;
3.75% of a whole tolerance is not a whole-target pass. Conversely, exceeding
an equal work share is not proof that the full physical target fails.
No work share or U normalizer was changed to obtain a favorable result.

For all12 moments, own=source baseline-minus-outgoing is evaluated in the
actual binary80 operation order by exact decoded arithmetic. Source sum
residuals are within the calculated composition-rounding errors; the own
mass remains positive. This is per-source partition conservation, NOT
conservation of the unassembled24-cell state or a compression-feasibility
certificate.

## Preserved evidence and next review

- `CURVED_INVERSE_PHYSICAL_LEDGER_20261021.json`: complete selected local
  components, work/whole tolerance distinction, old/new continuous-integral
  comparisons and all native contracts.
- `CURVED_INVERSE_SOURCE_CONSERVATION_20261021.json`: exact decoded12-moment
  source baseline/outflow/complement conservation.
- `curved_target_integral_v2_CASE`: new physical positive packets/moments,
  domains, oriented boundaries, input/source/binary/packet SHA ledgers.
- `curved_target_analytic_v1_CASE`: per-sector complex certificates and
  analytic remainders, with old/new positive inputs preserved.
- `curved_target_source_glue_v2_CASE`: source baseline, finite trace/strip
  records, target field ambiguity envelopes and extension proofs.
- Original general-evaluator v1 passes and preliminary source/glue v1 are
  retained. Integral v1 was a BUILD-only mixed-auto-declaration failure;
  no simulation ran. Final source/glue v2 adds an explicit extension proof,
  without changing quadrature points, physical parameters or normalizers.
- All original material/positive packets and old nominal moments are
  unchanged and SHA checked. No simulation or restart was discarded.

This completes the selected lens ->225/137 local-error bundle. The next
substantive decision is a SINGLE complete target or the prospective shared
trace representation contract. Do not use local results as a substitute for
actual native CUDA/MPI long Kepler stability.
