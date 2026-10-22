# Curved transfer: nominal geometry complete, cold error budget still fails

2026-10-21. This goal turn is PROGRESS. Original physical map/velocity/IC
coefficients, Re=infinity, no steering, and native CUDA/MPI completion goal
are unchanged. Controls run as native C++/MPI1/OMP1 on authorized lageunha.
These are the SAME finite h=2e-4 geometry controls, not time evolution.

## Implemented geometry

1. Reference-domain predicates now use oriented edge equations. Local
   TwoSum/FMA checks retain exactly representable binary edge arithmetic;
   non-exact operations still return outward intervals. No global interval
   arithmetic rule or physical map coefficient was changed. Stable exact
   circle/axis-line intersections avoid inverse barycentric normals and
   subtraction of almost equal squared radii.
2. Whole-panel simple-root tubes remain available. Source-clipped monotone
   branches use derivative exclusion on a fixed32-box inclusion grid plus
   the exhaustive source-filtered event guard, rather than an exterior fold
   being invented as an interior fluid event.
3. Simple-fold families require strict Krawczyk inclusion and separated f_r,
   f_theta_theta. One validation box is derived a priori from the center
   residual uncertainty/inverse derivative, not a box-width scan. Critical
   points near a chart edge are certified in a neighborhood extending past
   u=+/-1, while integration charts still partition exactly[-1,1]. A fixed
   angular/radial inclusion grid bounds the whole neighborhood/path.
4. For each confirmed family, |f_theta|>=C*sqrt(r-H), with H an UPPER bound
   on the complete protected critical-radius set, including arithmetic
   uncertainty. The entire box remains an event; H is not called an exact
   tangent. A shifted square coordinate regularizes adjacent integration.
   The reciprocal derivative is integrated analytically:

       integral_a^b r/sqrt(r-H) dr
       = 2*[H*(sqrt(b-H)-sqrt(a-H))
            +((b-H)^(3/2)-(a-H)^(3/2))/3].

   This also tightens uncertainty when a rectangle derivative is positive
   but very small at an endpoint. No denominator floor is introduced.
   Exact linear reference edges have a corresponding source-fold bound.
5. Within event boxes, whole-radius polynomial root bands cover every
   unresolved angular boundary. Integrate the certified known gaps with
   the original positive8x12/16x16 rules; retain unknown positive points
   and bound true unknown contribution PLUS nominal unknown contribution.
   Known-gap quadrature remains in the estimator. Whole annuli are not
   silently dropped or treated as exact. Event possible owners are checked
   across whole bands, including an assertion for every nominal unknown
   owner. Source-boundary receivers use eight certified thickened edge
   segments; power-boundary errors affect the two sites defining the plane.

## Completed final v5 nominal controls

| Same control | Panels / evaluations | Original cap | Known-region ratio | Unresolved / trace / order / fallback | Native time |
| --- | ---: | ---: | ---: | ---: | ---: |
| triangle225 | 142 /142 | 268 | .0038402727 | 0 /0 /0 /0 | 8.12s |
| triangle137 | 138 /138 | 266 | .0018763719 | 0 /0 /0 /0 | 13.18s |
| original tiny cell | 54 /54 | 300 | .7126285487 | 0 /0 /0 /0 | 2.24s |

Folders `curved_fold_transfer_v5_{225,137,lens}` all exit0 with
NOMINAL_GEOMETRY_COMPLETE_NOT_COLD_ACCEPTANCE. No adaptive subdivision was
needed in these final executions. The implemented adaptive driver retains
the ORIGINAL initial-event-cut count+256 evaluation cap and depth24, not
the number of previously accepted leaves+256. Every evaluated panel is
retained in evaluated_panels.csv; a boundary-budget-dominated state stops
radial refinement. No independent bisection witness was repeated.

The original tiny-cell area5.196152422706632e-8 is retained with error
1.7041213e-23. Native velocity-gradient enclosure excess is0 in all cases.
All original tolerances and rule orders remain. A Gauss discrepancy is an
estimate, not a rigorous quadrature-error theorem.

The positive packets contain37567,34498,32289 records respectively; every
record is finite with strictly positive reference weight/density and a valid
target/flag. Their total104354 records reproduce target outgoing mass to
9.24e-19/1.42e-18/3.80e-18 relative using70-digit reconstruction. The earlier
float-only audit rounded these differences to0; that was not exact equality.
`CURVED_FOLD_PACKET_ROUNDTRIP_V2_20261021.json` is the corrected comparison.

Receiver-budget refinement v4->v5 leaves all three positive-point packets
BIT IDENTICAL. It changes uncertainty attribution only. v4 packet SHA:
2250ef8b973...,1371c8f84f6...,lensa5cd4263.... Full SHA/audit evidence is
in CURVED_FOLD_TRANSFER_AUDIT_20261021.json and the immutable run folders.

## Failed intermediate controls remain

- v1_225: geometry/nominal pass, but a whole-event-annulus uncertainty bound
  was far too loose (source-budget ratio17840). Full-angle root covers
  replaced that bound, retaining known-gap quadrature and unknown mass.
- v2_225/137 passed nominal controls; v2_lens failed with8 unresolved roots.
- v3_lens resolves all roots but fails nominal/error gates. Using a minimum
  endpoint derivative over an entire panel inflated its native boundary
  bound to2.51e10. The integrated fold bound replaces that overestimate.
- No successful/failure status was rewritten. No Kepler job was stopped,
  reset, or relabelled as a long-time success.

## Positive target-U normalizers

`test_curved_target_lower_v1.cpp` certifies target-owned positive interior
reference triangles. J>0 and ||F-I||<1 give local injectivity. Disjoint
mapped bounding boxes permit adding their thermal lower bounds. Density
floor and pressure are READ from the original control, never imposed on a
state. The physical integrand is the preserved P0*J^(1-gamma)/(gamma-1).
No E-K subtraction or thermal correction occurs. The powl<=8ulp assumption
is explicit and the lower evaluation is widened downward.

- v1_225:3targets,48source boxes; v1_137:4targets,64boxes. All positive;
  minimum U_lower2.2146975276e-11.
- v1_lens: all4positive,175boxes; tiny target101 has lower6.9848668650e-13.
- v2_lens: largest-first search uses the SAME256-box/target cap. It proves
  stronger lower bounds for the three large targets but FAILS to certify
  tiny101 within that budget (native exit2). Its0 means no bound found,
  not zero physical thermal energy. Preserve this failure.
- Use the MAX of already valid per-target lower bounds from v1/v2; do not
  sum overlapping certificates or rerun to reproduce a missing bound.
  Region areas and pairwise mapped-box disjointness are independently
  reconstructed using the printed binary80-roundtrip coordinates.

## The cold gate remains unproved

v5 source-normalized maximum representation/native-boundary ratios:
225562.01/306.25;137128.92/294.66;lens1.034/113.65. Event ratios are below1.
These are UPPER BOUNDS exceeding a budget, not measured physical errors.

`audit_curved_target_budget_v1.py` and
`CURVED_TARGET_PARTIAL_BUDGET_20261021.json` attribute each identical
overlap record to source outgoing and receiver incoming. The source-own
error conservatively sums its outgoing errors; no global cancellation is
used. A fixed1/(2*N_triangle+1) share reserves target cold budget for every
possible source contribution. Positive target-U normalizers give maximum
partial known-bound/budget ratios3.88e5,1.90e5,2.02e7. Uniform shares and
interior-region lower bounds are conservative. These ratios do not prove
the actual solution has such an error; current bounds cannot certify it.

Missing terms are explicitly NULL/incomplete: source reference/quadrature,
field evaluation, analytic trace gluing, accumulation rounding. No full
target/cold acceptance, mass-matrix positivity, weight Newton, compression,
full384 remap or Kepler time evolution is claimed.

## Next planning question

Complete and tighten the error attribution/forward models. In particular,
regular target labels come from the compiled representative; native_point
is used for integrand fields. Clarify whether native position evaluation
uncertainty belongs in regular boundary motion, field evaluation, or both
before improving the numerical bounds. Event unknown-owner evaluation
already has its explicit native envelope. Do not silently delete an error
term or change the mathematical physical-map definition.

The pre-cancellation16384-operation ray envelope is much wider than observed
root-reference shifts. A stable evaluator with its actual forward error may
be needed. Exact analytic trace continuity remains unproved; its saved
nonzero gap cannot be discarded. No precision/order/cap scan, refit of a
physical state, heat addition or new PDE is authorized by these nominal
successes. The final native CUDA/MPI t6/8/10 and longest-GIZMO goal is ACTIVE.
