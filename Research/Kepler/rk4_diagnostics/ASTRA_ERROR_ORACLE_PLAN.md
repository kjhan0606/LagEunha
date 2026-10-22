# Astra: original-function oracle, stable evaluation and target error ledger

2026-10-21. One bundled evaluation after CURVED_FOLD_RESULTS.md. Driver
implements/runs. The nominal geometry and positive target-U lower bounds
are accepted only for225,137 and the same tiny cell. Cold accuracy remains
unproved. Reuse their saved packets; no radial refinement or full384 remap.

## Physical reference and error paths

Treat the originally stored P3/residual coefficients as exact input
constants, with the ORIGINAL mathematical basis and reference density
defining chi_*,v_*. Do not replace this reference by the rounded compiled
coefficients. Each original element function is continuous; exact shared
trace/global C0 remains a separate unproved property and gluing error.

| Actual code path | Required error |
| --- | --- |
| regular label/boundary from compiled representative | original->compiled conversion, representative evaluation, plane input conversion, root solve |
| native_point integrand | native evaluation versus chi_*,v_*,F_*,V_*, propagated to every moment |
| event unknown owner from native position | possible native owners/domain shifts, retained in event bounds |

If regular labels never use native position, the full native position
envelope need not also move that boundary. This is an attribution change
based on actual dependencies, not deletion of error. Native position errors
remain in potential, angular momentum and position moments; event unknown
owner protection remains. Keep a table identifying any duplicate paths.

## Replace the broad ray arithmetic model

Preserve the old16384-operation bound as evidence. Implement/measure the
actual stable local-cubic+Bv0 evaluator and propagate operation-level error
through its real graph: Horner/FMA/exact differences as actually used.
Include cubic and B conversion, source coordinates, power normal/constant/
weights, value/derivative evaluation, root bracket and residual. No observed
maximum times a chosen multiplier as an asserted bound.

An independent fixed binary128/high-precision oracle evaluates the ORIGINAL
mathematical basis from the same stored coefficients. It catches incorrect
forward models; it is not a new state or an alternative geometry oracle.
No coefficient refit, trace averaging, state projection or heat correction.

## Finish four missing terms on stored positive points

- Source: separate source-reference arithmetic from quadrature sensitivity.
  Existing rule differences/finite comparisons are estimates, not rigorous
  integration upper bounds.
- Field: propagate coordinate, velocity, J, density, thermal, potential and
  velocity-gradient error, including density and weight in K. Keep Udot.
- Gluing: retain the whole-edge position gap and velocity/derivative gaps;
  separate domain ambiguity from field influence and possible receivers.
  Do not set a small gap to0 or assert exact C0.
- Sum: use actual order/precision and absolute term sums. Global assembly
  cancellation is not a reason to lower target-local error.

The saved nominal packets/moments must stay unchanged during attribution
work; any change needs separate tracking.

## Target budgets

Uniform1/(2*Ntri+1) was a work allocation, not the final physical accuracy
criterion. A certified possible source-target incidence can replace it:
only whole-source exclusion proofs remove a contributor, never observed
zero flux. Count own outgoing subtraction, incoming and baseline once.
Fix the incidence/allocation BEFORE seeing errors. The final gate is the
SUM of all possible target contributions <= tau_jk(U_j^-).

The maximum of already proved U lower bounds is valid; do not add overlapping
certificates. Passing a limited subset is partial acceptance only while
other source contributions are uncomputed.

## Fixed bundle / stopping conditions

1. Specify oracle/paths and audit arithmetic on the three stored packets.
2. Add stable-evaluator boundary/field forward bounds.
3. Finish source/gluing/sum and incidence-based partial target budgets.
4. Check preservation of original nominal packets/moments.

Stop on an oracle value outside the proposed model; on a dominant bound
report its required accuracy; leave gluing/source-integration uncertified
when it cannot be bounded. Only complete terms meeting the partial budget
support the next ONE whole-target verification. No precision/order/cap
scan or favorable normalizer change. Weight Newton/compression/evolution
remain after whole-target acceptance. Native CUDA/MPI t6/8/10 and longest
GIZMO duration, no steering/Re infinity/conditional W0 remain the full goal.
