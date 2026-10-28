# Accepted operator budget over t2 to t2+0.1

Reused the completed onset reference and +/- velocity-seed branches.
No simulations were rerun. All 40 steps per branch were accepted without
retries; all 80 split HG applications have four radial-bin records each.
Summing all bins avoids material-membership ambiguity.

For the reference, circular-departure kinetic energy changes by:

| Contribution | Delta |
|---|---:|
| Full trajectory | +1.1516647161e-10 |
| HG, directly measured finite increments | -2.9562488392e-11 |
| RK, endpoint change minus HG | +1.4472896000e-10 |

HG simultaneously adds 8.0072170361e-11 of thermal energy. Its thermal
increment is not the negative of circular-departure energy: those are
different budgets. HG finite-increment identity absolute residual sums
to 1.99e-23. Both perturbed branches retain these contribution signs.

Thus HG directly damps the global velocity departure over this actual
early interval, while the combined RK evolution raises it. This does not
exclude indirect HG-heating feedback through subsequent pressure forces.
The logs cannot separate the RK increment into pressure, TCF and shear.
That requires actual stage attribution or a resolved trajectory quadrature.

Four particles change membership in r<0.5 between endpoints. The HG log
uses instantaneous radius; endpoint material differences use start IDs.
Do not subtract these regional quantities and call the result an exact
inner RK budget. Even unchanged endpoint membership would not by itself
exclude intermediate crossings. Global subtraction is unaffected.

Next: implement/passivity-check an accepted-stage force-work diagnostic
with a consistent material mask (or global sum), including the circular
reference's position derivative. Account for finite RK chain-rule error
explicitly; do not label a stage quadrature as an exact finite-step split.
Use the existing endpoint checkpoint for prospective measurements.

Evidence: canonical_onset_operator_budget_2026-10-18.json.
Analyzer: analyze_onset_operator_budget.py.
Physical long-term stability remains unresolved.
