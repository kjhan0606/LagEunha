# Entropy-generated pressure and perturbation work at t10

At fixed evolved t10 geometry and fixed cell masses define
P_iso=P_initial*(V_initial/V)^gamma and P_entropy=P_actual-P_iso.
The full TCF reversible pressure operator is linear in pressure, so its
force splits exactly (max residual 5.93e-21). Reconstructed actual pressure
force agrees with saved hydro force minus rounded Q at relative L2 1.03e-7.
Initial-ID material bins are used.

| Region | Actual pressure work on v-v_circ | Isentropic part | Entropy-excess part |
|---|---:|---:|---:|
| inner background | 2.80264e-5 | 3.16280e-6 | 2.48636e-5 |
| inner disk edge | 1.30449e-6 | 1.33396e-7 | 1.17110e-6 |
| disk core | 5.72813e-10 | 4.76478e-11 | 5.25165e-10 |

Positive work supplies perturbation kinetic energy at this instant.
Entropy-excess pressure accounts for about 88.7% of inner pressure work
and 89.8% at the inner disk edge in this algebraic decomposition. Inner
pressure acceleration RMS is .30327 actual, .06117 isentropic and .27663
entropy-excess. Vector RMS components do not add.

This supports a heating/pressure/perturbation feedback pathway at t10 and
qualifies the earlier finding that geometric seeds alone are modest along
prescribed circular orbits. It is not a unique causal attribution of growth:
Q/HG forces can damp perturbations while their heating changes future pressure,
and pressure work also exchanges energy with the circular background. No
counterfactual trajectory with entropy removed was evolved. Energy source and
transport across a selected material subset require the complete budget.

Do not remove entropy or pressure after each step to stabilize the test;
that would change energy conservation and thermodynamics. A viable candidate
must reduce erroneous excitation while retaining shock entropy production
and pressure-work compatibility, then pass the original failure interval.

Evidence: canonical_pressure_entropy_force_2026-10-18.json.
Script: audit_pressure_entropy_force.py.
Remote result: Codex_qbp_t10_decomposition_20261018/pressure_entropy_force.json.
No trajectory was modified. Long-term physical accuracy remains unmet.
