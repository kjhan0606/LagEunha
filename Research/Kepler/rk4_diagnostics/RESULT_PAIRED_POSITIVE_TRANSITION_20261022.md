# Actual paired RK4 step 5 passes; unchanged-space long-run cost fails

2026-10-22. Full native Ada MPI4/OpenMP/CUDA execution, not a Python evolution
or a replay of an accepted interval. Long Kepler artificial heating remains
unsolved. Original user repository identity reconfirmed: LagEunha original
source HEAD `8e25dfe0c2e41a7374c11e2e3102da2a8b028625`.

## Execution and preserved state

Full positive preprocessing and the independent whole-record audit are
completed; see `RESULT_POSITIVE_WHOLE_RULE_20261022.md`. All 1,138,432 triangles
and 2,040,256,512 original points are accounted for, with 249,691,129 positive
selected points. The point count reduction of 8.1711x is not a step speedup.

H200 job 433999 was cancelled while pending, with authoritative Slurm runtime
zero and no start. Existing exclusive-claim dispatcher launched exactly one
Ada runner. Actual native PIDs 3028401--3028404 used the Ada GPU, about 5.49 GiB
per process at the observation; all four are now terminal/absent. Dispatcher
2351718 and runner 3026814 are also terminal/absent; dispatcher recorded exit 0.
The exclusive claim remains and must not be removed for a duplicate retry.

Authority: `material_full_positive_transition_v2_ada_mpi4/contract.json`,
`native_exit_code.txt` (0), and dispatcher `ada_runner_exit_code.txt` (0).
Final status: `ACTUAL_PAIRED_STEP_PASS_LONG_PRODUCTION_STILL_UNPROVED`.
Final integrity reports contain 68 source, 2 binary, 7 accepted-input,
6 prepared-k1 and 9 positive-rule records, all marked OK.

The same accepted trajectory advanced one genuine step:

| Item | Actual value |
|---|---:|
| Time before | 8.2355880701376864e-5 |
| Accepted h | 3.027789731668267e-5 |
| Accepted time | 0.00011263377801805954 |
| Total accepted steps | 5 |
| New accepted steps in this run | 1 |

Current authoritative checkpoint prefix:
`material_full_positive_transition_v2_ada_mpi4/run/continuous_step5/accepted_r`.
The original step4 checkpoint and its saved high k1 remain preserved. That k1
belongs to the step4 base and was reused for this transition's high comparison;
it must not be retargeted to the newly accepted step5 state.

## Actual accuracy and restart evidence

All four `stageN/quadrature_difference.json` gates and the endpoint gate pass
at the identical actual low-rule stage states. There is no separately evolved
high-rule trajectory. Stage1 explicitly records
`high_k1_force_and_solve_reused=true` and
`high_missing_source_summaries_only=true`; no saved high force/solve was replayed.
Endpoint comparison records no repeated high force or solve.

- Maximum source pressure difference over the four stages: 3.5667e-16.
- Maximum point-density relative bound: 1.9583e-16.
- Maximum pressure-rate difference: 1.7084e-13.
- Maximum mass-metric acceleration difference over velocity: 1.2948e-12.
- Cumulative numerical pressure/density/velocity error:
  3.4800e-16 / 1.9811e-16 / 3.9203e-17.
- This step's dt-scaled numerical allowance: 3.0277897316682677e-9.
- Endpoint delta U = -8.9512e-22; delta Phi = -3.4263e-17.

These are errors relative to the original rule at this very early state, not
errors relative to the physical Kepler solution at t6--10. Numerical integral
differences remain separate from physical heat; no physical state correction
or energy rebase was made.

Actual `step.json` passes its finite-step gates: E defect/U = -1.4095e-12,
L defect/(2K) = 1.2632e-17, thermal increment defect/U = -2.2671e-18,
authoritative source mass relative difference = 0. The point-vs-M kinetic
discrepancy/U = 1.6544e-12 is recorded, not silently erased.
`numerical_restart_roundtrip.json` version3 reloads the actual checkpoint:
coefficient hashes, physical and numerical ledgers, rule generation and active
mode all match, with no accepted interval replay.

Original-node corrected timestep bound at the endpoint is
0.0004322978132696750. Corrected sampled J lower bound is about0.999999998545;
this is not a global injectivity proof. Existing enriched map guard remains
incomplete. No shock/recipient thermal transfer/W0 claim follows from this step.

## Measured cost decision — both alternatives fail

Authority: `run/continuous_step5/selected_rule_cost.json`.

| Cost item | Measured or explicitly optimistic value |
|---|---:|
| Actual paired step wall | 492.388295 s |
| Separately timed one-time high-reference comparison | 361.841091 s |
| Selected step excluding that reference comparison | 130.547204 s |
| Optimistic sustainable h at current bound, with safety | 0.00034583825061574003 |
| Optimistic remaining wall to t1, excluding future guard | 377438.007 s = 4.368495 days |

The future h used in the projection is larger than the actual transition h;
it is not a measurement that long evolution can sustain that h. Future map
guard cost, mesh distortion and lifetime can worsen this estimate.

Previous same-Ada step4 total wall was442.167919s, including its first cache
comparison. Even treating that entire old wall as the baseline and excluding
all one-time reference work from the new step gives the favorable improvement
ratio only3.3870x, below10x. This is not an H200 throughput measurement.
The optimistic t1 estimate also exceeds24h. Therefore both adopted alternatives
in `ASTRA_LONG_EVOLUTION_COST_DECISION_20261022.md` fail. The native code stopped
at the preserved accepted checkpoint with
`measured_cost_exceeds_24h_t1_budget`.

**End further performance optimization of the unchanged 10,816,441-DOF space.**
Preserve the positive-rule, cache, resident-DD and native/RK4/restart work.
No blind24h continuation, new tolerance/order/dt/rank sweep, old interval replay
or additional primitive guard experiment is released by this result.

The next major bundle must address conservative field representation scale and
material shear lifetime, with explicit E/L/thermal and geometric-weight/mass
reduction contracts. One major Astra planning/evaluation request has been sent
using these actual results, as delegated by the user. Driver remains coder.
No other agent or external reviewer was started.

## Still required for the actual goal

t1 and t6/8/10 physical accuracy; pressure/density/orbit/rotation behavior;
global map validity or a compatible conservative representation transfer;
weight-dependent mass with approved W=0 Voronoi/fixed-mass reduction; shocks;
boundary-consistent comparison; and the longest correctly normalized GIZMO
paper Kepler evolution are all incomplete. Kepler remains Re=infinity, with
resolved rotation, no steering, no damping and no artificial thermal/spin
reservoir. Original user code/results and shared global context are preserved.
