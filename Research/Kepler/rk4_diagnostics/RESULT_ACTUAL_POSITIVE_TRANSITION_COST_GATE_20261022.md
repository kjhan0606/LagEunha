# Actual paired RK4 transition: accuracy passed, long-run cost failed

2026-10-22. The long-term Kepler objective remains active and unachieved.
The completed calculation used the existing native C/CUDA/MPI/OpenMP path
on one RTX5000 Ada, four MPI ranks and eight OpenMP threads per rank.

## Authoritative accepted state

The single downstream run completed with native exit 0 and actual runner
exit 0 recorded by the dispatcher. Its four native processes and dispatcher
are now terminal. H200 job 433999 had been cancelled while pending before
this sole Ada execution. No H200 performance or physical result exists.

The retained trajectory advanced once from accepted step4 at
`t=8.2355880701376864e-5` to accepted step5 at
**`t=0.00011263377801805954`**, with `h=3.027789731668267e-5`.
No accepted interval was replayed. Its checkpoint prefix is:

```
material_full_positive_transition_v2_ada_mpi4/run/continuous_step5/accepted_r
```

All four version3 shards were actually reloaded by the native reader, which
checked payload checksums and equality of coefficients, physical endpoint,
time/controller, active rule generation and physical/numerical ledgers.
The accepted first high-rule RHS force/solve was reused: stage1's
`high_k1_force_and_solve_reused=true` proves that fact. The separate step
field `cached_initial_k1_reused=false` means the low-rule rate was computed
for the paired transition; it does not mean the saved high-rule solve was
repeated.

## What accuracy evidence proves

All four actual same-state RK4 low/high comparisons, the endpoint comparison,
the original cold-energy/work/torque/mass gates, and restart readback passed.
The authoritative binary80 M, basis, generating entropy/context and initial
physical fields were retained. Rule-transition integral discrepancies stayed
in a separate numerical ledger; there was no heat or state correction.

| Quantity | Result |
| --- | ---: |
| Finite-step energy defect / U | -1.4094657634e-12 |
| Thermal increment defect / U | -2.2670779667e-18 |
| Angular defect / (2K) | +1.2632178720e-17 |
| Relative mass difference | 0 |
| Endpoint mass-metric vs point K error / U | 1.6544129928e-12 |
| Low-stage true M-solve residual, maximum | 8.1281904781e-19 |
| Endpoint minimum sampled J | 0.9999999985446771 |
| Endpoint accumulated extra pressure error | 3.4799907310e-16 |
| Endpoint accumulated extra density error | 1.9810618878e-16 |
| Endpoint accumulated extra velocity error | 3.9202357750e-17 |
| Additional-error allowance at this endpoint | 3.0277897317e-9 |

These are observations near the initial state, not a time-uniform quadrature
bound, global map-injectivity proof, or evidence that the historical t6–10
heating failure has been resolved. Total energy conservation alone would
also be insufficient to establish the latter.

## Measured cost rejects this fixed-space production strategy

The full rule originally used 2,040,256,512 positive points per pass. Its
completed positive-rule preprocessing took 18,966.2571 seconds (5.27 hours)
and retained 249,691,129 points, an 8.1711-fold point-count reduction.
Point count is not end-to-end speedup.

| Cost | Measured result |
| --- | ---: |
| Actual paired RK4 step including reference work | 492.388295 s |
| Selected-rule step, separately timed reference excluded | 130.547204 s |
| Improvement versus adopted 225 s baseline | 1.7235 times |
| Optimistic sustained h used for projection | 0.00034583825061574003 |
| Optimistic remaining t1 time | 377,438.007 s = 4.3685 days |
| Rank0 low-stage mass solve sum | 70.6411822 s |
| Rank0 low-stage preparation sum | 46.0810869 s |

The rank0 component sums are not independent whole-run timers and should
not be added to reference work or other overlapping timers. Resident halo
timings are cumulative counters; in particular stage4's halo number is not
the duration of that stage alone.

The projection already allows the current small controller step to grow to
0.8 times the current original-node dt bound. It excludes unmeasured future
guard overhead and worsening material-mesh deformation. It is therefore an
optimistic estimate, not a completion schedule.

Both adopted alternatives fail: the speed improvement is below tenfold and
t1 exceeds 24 hours. Native stop reason is
`measured_cost_exceeds_24h_t1_budget`. In accordance with
`ASTRA_LONG_EVOLUTION_COST_DECISION_20261022.md`, do not start a blind
continuation, another quadrature/precision/dt sweep, or further optimization
of this unchanged 10,816,441-scalar-DOF material space. Preserve the accepted
trajectory and useful GPU work. The next strategic decision must change the
representation/cost structure and address conservative physical-field and
thermal-context transfer, with an explicit route to the original t6–10 test.

## Retrospective audit

`audit_positive_transition_completion_20261022.py` completed successfully.
It wrote `material_full_positive_transition_completion_audit_v1/contract.json`
and its evidence manifest. It independently scanned all 71,080 low-stage
cell summaries at each of four stages and all 71,080 endpoint cell summaries,
checked finite/positive fields, ID coverage, and local energy/angular gates.
All four ranks completed all required low/high passes. It verified reported
final integrity counts (68 source, 2 binary, 7 accepted-input, 6 saved-k1,
9 positive-rule and 29 immutable-payload records), froze SHA256 of the actual
new accepted shards and checked their version/partition/ABI headers.

This audit does not independently recompute GPU point forces or repeat the
18 GB full-rule-file checks. Those checks were performed by the completed
native reader/wrapper and their records are preserved. The audit launches
no physical evolution and alters no checkpoint.

One bundled strategy review was sent to the existing Astra planning agent
after the measured cost failure. No new agent or test was started. The
optional compression-diagnostics build remains unexecuted and is not a
substitute for the next representation decision.

## Explicitly unfinished

Historical t6/8/10 stability and heating suppression; continuous/global map
validity; long material-space lifetime or conservative reexpression;
W0 dynamic-Voronoi/fixed-mass reduction; matched-paper IC/boundary and orbit
normalization; and the longest GIZMO-paper evolution are all unproved.
The preserved accepted step5 is a valid small-time state, not a revision of
the long-term success criterion.
