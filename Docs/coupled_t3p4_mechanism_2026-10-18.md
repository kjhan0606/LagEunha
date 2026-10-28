# Coupled response through t3.4004059

Existing evolution was preserved. Completed snapshots were analyzed on
lageunha; no production force formula was changed.

## Matched-time gravity controls

Squared circular-departure response relative to the original t2 seed:

| Time | Hydro | Gravity only | Hydro/gravity |
|---|---:|---:|---:|
| 3.0004059 | 16.45961 | .932439 | 17.65221 |
| 3.2004059 | 30.66060 | .997576 | 30.73509 |
| 3.4004059 | 79.62488 | 1.028640 | 77.40794 |

Gravity-only RK4 caps .0005/.00025 change the diagnostic growth ratios by
at most1.05e-9 over these checkpoints. The large hydro departure response
is not reproduced by independent orbital evolution from the same seed.
This comparison is not an additive causal energy budget or eigenmode proof.
Hydro amplitude/dt convergence at these later endpoints remains pending.

## Localization and reversible pressure contribution

At t3.4, initial-radius material bins .3-.5 contain95.22794% of departure
response and84.10116% of acoustic response. The top10 cells contain30.12%
and37.75%, respectively. The center .1-.2 now has departure squared growth
1.38082, whereas0-.1 remains damped(.02272). Thus the earlier statement
that all r<.2 bins damp must remain scoped to t2.8.

Inner r0<.5 instantaneous pressure contribution sum(delta_u dot delta_F):

| Algebraic component | Power |
|---|---:|
| Total reversible pressure | +3.36941784e-7 |
| Pressure-value component | +2.27435185e-7 |
| Explicit geometry-operator component | +1.09506599e-7 |
| Isentropic with each cell's seed-t2 entropy | +3.29916974e-7 |
| Entropy excess since seed | +7.02481032e-9 |

The explicit geometry component is now32.5% of pressure power, compared
with2.62% at t2.8; the seed-entropy isentropic component is97.9% of total.
These are two distinct algebraic decompositions, not independent causes.
Pressure values themselves depend on geometry. The seed already includes
pre-t2 entropy history; this is not a claim that all AV heating is negligible.
TCF/HG, gravity, and the position derivative of circular reference velocity
are absent from this pressure-only budget. Positive power does not prove
negative stiffness. No stabilizing modification is justified by this alone.

Independent Voronoi volumes match both snapshots within3.464e-14 for18385
cells at current r<2.4; IDs/masses match the seed and all weights are zero.
Algebraic force split closures are below2.4e-12. This validates the geometry
used by this diagnostic, not the complete force or backend implementation.

## Next decision

Finish the running fine-dt comparison at t2.6 and retain the original t4
continuation. The growing geometric contribution motivates examining the
full coupled pressure/geometry response, rather than an AV-only remedy.
Do not replace historical t6-10 failure and long-time physical gates with
these short mechanism diagnostics. CUDA hydro, long MPI physics, and
time-dependent Laguerre weights remain outstanding parts of the objective.

Evidence: canonical_coupled_t3p4_comparison_2026-10-18.json,
canonical_coupled_t3p4_regions_2026-10-18.json,
canonical_coupled_t3p4_pressure_2026-10-18.json,
canonical_coupled_t3p4_geometry_2026-10-18.json.

## Actual hydro RHS and HG generator added

Two terminal MPI1/OMP4 production probes from the saved +/- checkpoint743
completed with exit0 and no accepted RK step or HG application. Snapshot
hashes match, initialization x/v/m/P/V match, and each emits221081 HG modes.
Mode reconstruction uses the existing production nonlinear HG coefficient
law, including frozen-support treatment. Its kinetic/thermal balance closes
within2.64e-22. The analytic circular-reference Jacobian agrees with centered
directional finite differences to2.70e-10 at the smaller increment.

Inner r0<.5 instantaneous contributions to the derivative of half the
mass-weighted squared central departure response:

| Term | Power |
|---|---:|
| Gravity | +1.06847356e-6 |
| Moving circular-reference velocity | -9.66742974e-7 |
| Their sum | +1.01730589e-7 |
| Reversible pressure | +3.36941784e-7 |
| Production hydro RHS minus reversible pressure | -7.93245283e-8 |
| HG modal generator | -4.37110299e-8 |
| Total instantaneous generator | +3.15636815e-7 |

The hydro remainder includes TCF and any other hydro RHS terms; it is not
declared pure TCF without a term-by-term identity. Both this remainder and
HG damp the inner response at this instant. Their signs differ at the edge
(.5-.6), so no global universal damping claim follows. The large gravity
term must not be quoted without its cancelling circular-reference term.

This is a continuous-generator diagnostic. A finite ordered HG split map
need not equal its zero-step generator; an accepted-step response-budget
closure is still needed. Positive departure power is not creation of total
energy, nor proof of an unstable eigenmode or defective pressure Hessian.
Evidence: canonical_coupled_t3p4_kinematic_2026-10-18.json and
canonical_coupled_t3p4_generator_2026-10-18.json.

## One-step consistency check: discrepancy remains

Four actual production HG/RK4/HG restarts from checkpoint743 completed:
plus/minus at dt1e-4 and5e-5, each exactly one finite accepted step, no
retries, exit0. Existing original and refined continuations were untouched.
The forward change of half the squared central departure response gives:

| Fixed initial region | Generator | Measured dt1e-4 | Measured dt5e-5 |
|---|---:|---:|---:|
| r<.5 | 3.15636815e-7 | 3.13884875e-7 | 3.13787542e-7 |
| .3-.5 | 3.03275995e-7 | 3.03136799e-7 | 3.03697435e-7 |
| r<.2 | 3.33587737e-10 | 3.59196135e-10 | 3.59360333e-10 |
| .5-.6 | 2.62660035e-10 | 2.61766954e-10 | 3.01150612e-10 |

Inner mismatch is -.5550% then -.5859%; central mismatch is about+7.7%.
The edge difference increases to14.65%. These do not show the expected
simple forward-derivative convergence. This is an unresolved consistency
finding, NOT a passed generator/finite-map closure. Growth stays positive,
but precise force attribution requires investigating this discrepancy.

Next isolate per-branch acceleration errors and central differencing
sensitivity before attributing the mismatch to HG, topology changes, or
roundoff. Inspection of hg_apply_mode confirms the reconstructed zero-step
coefficient law algebraically; that alone does not validate its agreement
with the assembled finite-step map. Do not change physical terms based on
this unclosed diagnostic. Evidence:
canonical_coupled_generator_step_2026-10-18.json.
Runner/analyzer: check_coupled_generator_step.py.

## Acceleration localization and OpenMP control

Using the same saved outputs, each branch's inner total-acceleration forward
difference error decreases from0.02927% to0.01464% when dt halves. This is
consistent with forward differencing of the dominant smooth orbital force;
it does NOT establish convergence of the much smaller hydro response.
The central-perturbation acceleration relative L2 error instead changes
from3.513% to4.361%. Its inner projected discrepancy remains approximately
-1.95e-9. Cells154 and156 at initial radius.21875 retain conspicuous
negative contributions at both timesteps. Other leading cells change:
507/487 at dt1e-4,312/734 at dt5e-5; the edge discrepancy at the smaller
step is concentrated in848,849,863. This localization does not establish
a topology-change or rounding cause.

Two independent OMP1 one-step runs at dt5e-5 completed with exit0 and
reproduce OMP4 bitwise in every reader field: id,x,y,vx,vy,mass,volume,den,
pressure,ie,w2, for both +/- branches. Full binary hashes differ, so do not
claim byte-identical snapshot files. This rules out thread-count dependence
as the explanation of this particular observed discrepancy; it is not a
general race-free proof. No MPI or production CUDA-force scope is implied.

Next inspect actual stage force and geometry updates at the identified
cells, comparing the terminal probe with the first actual split/RK stages.
The zero-step HG coefficient formula matches the inspected hg_apply_mode
algebra, but assembling/staging its fields still needs verification.
Evidence: canonical_coupled_acceleration_discrepancy_2026-10-18.json,
canonical_coupled_generator_threads_2026-10-18.json.

## Actual stages locate the discrepancy in changing RK force evaluations

Two OMP1 dt5e-5 runs dumped the six actual force calls: pre-HG assembly,
four RK stages, final HG assembly. Both exited0 with one accepted step.
Every reader field matches the no-dump OMP1 controls bitwise. Initial hydro
force equals the terminal probe exactly in all audited interior bins.
Actual RK velocity change matches the weighted four-stage hydro+gravity
RHS to mass-weighted absolute L2 about3.01e-13 inner (1.69e-12 edge).
Thus no incorrect RK stage weights are demonstrated here.

Inner response discrepancy power decomposes into first HG+5.53e-13,
last HG-7.69e-12, and RK relative to initial hydro+gravity-1.94495e-9.
The dominant difference is the evolving force sampled by RK, not HG's
initial-generator reconstruction or the RK weighted-sum implementation.

For the minus branch, cell154's first-to-second RK force changes by
(1.51710e-9,1.66749e-10). Reconstructing just -sum(A*Q*e) from rounded
face diagnostics gives(1.54384e-9,1.59142e-10), leaving the small remainder
(-2.67460e-11,7.60680e-12). The same localization appears at156 and848.
At154, the six neighbor IDs remain unchanged. On face154-155, Q changes
from3.214168e-7 to2.582168e-7; QBP delta changes from.0009372032 to
.0007508496, and chi from.002359157 to.005087549. Cell155's reconstructed
velocity gradient changes noticeably while its face-neighbor set also
remains fixed. The plus branch does not show the analogous large jump.

This identifies a Q/reconstruction transition behind the local diagnostic
discrepancy, not the root cause of the full long-time physical failure.
Unchanged face neighbors do not prove the larger reconstruction stencil
is unchanged. Rounded face records limit exact force closure. Next inspect
gradient stencil selection/rank and QBP switching, especially cell155.
Do not confuse net AV damping of the response with smooth dependence of
its local force on the state, or remove AV based on this finding alone.

Evidence: canonical_coupled_stage_decomposition_2026-10-18.json,
canonical_coupled_selected_stages_2026-10-18.json,
canonical_coupled_face_transition_2026-10-18.json.
