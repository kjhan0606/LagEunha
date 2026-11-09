# Cylinder local entropy and volume consistency (2026-11-08)

The current priority remains the first astronomy draft, with Cylinder first
among unfinished 2D tests. This diagnostic retains evidence for the later
fluid-mechanics paper; no second manuscript is started. No new physical law,
coefficient, energy repair or canonical-method change is introduced here.

Let S denote Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/
cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1,
O=S/contact_failure_t2838_v1/positive_gap_v1, and A=O/entropy_trace_v1.

## Implementation and identity

A/native adds read-only per-cell Udot+P Vdot and entropy-rate logging to the
existing native C/CUDA adapter. Only libexam.a:exam_gpu_extract.o changes;
every other archive member, the physical update and timestep selection are
preserved. The existing contact/stage and volume samples are reused.
A/parity_exact.json proves all 13 extracted physical fields of all126696
cells are exactly equal to the prior same-source/same-cap t2.839 calculation.
Existing selected_step_budget.py is used without modification, with its hash
recorded. The effective gamma is1.399999976158142 (the float representation of
the input1.4), verified against actual RHS P. The preceding nominal-gamma
entropy audit has a negligible numerical difference and the same conclusion.

## Measured local mechanism

For retained ID63222 during t2.837825--2.839, all sampled first/predictor
RHS values H=Udot+P Vdot are positive (minimum0.00181839414). Nevertheless
the measured logarithmic entropy proxy change is negative:

    Delta ln(P/rho^gamma) = +0.02972412612 -0.00000000227 -0.04943228268
                          = -0.01970815884.

The three terms are the same-stage U/V quadrature proxy, the discrepancy of
actual U from its quadrature, and (gamma-1)ln(Vactual/Vquadrature).
This is an exact algebraic identity for the recorded states. It establishes
that volume inconsistency dominates this window's local entropy loss; it
does not identify the origin of the entire inherited65% loss or establish
positive entropy production in every other cell. No direct contact event
occurs for ID63222 in this interval.

At the largest residual step, t2.8378701039288883--2.837895242568985, both
sampled Vdot values are positive, yet the rebuilt volume decreases.
The defect is -1.399929578e-6 (4.40344% of this cell's initial step volume).
Other15 material IDs undergo contact in this step, so absence of a contact
for the selected cell does not rule out neighbour/contact contributions.

## Actual smaller-step comparison

All three cases use the identical saved state at2.837825 and end at2.837925.
The caps are below the previous problem step's actual dt. All126696 retained
masses are exact; U and V remain positive and every run exits0.

| dt cap | Accepted steps | Maximum local volume-step discrepancy | Selected net entropy-proxy change |
|---|---:|---:|---:|
| 1.25e-5 | 20 | 1.42844% | -0.87716% |
| 6.25e-6 | 26 | 0.462843% | -0.30455% |
| 3.125e-6 | 39 | 0.248512% | +0.07385% |

For the finest case the interval's net decrease disappears, but individual
steps still have negative entropy changes. The geometric log discrepancy
decreases in magnitude0.0115417 ->0.00570748 ->0.00186383; actual thermal-U
discrepancies are below4.7e-10. This supports a finite-step geometric/contact
issue. Adaptive substeps and different peak locations remain, so this is not
a formal uniform-step convergence-order measurement or a completed fix.
The last two endpoint densities of ID63222 still differ0.942%; whole-domain
mass-weighted relativeL1 is only1.32e-6 and masks that local error.

Reducing every timestep to this size is not adopted as the practical solution.
The speed objective and accurate wall geometry/pressure work must be met
together. Existing endpoint/stage volume-coordinate inverse implementations
remain available as diagnostic controls, with their earlier cost and boundary
accuracy limitations; they are not newly invented or already qualified here.

## Evidence and state

A/assessment/completed.json and per-run JSON files contain the exact stage
decomposition, matched RHS records and endpoint comparisons. Source hashes,
input lineage, effective params, output plan and process result are retained
for every run. Shared-GPU times are not a matched-cost benchmark.

The additional t5 continuation is terminal by intentional SIGTERM at its last
accepted t2.8708767579956156; the saved restart is2.839. Old t5 results and all
checkpoints remain preserved. The active work is short root-cause diagnosis,
not another long provisional trajectory, 3D evolution or manuscript2.

## Neighbour-contact sensitivity at the exact maximum-error state

The recorded 16-step prefix is reproduced exactly through2.837895242568985;
accepted t/dt/energy records, selected stages and volume budgets match the
original longer parity run. This state is saved separately. The existing
jump_volume() function is reused in the current gap-aware native source,
without the obsolete standalone reader or its rounded-radius rejection.
A frozen first-RHS reader reproduces every saved owned physical field exactly
and exits before any evolution step.

Although ID63222 has no own contact, all15 other contact IDs and their radial
image/periodic dependencies give endpoint first-order volume sensitivity
-1.393151106e-6, versus actual volume defect-1.399929578e-6: a ratio0.995158.
The second-stage contact contribution dominates; the weighted first-stage
contribution is only-6.32345e-11. This strongly localizes this event to the
coupling of neighbour contact position changes and rebuilt cell volumes.
This first-order attribution is not an exact finite-jump decomposition.
Its independent finite-difference check is recorded below. It does not prove
that every inherited thermodynamic error has this origin.

Evidence: A/contact_attribution_v1/replay_completed.json and completed.json,
with the exact endpoint checkpoint, contact displacements, native source,
build hashes and frozen process result. All new diagnostic runs have finished.
The subsequent sections record the independent derivative check and the
existing inverse comparison. Both preserve fixed masses, the original wall
law and energy accounting. Do not repair U
post hoc or promote an expensive capacity inverse without cost evidence.

## Independent native-geometry derivative check

The actual native checkpoint codec now creates private position-only clones
along the measured contact displacement. At scale0, all three checkpoint
files and the frozen owned state are byte-exact. For nonzero scales, all
particle bytes outside x/y, the footer and exterior state remain unchanged;
only the required checksums change. Native geometry refresh supplies the new
V; mass, U, W, velocity and identities remain unchanged.

Forward differences with scales0.1,0.05,0.025 differ from the analytic volume
sensitivity for ID63222 by0.0384865%,0.0192564%,0.00969625%, respectively.
All perturbed actual GPU geometries pass the existing positivity/support
checks. This independently qualifies the local derivative at this nearly
coincident contact event. It does not turn the first-order attribution into
an exact decomposition of the finite contact jump or a global method result.
Evidence: A/contact_attribution_v1/finite_difference_v1/completed.json.

A direct thermal correction -P deltaV is not a compatible remedy: nonuniform
P means its sum need not vanish even when the volume changes sum to zero.
The comparison below reuses the existing volume-coordinate inverse, solving
for W so that geometry realizes the volumes integrated by the same pressure
work, while preserving positions, momenta, U and the original contact ledger.
It is the previously implemented capacity method, not a new model; its cost
and remaining wall/entropy limitations must be measured before any adoption.

## Existing volume inverse: terminal comparison and storage diagnosis

The existing volume-coordinate inverse is now reconnected to the current
positive-gap native geometry. The opt-in main object is reused from the earlier
qualified inverse build; only exam.o and exam_gpu_extract.o differ from the
entropy-probe parent archive. With the inverse disabled, the 16-step control
reproduces all 13 extracted physical fields of all 126696 cells exactly.
This qualifies the adapter, not the enabled inverse's physical accuracy.

Both enabled comparisons start from the saved t=2.837825 state and terminate
with native exit 99 before accepting a step: endpoint mode in 82.076 s and
predictor-plus-endpoint mode in 77.444 s. They perform 108 and 110 nonlinear
volume evaluations respectively. Reducing the trial step does not make the
strict relative-volume threshold 1e-11 attainable by the recorded iteration.
These shared-GPU times are failure costs, not production performance results.
They do not reach the later t=2.8378952426 maximum-contact-error event.

Read-only native diagnostics replay the first two inverse updates and stop
before the third; no fluid step is accepted. The first three volume records
are exactly equal to the failed native histories. At material ID63222,
W is about 0.0132154 and its FP64 spacing is 1.734723476e-18. The difference
between the requested Newton increment and the actually stored increment
is amplified by the near-contact volume Jacobian.

| Mode and correction | Predicted next signed relative residual from stored W | Actual next residual magnitude | Magnitude difference |
|---|---:|---:|---:|
| Endpoint iteration 1 | -3.476911027e-9 | 3.476910994e-9 | 3.34e-17 |
| Endpoint iteration 2 | +1.753299398e-9 | 1.753299475e-9 | 7.70e-17 |
| Predictor iteration 2 | +2.374696333e-9 | 2.374696598e-9 | 2.65e-16 |

The nominal linear correction has selected-row relative residual about 1e-24.
Using the stored update instead predicts the measured stagnation magnitudes.
The later signed predictions are not inferred from absolute-residual logs. The row-wise half-ULP rounding
bound is about 5.23e-9. A 100-digit Decimal evaluation of that row independently
confirms the prediction. The preceding, larger predictor correction retains
a nonlinear remainder about 4.16e-11; not every iteration is purely linear.
This diagnoses the observed stopping mechanism. It is not a proof that no
other globally representable W vector could satisfy the strict threshold.

The initial inline diagnostic expression (W+delta)-W returned delta under the
optimized build and must not be used as a stored-increment measurement.
The authoritative analysis instead subtracts independently dumped binary64
next-W and W; next-W matches the subsequent native W for every cell exactly.
Initial assessment.json files are retained but superseded explicitly by
assessment_stored.json and independent_receipt.json. No numerical conclusion
relies on the invalid inline subtraction.

Evidence: A/volume_coordinate_v1/control_exact.json, endpoint and stages
process_result.json, and representability_v1/{completed.json,
independent_receipt.json,*/assessment_stored.json}. All these runs are terminal.
The original source checkpoints, wall law and boundary-energy ledger remain
unchanged. No canonical physics update, thermal repair, tolerance relaxation,
long trajectory, 3D evolution or second manuscript has been introduced.

The next decision must address the precision of the geometric coordinate or
use a justified finite-precision convergence criterion, rather than forcing
smaller timesteps through a storage limit. Any acceptance criterion must keep
the residual and its accumulated thermodynamic effect visible and must never
hide the earlier percent-level volume defect. A gauge change confined to owned
weights is not admissible with fixed exterior weights. Cost and the still
unqualified wall thermodynamics remain required evidence before adoption.

## Subsequent completed capacity comparison

The bounded precision classification now permits the actual existing inverse
to evolve from2.837825 through2.839. Both endpoint and stage variants complete35
steps with fixed masses. The selected entropy loss is removed in this interval,
while every accepted-step volume discrepancy remains below3.48e-9. Costs,
limitations, intermediate-stage interpretation and exact evidence are recorded
in Docs/CYLINDER_CAPACITY_PRECISION_2026-11-08.md. This supersedes the previous
"zero accepted steps" status only for the new explicitly classified candidate;
the original strict inverse failures remain preserved.
