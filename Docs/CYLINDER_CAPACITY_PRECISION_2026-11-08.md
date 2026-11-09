# Cylinder capacity inverse and finite-precision stopping (2026-11-08)

The first astronomy draft remains the priority. This is evidence from existing
Cylinder code for that work and for a later fluid-mechanics paper, not a second
manuscript or a claim of complete 2D/3D qualification.

Let S denote Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/
cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1,
A=S/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1, and
Q=A/volume_coordinate_v1/precision_stop_v1.

## Existing equations and new stopping classification

The existing inverse varies W to realize the volume integrated by the existing
pressure-work rates. It does not alter m, x, momentum, U, the contact law or
the kinetic-energy loss ledger. All native C/CUDA force and geometry objects
are reused. The precision-stop build changes only exam_gpu_extract.o in the
parent archive. No canonical physical method has been adopted.

The original relative-volume threshold 1e-11 remains normal convergence.
A distinct precision-limited stop requires a previous actual Newton update:

    r_pred = r_previous - J_previous * deltaW_stored.

Here r=Vtarget-Vactual. For every material ID, the next stored W and target
must equal the saved values, and stage, time, dt and iteration must match.
The signed actual-minus-predicted residual divided by Vtarget must be <=1e-11.
The row-wise nominal linear residual and the triangle bound for the actual
storage-rounding error are also checked. Both actual volume residual and the
half-ULP sensitivity bound must remain <=1e-8. This last bound was declared
before the tests as a diagnostic error budget, not fitted to each outcome.
It is not a physical coefficient or a substitute for accuracy convergence.

The half-ULP bound uses next-W spacing, including an exponent-boundary test;
self-image terms cancel and prescribed-exterior update error is zero.
Precise floating-point helpers prevent compiler reassociation of the stored
increment. New trial/stage and completion/failure invalidate predictions.
The standalone helper now also rejects an infinite target; the native caller
already checked it. The usual return status preserves the final device geometry.

Saved-state checks cover every row, cell reordering, stale W/target/stage/dt/
iteration, trial reset, excessive error, tighter budget, and nonlinear predictor
remainders. Disabled control reproduces all 13 physical fields exactly.
The subsequent precision-limited evolutions are separate from those checks.

## Actual evolution, same saved source and final time

Both modes start from t=2.837825 with 126696 retained cells, using the current
original-wall C/CUDA executable and dt cap6.25e-5. They reach2.837925 in17
accepted steps, then restart from those saved endpoints and reach2.839 in18
more. Selected restart state is exact in x,y,vx,vy,m,U,V,W. Both native
processes exit0; no earlier trajectory is restarted from its IC.

| Method | Steps to2.839 | Selected maximum step volume discrepancy | Selected entropy-proxy change |
|---|---:|---:|---:|
| No capacity inverse, previous same-cap run | 34 | 4.40344e-2 | -1.95152% |
| Endpoint inverse | 35 | 3.47708e-9 | +2.23556% |
| Predictor + endpoint inverse | 35 | 2.63914e-9 | +2.23347% |

For the two inverse runs, the listed maxima also bound all retained cells in
the accepted-step volume audit. The selected cell's every accepted entropy
increment is positive. Its cumulative absolute geometric contribution to
log(P/rho^gamma) is7.47e-9 and7.12e-9, below the predeclared1e-6 interval budget.
The conservative sum of all-cell maximum volume residuals also stays below
1e-6. Every retained mass remains exactly fixed, U/V stay positive, and the
existing contact-adjusted energy ledger remains at roundoff scale.

At2.839, endpoint versus stages differs by0.0292% in selected density and
0.0388% in selected pressure; mass-weighted global L1 differences are3.67e-8
and5.07e-8. Both differ much more from the uncorrected local trajectory.
That difference is not an error against an independent physical reference.

Evidence: Q/assessment/completed.json, Q/continuation_assessment/completed.json,
per-run process_result.json, audit.json, exact checkpoint lineage, and source/
binary hashes. The figure Q/continuation_assessment/thermodynamics.png (and PDF)
uses actual accepted states without smoothing; plot.py reproduces its layout.

## Scope and remaining physical tests

This removes the measured finite-step geometric source of entropy loss in
this selected interval. It does not restore the roughly65% entropy loss already
present in the incoming saved state, prove nonnegative entropy production for
every cell, eliminate physical/numerical viscous heating, or qualify wake
morphology and stationary-wall work. A closed contact-adjusted energy ledger
is still distinct from physically correct work at a stationary wall.

Endpoint-only can retain second order for smooth dynamics and a regular
volume-weight map: its predictor constraint error O(h^2) changes the final
update by O(h^3). Existing noncontact-annulus tests support that conditional
claim. Contact jumps, ill-conditioned geometry and limiter changes can violate
those assumptions. Enforcing the stage volume directly supplies stage2 with
the intended density/pressure, whereas an endpoint correction cannot undo
forces and work already evaluated at an incorrect predictor state.
Neither option by itself establishes a global SSP/entropy proof. Earlier
whole-domain wall errors showed about first-order temporal behavior and remain
unresolved. Close agreement in this interval does not settle that question.

The current dimensional goal remains the same physical volume in geometry,
density, pressure force and pressure work. One-dimensional moving faces and
multidimensional capacity inverses are different realizations; their common
conservation and limiting behavior still require cross-test qualification.
Finish Cylinder and the other2D cases before full3D fluid evolution.

## Cost and reuse experiment

The two-segment endpoint/stages runs took101.66/159.17s including two process
startups on the shared Ada. Their inverse calls consumed44.26/93.03s. These are
diagnostic costs, not matched-accuracy GIZMO timings. The <=2x GIZMO objective
has not been demonstrated.

On three saved actual nonsymmetric anchored Newton matrices, separating the
solver phases gives about0.34--0.38s for repeated structural analysis,
0.20s numerical factorization and0.03s triangular solve. Reusing the same CSR
analysis while refactorizing each new numerical matrix reduces the latter
calls from0.60--0.64s to0.25s. Raw-face true residuals remain below1e-15, there
are no perturbed pivots, and every stored next-W is identical to fresh analysis
for these matrices. All original full rows, anchors and self-image cancellation
are retained; no SPD approximation or fixed gauge row is substituted.

This small saved-matrix result motivates the private native reuse comparison.
It is not an end-to-end speedup or universal equivalence claim. Pattern changes
invalidate the factorization; each new coefficient matrix is numerically
refactorized, and error/residual failures release the cache.
Evidence: Q/factor_profile_v1/completed.json and phase timing logs.

## Actual native symbolic-reuse comparison

The same two native modes were tested to2.837925 with cached structural
analysis; numerical factorization is still performed for every new matrix.
Both complete17 steps with exactly retained masses and the original energy
ledger. Only exam_gpu_extract.o differs from the prior private archive.

| Mode | Structural reuse / solves | Inverse seconds, fresh -> reuse | Observed total seconds, fresh -> reuse |
|---|---:|---:|---:|
| Endpoint | 30/35 | 22.13 ->10.25 | 49.64 ->38.40 |
| Predictor + endpoint | 73/78 | 49.08 ->20.65 | 81.40 ->55.28 |

Each run rebuilds analysis5 times when the exact CSR pattern changes. Actual
raw-face linear residuals remain <=2.78e-13 with no perturbed pivots. The
precision-stop signed closure, actual-volume ceiling and accumulated-volume
budget all pass. The cached step-volume maximum is3.48e-9 for endpoint and
4.37e-9 for stages, still within the declared1e-8 ceiling.

These are numerically close trajectories, not bit-identical results. Maximum
endpoint density/pressure differences from fresh solves are3.79e-10 for
endpoint mode and7.38e-10 for stages. Stage-mode adaptive times differ by up
to3.92e-12, although final times agree exactly. An additional analysis assertion
requiring schedule differences below1e-12 failed and was replaced by explicit
reporting of those differences, not a claim of identical schedules. The
predeclared physical volume/precision bounds were unchanged.

Total times were measured with different concurrent GPU workloads; only the
solver decomposition directly identifies the saved repeated-analysis cost.
No matched-accuracy GIZMO speed claim follows. Evidence is
Q/symbolic_reuse_v1/{completed.json,stages_completed.json,assessment.json},
source/binary manifests, full logs and separate phase timings.

Both longer uncached inverse trajectories remain saved at2.839; both cached
qualification trajectories are saved at2.837925. All current runs are terminal
with exit0, and no background production or3D run was launched. The next
physics decision is to qualify intermediate-state pressure/work and stationary
wall accuracy, using prior event history and suitable saved states. Continuing
an inherited low-entropy trajectory alone cannot establish accuracy from IC.
Use the cheaper solver as a numerical implementation improvement; do not select
endpoint-only over stages solely because it costs less. Canonical-method
adoption, the remaining2D suite, full3D and the GIZMO performance gate remain open.

## Earlier positive-volume contact event: independent short comparison

The same cached native executable was tested from the already saved
t=0.9749999999999952 checkpoint to0.976499999999995, with126696 cells and
three accepted steps of cap0.0005 in each mode. The selected material ID is
64062, distinct from the later negative-volume event's ID63222. The actual
cylinder remains centered at(2.5,2), radius0.5. No IC restart or second-paper
campaign was needed; this is a regression check for the existing Cylinder work.

| Mode | Maximum all-cell step volume residual | Selected entropy-proxy increase |
|---|---:|---:|
| Current executable, inverse off | 1.68330e-2 | 0.689255% |
| Endpoint inverse | 1.20529e-13 | 0.00998745% |
| Predictor + endpoint inverse | 1.47426e-13 | 0.0107476% |

The off control differs from the older executable's final density/pressure by
at most3.51e-13 relatively; it is not byte identical because the current
boundary representation differs. The new comparison uses the same current
executable, IC/BC, saved start and final time for all three modes. Every retained
mass remains exactly fixed; U/V are positive and no generator penetrates.
The original contact-energy ledger closes below6.7e-15. All inverse stops in
this earlier interval satisfy the original strict tolerance: no precision-
limited stop is used. No U, mass, position or momentum repair is applied.

The selected geometric contribution to delta log(P/rho^gamma) falls from
0.00676713 to about1e-14; accumulated absolute values are below7.6e-14.
This shows that the same volume/pressure-work consistency removes a measured
positive geometric entropy defect as well as the previously diagnosed negative
one. It does not prove that all remaining entropy production is physical.
The thermal discrepancy remains about8.4--8.5e-5 in log units, and stationary-
wall physical energy defects remain approximately-2.45--2.47e-7, versus
-2.52e-7 without the inverse. Closed ledgers therefore still do not qualify
stationary-wall thermodynamics or morphology.

Endpoint versus stages differs by0.0658% in selected density and0.0914% in
selected pressure at the common endpoint. This short comparison does not
select a universally preferred mode or establish temporal order. Endpoint
and stages shared CPU cores/GPU concurrently, so their wall times are not
performance evidence. The inherited ID63222 entropy-loss history is separate.

Evidence: S.parent/contact800_t0975_v1/capacity_compare_20261108_v1 contains
the run.py and assess.py scripts, source/checkpoint/executable hashes, complete
logs, effective parameters, three exit0 records, assessment/completed.json,
and terminal_audit.json. No source modification, manuscript update, production
restart, full3D run or second manuscript was performed for this comparison.

## Inherited entropy-loss interval from saved history

The current common_swept_cfl lineage's retained ID63222 has a large
decline between available saved observations at t=2.5415 and2.7915. Its
A/A_IC falls from1.641003916 to0.592560931: A decreases by63.89% relative to
the earlier endpoint. Its mass is exactly unchanged and it is already wall-
adjacent at both times (radius0.5000000916 and0.5000077924 around(2.5,2)).
Earlier saved values are approximately1 through1.041, then1.625--1.641 from
1.291 through2.5415. A literal monotonicity assertion on the early plateau
failed; every adjacent saved-state increment is now reported explicitly.
Sparse saved states cannot establish the exact first
negative entropy increment or its cause, and monotone snapshots do not prove
every intervening step was correct.

The nearest existing source before this large loss is
S/common_swept_cfl_v1/continuation_t5_v1/cylmat.007029.dat at2.5415.
Its stream/exterior companions exist and their hashes are recorded. Do not
restart from IC or substitute the earlier cap4 branch's .931/1.0205 or its
different2.6715 checkpoint. The .975 ID64062 test above is a separate local
sign regression, not the source of ID63222's inherited loss.

Root independently extracted the IC and both bracket endpoints with the
existing compiled native reader. The full saved-history table uses the696-byte
ABI and matches those extracted fields exactly. A/A_IC consistently uses the
effective stored gamma=1.399999976158142 at all times; the first analysis that
mixed the earlier nominal-gamma IC normalization is preserved separately.
Evidence: A/historical_bracket_v1/{completed.json,history.json,assess.py,
record_history.py}, with source/extractor/checkpoint hashes. No new evolution
over this interval has been launched. Its full entropy loss has not yet been
attributed to the later2.8379 event mechanism.
