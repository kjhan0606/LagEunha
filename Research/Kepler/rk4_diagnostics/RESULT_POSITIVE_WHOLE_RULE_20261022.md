# Full positive cubature preprocessing and whole-record audit

2026-10-22, verified about 13:04 KST. This is a completed preprocessing and
stored-record audit bundle, not a new physical evolution interval.

## Actual full producer result

`material_full_positive_rule_v1_ada_mpi4` finished with native and runner exit
codes 0. Its original MPI launcher 2166274 and ranks 2166305--2166308 are now
absent from actual `ps`. No restart or duplicate producer was launched.

- All 1,138,432 original triangles processed and compressed.
- Original positive points: 2,040,256,512.
- Selected positive points, including any original-rule fallbacks: 249,691,129.
- Point reduction: 8.171121337674755x. This is NOT a measured step speedup.
- Actual producer seconds: 18,966.2571 (about 5 h 16 min).
- Anchor integral discrepancy: delta U/U = 2.2939011722919333e-18;
  delta Phi/U = -3.882977377104479e-14.
- `candidate_ready_for_actual_RK_oracle=true`; `production_rule_committed=false`.

Authority: producer `contract.json`, `native_exit_code.txt`,
`runner_exit_code.txt`, and completed `run/positive_rule_r[0-3].bin`.
These numerical integration weights are separate from geometric Laguerre
weights. No mass, entropy, state coefficient, accepted time, or thermal ledger
was reset to accommodate the rule.

## Independent whole-record audit

Executed the existing `run_material_positive_whole_audit_v1.sh` ONCE on
authorized lageunha. Output: `material_positive_rule_whole_audit_v1_mpi4`.
Foreground tool session 69972 completed with exit 0; native audit exit 0.
Final status: `WHOLE_RECORD_AUDIT_PASS_PAIRED_PHYSICS_STILL_REQUIRED`.

The native MPI4 inspector checked all 1,138,432 triangles and all 2,040,256,512
original points accounted for, yielding the same 249,691,129 selected points.
Header/record/support/all-row and arithmetic-bound checks passed. Its
actual-anchor contractions reconstructed from saved signed moment residuals
and original mass data agree with saved diagnostics: maximum difference 0,
maximum contraction roundoff bound 2.2611746473859075e-26.
Actual native audit seconds: 19.1464664; wrapper hashing time is additional.
All source/binary/immutable input checks and final completed-rule SHA256 checks
passed; no completed rule changed during the audit.

This does not query new physical fields or replay a force/mass solve.
The force proxy threshold is checked from its saved report only; actual
same-stage full-rule force and acceleration comparisons remain mandatory.
No rule or new physical step has been adopted by this audit.

## Existing downstream dependency

The existing dispatcher PID 2351718 was confirmed live by actual `ps`.
It submitted H200 Slurm job **433999**, requesting one H200 card and explicit
128 GiB memory under the adopted fixed policy. Actual `squeue -j 433999`
at about 13:04 KST: `PENDING (QOSGrpGRES)`. The snapshot alone must not be
used as evidence of future liveness.

Do not submit another physical job. Existing pending-only 600 s cancellation,
authoritative terminal confirmation, and exclusive physical claim control the
Ada fallback. A start race must preserve the running H200 job. Never remove
the claim or restart a physical transition on an observation timeout.

Next actual operation is paired RK4 step 5 from audited step 4, using saved
high-rule k1 and original high evaluations only at the actual low stage
states. Inspect all four stages, endpoint, finite physical gates, restart
roundtrip, quadrature budgets, and measured selected cost. Native exit 0 or a
wrapper PASS alone does not authorize long production.

Apply `ASTRA_LONG_EVOLUTION_COST_DECISION_20261022.md`: measured sustained
step cost must improve at least 10x, or the optimistic measured t1 estimate
must be at most 24 h. If neither holds, end further optimization of the
unchanged 10.8M-DOF representation and address conservative representation
and mesh lifetime. Do not substitute tolerance/order/dt sweeps or blind long
continuation.

## Physical scope remains incomplete

Execution update, about 13:15 KST: job 433999 was cancelled while pending at
13:12:45; actual Slurm accounting confirms CANCELLED, runtime 0, no start.
Existing dispatcher acquired the exclusive Ada claim at 13:13:50 and launched
runner PID 3026814. Actual ps confirms that runner and native MPI launcher
3028372 (timeout 3028371, Hydra 3028377/proxy 3028395) live. Ada status is
`RUNNING_ONE_ACTUAL_PAIRED_RK4_TRANSITION`. This is startup/in-progress evidence,
not an accepted step or an oracle/cost PASS. Do not duplicate or remove claim.

Accepted state remains four steps at t=8.2355880701376864e-5. No t6/8/10 or
long GIZMO-paper stability result exists. Actual whole-step speed, future
quadrature accuracy, enriched global map guard, conservative thermal
transfer/W0 reduction, shocks and boundary/paper normalization remain open.
Artificial shear heating is not established as solved.

Original user repository identity was reconfirmed at
`Kepler_stencil5x5_Nx256_A100_call33720_20261020/src`,
HEAD `8e25dfe0c2e41a7374c11e2e3102da2a8b028625`.
Original code/results and shared global context were not modified.
