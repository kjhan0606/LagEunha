# Positive-rule native connection and dependency — 2026-10-22

Goal ACTIVE and unachieved: suppress artificial Kepler shear compression/heating,
reach t6/8/10, ultimately the correctly normalized longest GIZMO comparison.
The preceding question-only reply was a status explanation, not new physical
validation. This turn changes the native connection/build and adds a live,
single downstream dependency; it does not claim a longer accepted trajectory.

## Project and preserved state

Supplied cwd is `/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce` and is not the
original repository. Original source remains
`/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src`.
Read-only git through syntax revalidated HEAD8e25dfe0c2e41a7374c11e2e3102da2a8b028625, remote
`git@github.com:kjhan0606/LagEunha.git`. No original-source mutation, git repair,
commit, push or shared-context update occurred. Prototype code is in this
owned rk4_diagnostics tree, with separate build snapshots.

Accepted actual trajectory is unchanged at t=8.2355880701376864e-5, four RK4
steps. Reuse checkpoint
`material_full_cached_continuation_v1_ada_mpi4/run/continuous_step4/accepted_r`
and actual next high-rule k1
`material_full_resident_rhs_v1_ada_mpi4/run/next_RK_stage1/prepared_k1_r[0-3].bin`.
The saved k1 force and mass solve must not be recomputed for setup.

## Running preprocessing — verified live, not restarted

Original MPI launcher2166274 and GPU ranks2166305–2166308 were revalidated
live on LagEunha. Producer remains
`material_full_positive_rule_v1_ada_mpi4`. At report preparation, committed
prefixes contain 176184/1138432triangles (15.48%), all compressed:
315739648original points ->38486736selected points (8.20x in this prefix).
Projected total wall is5.19–5.40h
from preprocessing start. Prefix performance is not a completed-rule result.
Fields/time/M/generating entropy remain unchanged during preprocessing.

## Native connection built, not yet physically exercised

Final current build: `material_full_positive_transition_v2_ada_mpi4`.
CUDA and C++ compiled together; native original C entry and linkage retained.
Build status is BUILD_PASS_NO_PHYSICAL_EXECUTION. The prior v1 and v2 first
build are compilation-only history; no accepted evolution was replayed.
The inherited ifx option warning is recorded in build_link.log; linker exits0.
The enlarged MFEvolvedBlock ABI includes max_J, so older CUDA.o must not be
mixed with these host headers.

Implemented in the current compiled snapshot:

- Native whole-rule loader refuses live/partial output. It checks complete
  producer/native status, every ordered source record and checksum, all signed
  moment residuals/rounding arrays, positive original-point subsets, every
  retained native contraction gate, and global original/selected point counts.
  Fallback triangles retain the original full positive rule.
- Raw/cache CUDA evolution select the same positive reference subset. Initial
  mass and original fixed binary80 M stay unchanged. CPU shadow uses the same
  selected rule. A rule-mode switch invalidates the cache.
- Four actual RK4 stages compare low and high operators at identical x/v stage
  coefficients. Only low rates advance the state. First high k1 comes from the
  retained real calculation; missing all-source summaries are metadata-only
  original-node evaluations without force/CG. Stages2–4 compute the high RHS at
  the actual low-stage states, without a second high trajectory.
- Every stage and endpoint retain all-source U/Phi, pressure/density,
  volume/rate, delta/delta_dot and mass-metric acceleration-difference CSVs.
  The acceleration norm uses original M, hence also measures mass-implied load
  discrepancy in its dual norm. Unsaved raw assembled high k1 loads are not
  claimed to have been recovered exactly from rounded accelerations.
- Point pressure/density discrepancy estimates use analytic bounds on every
  original reference node, in addition to source averages. These are numerical
  comparison estimates, not a proof of continuous-map positivity/injectivity.
  The T10 additional quadrature budget remains1e-3; stage rate and cumulative
  limits scale with physical dt. It is not a stability/error-growth theorem.
- Original-node CFL information is retained during this actual comparison,
  conservatively rescaled for the changed corrected-j sound speed. A subset
  maximum alone cannot authorize a larger dt.
- U/Phi/mass integral changes at the quadrature transition are separate signed
  numerical ledger entries. No physical heat/coefficient/entropy correction is
  added. A rejected step restores coefficients, time, physical/numerical
  ledgers AND the previously committed rule mode.
- v3 restarts bind the exact immutable rule generation and numerical ledger.
  Accepted v1/v2 imports remain supported. The actual accepted v3 checkpoint
  will be reloaded and compared without replaying its evolution interval.
- Selected-step cost and one-time high-rule reference work are timed separately.
  Cost extrapolation explicitly excludes future guard cost and is optimistic.
  Both actual paired wall and the selected-cost estimate are retained.

## Live single dependency

Dependency monitor PID2296727 was installed for original producer PID2166274,
start_ticks1083750537. Files: dependency_pid.txt, dependency_state.json,
dependency.log under the v2 directory. Before treating it as running, inspect
its actual /proc or ps handle; a JSON state alone is not a live-process proof.
It polls the same original producer at60s cadence, never duplicates/restarts it.

On actual ready completion it invokes the already built native runner once.
The runner independently checks producer/source/binary/immutable/accepted/k1
hashes, hashes the completed four rule shards, rejects a busy GPU, and parses
all records natively. The actual purpose is ONE paired transition step, using
Ada1/MPI4/OpenMP8. Wrapper wall cap7200s; native cost-projection budget86400s;
planned new outputs <=4GiB with >=8GiB free required. No point atlas is written.
Failure, partial candidate, missing original handle or failed transition ends
this dependency without an automatic retry.

A native exit0 alone is not sufficient: contract.json requires step5 accepted,
all4stage comparisons, endpoint comparison, finite physical gates and v3
restart roundtrip. Then inspect actual cost/physical evidence. The code still
stops after this comparison for the unimplemented long-production map guard,
or earlier on measured cost. This is not a long trajectory release.

## Next decision and remaining scientific work

Continue the adopted ASTRA_LONG_EVOLUTION_COST_DECISION_20261022.md bundle;
no new routine Astra call was made. If accuracy/cost fails, retain valid work
and stop optimizing this10.8M-DOF representation as specified there; a strategic
space/mesh-lifetime redesign is next, not tolerance/order/dt sweeps.
If it passes, remaining full-node/continuous-map safeguards and sustainable
cost must be resolved before the t1/24h trajectory and t6/8/10 validation.
Spectral enriched-space CFL, material-topology shear lifetime, conservative
repartition/recipient thermal compatibility, W0 dynamic Voronoi reduction,
shocks and paper boundary/normalization comparisons remain unproved. Re=inf,
no steering, no heat reset, no independent spin-energy reservoir remain.

Created 2026-10-21T23:35:08.182501+00:00.

## Diagnostic continuation

See RESULT_SHEAR_AND_COMPRESSION_DIAGNOSTICS_20261022.md for the quantified
exact material shear, central-gradient dt comparison, and original-rule
step4bulk/subcell compression decomposition. No physical trajectory or live
producer/runner changed.

## Latest dispatch update

Ada-only dependency2296727 is superseded; current dependency PID2351718
is material_full_positive_dispatch_v1. See RESULT_H200_PAIRED_DISPATCH_20261022.md
for the H200-first single execution, terminal-pending cancellation/atomic-claim
Ada fallback, and actual verification. Original producer2166274 is unchanged.
