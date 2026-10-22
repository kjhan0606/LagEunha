# Independent positive-rule record inspection — 2026-10-22

Goal ACTIVE. Actual accepted trajectory remains four steps,
t=8.2355880701376864e-5. No t6/8/10 heating solution is claimed.

The source/actual-state revalidation found newer completed cached and resident
mass work than the earlier conversational summary. RESULT_FULL_EVOLUTION,
RESULT_RESIDENT_MASS_AND_POSITIVE_RULE, RESULT_POSITIVE_RULE_CONNECTION,
RESULT_H200_PAIRED_DISPATCH and ASTRA_MATERIAL_SHEAR_LIFETIME_AND_GUARD are adopted.
The obsolete attempted cache patch was rejected atomically; it changed none of
the existing operator/runner files. No existing calculation was restarted,
no prepared high k1 was replayed, and no original repository file was changed.

## Actual new independent evidence

`audit_positive_rule_records.cpp` is a separate native saved-record inspector.
Its frozen source/build is `material_positive_record_prefix_audit_v1_mpi4`.
The build and MPI4 execution on authorized lageunha ended0; SSH tool4826 is
CLOSED0. Four ranks used OMP1 and inspected one already fsynced prefix, exactly
2048records/rank. Runtime0.647582997sec. This read no uncommitted suffix and did
not evaluate field quadrature, solve mass, propose weights or advance a state.

- All8192records preserve exact actual accepted-state and prepared-k1 bindings,
  original source IDs/shape count and every original constraint row.
- Every raw checksum includes actual binary80 padding bytes. Signed row
  residuals are finite; arithmetic bounds are finite/nonnegative.
- Every retained weight is finite/strictly positive and references a unique,
  ordered original positive point. Full-rule fallback handling is explicit.
- Original14680064points correspond to1750827selected points in this prefix.
- Independently reconstructed actual-state dK/U, dL/(2K), dWork/U, dU/U and
  dPhi/U agree with native stored contractions exactly in binary80 (maxdiff0;
  max estimated contraction rounding allowance1.20359223582696324896e-26).
- The force-proxy threshold is only checked against its saved value; its actual
  force/acceleration validation still requires the already planned paired RK4.

This is a prefix audit, NOT a whole-rule or physical evolution pass, and it
adopts no candidate. Signed dU and dPhi sums are reported without resetting
physical heat, M or entropy. The whole mode (`prefix_limit=0`) refuses live or
partial output, requires actual terminal complete producer status, inspects
all1138432records including fallbacks and EOF, and checks global point totals.
It can run after producer completion without repeating any physical queries.

## Live pipeline remains intact

At the saved metadata snapshot: 431928/1138432triangles (37.94%),
774032896original points ->94722289selected points (8.172x
within the committed prefix). Projected total preprocessing wall remains
approximately5.2–5.4h FROM START, not remaining time. This projection and the
prefix compression factor do not establish completion or sustained accuracy.
Native producer launcher2166274/ranks2166305–2166308 were separately verified
live by ps earlier this turn; dependency2351718 was also verified LIVE by the same explicit ps command
(elapsed50:20). The launcher elapsed02:00:51 and all four native rank command
arguments match the original producer. See live_process_verification.json.
Producer and dispatcher are untouched, with no duplicate monitor or job added.

Continue existing H200-first single claim policy after FULL ready preprocessing;
Ada fallback only after actual terminal pending-job cancellation with no claim.
One physical paired next step checks all four stage states, using saved highk1,
plus endpoint and restart. Independent whole-record audit remains to run at
completion. Neither this audit nor prefix8x compression authorizes long runs.
If paired accuracy/cost fails, preserve the state and follow the adopted
strategic cost decision instead of further order/dt/tolerance sweeps. Enriched
map guard/CFL/topology lifetime, recipient thermodynamics, W0, shocks, boundary
matching and normalized longest GIZMO comparison remain incomplete.

## Whole audit invocation prepared

`run_material_positive_whole_audit_v1.sh` now makes the complete independent
read-only audit one foreground command on authorized lageunha after actual
producer/runner success. Shell syntax and embedded Python AST pass. The
wrapper refuses live/partial output before creating an output directory,
checks source/binary/immutable/accepted/k1 and complete-rule hashes, invokes
the already built native inspector with limit0, checks all1138432records and
original2040256512point accounting, and records the final complete-rule hash.
It adopts no rule or physical state. It is not a watcher, GPUjob, new proposal
or rerun of the8192prefix inspection. The wrapper was NOT executed while the
producer is live. See whole_runner_readiness.json and WHOLE_AUDIT_COMMAND.txt.
Original producer/frozen paired runners/dispatcher are unchanged.
