# Common-method necessity and current Cylinder evidence

The method must be physically justified, dimensionally consistent and simple.
Do not retain a new term merely because it improves one test. The first
astronomy draft and its validation precede the second fluid-mechanics paper.

## What is common, and what is not yet connected

The volume inverse solves J deltaW = V_target - V(x,W), where J=dV/dW and
V_target uses the existing SSPRK2 volume rates. It is a geometric time-update
compatibility condition, not a new EOS, viscosity, internal-energy repair or
force law. The current native adapter explicitly requires single-rank 2D
Cylinder and folds actual aliases/images/exterior anchors into the full J.
It does not prove exact material-face trajectories or full physical accuracy.

The algebraic face coefficient A/(2d) is dimension independent. The actual
adapter is not: it uses x/y coordinates and cylinder boundary ownership.
General2D periodic use needs the constant-W gauge and compatible target sum.
Existing `Exam/gfs_material_direct_basis2d.h` and the frozen direct-weight
helper in `material_finite_slack_20261107_v1/certificate_direct_weight/`
already fix a gauge and verify the original full directed operator. The
static Gresho `actual_equal_volume/result.json` reaches 8.51e-13 volume error
in three Newton corrections; it does not validate SSPRK stage targets.
Do not reuse the old volume Newton's target rescaling, floors or W clipping.
No current volume-inverse KH/RT/Sedov/Noh/3D qualification was found.

## Continued Cylinder result

Define A and F as in `Docs/CYLINDER_LOSS_PREVENTION_2026-11-08.md` and
C=F/continuation_t27915_v1. All three runs completed t=2.7915000000000014.
The authoritative combined data and full-range unsmoothed figure are in
`C/combined_v1/`. Endpoint-mode cumulative maximum volume discrepancy since
t2.5415 is 3.202e-11 and its cumulative sum of per-step maxima is 4.893e-10.
ID63222's A_final/A_start is 1.007793. Stationary-wall work remains unresolved.

The current off control first loses entropy proxy, then gains it strongly;
its final A/A_IC is 2.51142, whereas the old historical trajectory at the same
time had .59256. These different current/historical trajectories must not be
conflated. The current off control has finite-step volume residuals as large
as4.09478. Scalar energy-ledger closure and final entropy sign alone would
miss this failure. Both current enabled modes keep the volume budgets and
avoid that large excursion. This is still not full Cylinder qualification.

The last interval's lowest retained-ID delta logA is about-.01308 at ID73672
in both enabled modes. Do not start an entropy repair from this scalar sign.
At t2.65, the earlier minimum ID60066 is a downstream cell outside the
prespecified upstream cohort. A new integration of the already saved fine
and peer fields over its actual polygon finds native rho~.12560 versus
.06516/.09605 and P~.034677 versus .012701/.021629. All three native modes
are almost identical there; volume correction has not fixed that wake error.
The independent reference is itself resolution-sensitive there. This is a
spatial comparison, not a matched material entropy history. Raw evidence:
`F/downstream_entropy_reference_v1/`, no new reference evolution was run.

Endpoint mode continues, unchanged, from its verified t2.7915 state to the
next existing independent-reference time t=3.00025. See
`F/continuation_t300025_v1/{plan.json,endpoint/process.json}`. Do not treat
a live job or a completed time interval as a global accuracy verdict.

## Computational change tested and rejected

The existing same-CSR cache was already on. In the completed t2.5415->2.65
comparison, endpoint/stages inverse symbolic analysis costs105.3/190.2s;
numerical factorization costs97.5/215.4s. Extra inverse work, rather than just
the mesh builder, is a major cost. Shared Ada totals are not a GIZMO benchmark.

A private offline prototype retains zero-valued structural slots out to two
graph hops, while copying every actual matrix coefficient exactly. It always
refactorizes numerically and checks the unchanged original raw-face residual.
It reuses a structural envelope only when all actual entries fit; otherwise
it rebuilds. It adds no physical parameter. This differs from the previously
rejected active-set union of physical weight constraints.

Six saved actual 126,696-row Newton systems were replayed on one CPU core:

| Solver strategy | Total seconds | Numeric factorization seconds | Peak RSS KiB |
|---|---:|---:|---:|
| Fresh symbolic analysis each call | 3.6391 | 1.2094 | 224984 |
| Existing exact-pattern reuse | 1.8441 | 1.1740 | 235420 |
| Two-hop structural envelope | 5.0875 | 3.9403 | 541400 |

The six systems happen to share a pattern; both caches reuse five times.
Thus this is not a representative long topology-change benchmark. However,
the envelope is slower even than fresh analysis on these systems, and its
numeric solve cost alone is too large. No native adoption follows. All
original-face residuals are below1e-15, and all stored next-W values match
the fresh solver exactly. Correct algebra does not justify retaining a slower,
more complex implementation. Preserve it as a negative performance result.

Evidence: A/volume_coordinate_v1/precision_stop_v1/
`symbolic_envelope_20261109_v1/` (sources, input hashes, commands, logs,
coefficients' results and completed.json). Canonical physics/source paths
and the running native executable were not changed by this prototype.

## Qualified private numerical reuse

The existing full-J volume equation now has a private previous-LU
preconditioner qualified on near-contact, ordinary and ownership-event native
intervals. It adds no physical term and retains the current raw residual and
nonlinear volume checks. Reuse requires the same positive CSR graph and
ordered ID/anchor basis. Ordinary and ownership endpoints are bit-exact;
near-contact rho/P differences stay below4.6e-9. This is a computational cost
reduction, not general2D/3D or GIZMO goal completion. The running3.50025
continuation retains its original executable. See
Docs/COMMON_SOLVER_REUSE_2026-11-09.md for tests, rejected prototype, actual
cost and the earlier KH GIZMO comparison's accuracy limitations.


## Later common 2D update (2026-11-09)

The Cylinder-only statements above describe the earlier evidence. Canonical source now includes an opt-in closed-domain inverse sharing the pressure adapter's resolved-face construction and certifying original equations in the unchanged physical-volume tolerance. Periodic KH/Sedov and reflected-wall RT pass the short native coupling tests; predictor-stage and temporal-refinement results, their limits, and the newly launched full KH128 trajectory are recorded in `Docs/COMMON_CLOSED_VOLUME_2026-11-09.md`. These results do not qualify all2D morphology, MPI,3D or matched-accuracy GIZMO speed. Endpoint-only is retained for the next full-history comparison because adding predictor corrections did not resolve the remaining field behavior and increased local cost. Signed excess-work diagnostics must not be mislabelled purely irreversible heating.
