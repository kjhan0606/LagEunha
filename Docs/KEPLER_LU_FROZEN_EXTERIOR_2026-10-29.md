# Recovering the historical n64 Kepler moment preflight

The old `origin/codex/physics-validation-20261017` research branch recorded
a first-RHS failure in an auxiliary Kepler moment model. Its n16 failure was
an **active–active interior face** with a reconstructed normal-velocity
gap 0.0355921 exceeding the positive-state vacuum threshold 0.0354348. The
analytic reference velocity was continuous there. Separately, the archived
n32 run stopped on an initial-pressure accuracy assertion; it did not reach
an independently accepted fluid step. These are distinct failures.

The preserved n64 run also failed on a vacuum exception. Re-reading its
hash-verified saved initial state using the original diagnostic identifies a
different first failure: a **frozen outer cell and external wall ghost** at
r=3.23276. The reflected states have normal velocity gap 0.414444 versus
the vacuum threshold 0.408248; the same failure occurs before and after the
first gravity half-kick. See [original n64 diagnostic](KEPLER_LU_FROZEN_EXTERIOR_FAILURE_2026-10-29.json).
This identifies the first failure, not all potentially failing faces.

The old model treats frozen cells as prescribed reservoir states, advances
only active cells, and books active–frozen momentum/energy/angular exchange.
An isolated [candidate patch](KEPLER_LU_FROZEN_EXTERIOR_CANDIDATE_2026-10-29.patch)
omits only an external wall face whose owner is frozen. It preserves every
active–frozen face and every active exterior face. No saved IC, original
research source, or running native solver was changed. In six small tests
(no/mixed/all frozen, scalar/batched flux), the active RHS and reservoir
ledger agree with the old implementation to max 2.09e-16; four existing
frozen-boundary unit tests pass. The [comparison script](test_kepler_lu_frozen_exterior_candidate_20261029.py)
uses the archived source and isolated candidate at their recorded paths.

With this candidate, the original n64 initial state passes the first RHS both
before and after the gravity half-kick; see [RHS diagnostic](KEPLER_LU_FROZEN_EXTERIOR_RHS_2026-10-29.json).
One actual split fluid/gravity step from that saved state then completed for
4,084 cells, with 1,878 frozen and `dt=1e-5`. Maximum stage net-mass flux was
6.52e-14, frozen state remained exact, active generator steering defect was
zero, and all cell thermal energies computed from the stored moments were
positive. The separately computed gas-energy minus boundary/gravity budget
residual was -7.98e-18; the central angular-momentum minus boundary budget
residual was -5.14e-17. Source hashes, checks and the preserved checkpoint
are in the [step audit](KEPLER_LU_FROZEN_EXTERIOR_STEP_2026-10-29.json)
and [auditor](audit_kepler_lu_frozen_exterior_20261029.py); the isolated
[runner](run_kepler_lu_frozen_exterior_candidate_20261029.py) resumes the
existing initial state instead of regenerating it.

This is a **one-step historical Python research prototype**, not the
user's CUDA/MPI production solver and not a long Kepler result. The old
frozen-outer-cell boundary is not the present no-steering fixed-wall setup;
the candidate does not change that physical contract. The n16 interior
vacuum, n32 accuracy failure, endpoint moment support/realizability,
subsequent steps, topology changes, timing efficiency, genuine vacuum
Riemann support and physical solution accuracy remain unproven. One step
took 69 s on a shared CPU host, so this Python implementation cannot be
promoted as an AREPO-speed simulation. No floor or state repair was added.
