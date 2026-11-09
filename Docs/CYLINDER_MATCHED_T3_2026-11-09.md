# Matched Cylinder endpoint at t=3.00025

This is first astronomy-draft validation. The second fluid-mechanics paper
starts only after the first draft exists; its evidence is retained here without
starting another manuscript or an unrelated campaign.

## Executed method and lineage

The existing private endpoint-volume inverse continued its saved t=2.7915
checkpoint to3.00025 with the same C/CUDA binary, fixed retained cell masses,
SSPRK2, original radial reflection/contact-energy accounting and positive-gap
representation. No LU-preconditioner, new force, thermal repair or limiter was
added. The inverse is active from2.5415 in this lineage, not from the initial
condition. Older baseline and failed trajectories remain separate evidence.

428 accepted steps and9 retries completed with exit0. There are126,696 owned
cells;125,416 material IDs are shared between the two endpoints, with unchanged
masses. U,V remain positive. All126,696 actual polygons are convex, with no
fluid intrusion into the analytic circle and no boundary joining gap.
The maximum local finite-step volume discrepancy is9.8437081e-12;
the maximum step energy-ledger residual is1.3374718e-14.

The selected earlier-loss cell63222 has entropy-proxy ratio1.0017430 over
this interval. This is a local history result, not global entropy validation.
2,028 retained IDs have delta log(P/rho^gamma)<-1e-6; the most negative change
is-0.0097894. Negative cell-average entropy changes alone do not establish a
physical entropy violation or its cause. No entropy repair is authorized.

## Matched independent fields and shape

Native and both SU2 references have exactly the same requested time3.00025.
Both reference ancestries complete to5 without reported nonphysical recovery.
They use h=.025 and.05; two-mesh sensitivity is not an exact-error bound.

The bow-shock centreline locations are x=1.744 native and1.736 fine reference;
the difference.008 is below the initial native spacing.025. Across seven
horizontal cuts, measured differences are.002--.012. Pressure drag is
1.36858059 versus1.36977869, a0.08747% difference. Reference-mesh drag
sensitivity is about0.596%. These do not certify the wake or all wall pressures.

Actual-native-cell-average field comparisons:

| Region | Density relative L1 | Pressure relative L1 | Reference density/pressure sensitivity |
|---|---:|---:|---:|
| Prespecified upstream1834 material IDs |2.283%|0.739%|1.800% /1.721%|
| Wake7544 complete cells |3.049%|4.467%|2.091% /3.019%|

The same native polygons are used for both reference integrations. Pressure
is the volume average of physical P, not a residual after subtracting averaged
bulk kinetic energy. Projection changes and reference sensitivity must not be
subtracted from the listed errors. The unsmoothed native density figure still
shows narrow wake structures; the stored pressure cuts also differ locally.
These differences are not assigned to KH, noise or boundary error by image
appearance alone. Global drag agreement is not a reason to discard them.

For the sampled interior region, the common-fluid density/pressure/velocity
differences are1.478%/2.162%/1.466%. It excludes.0043 of physical fluid area
behind the approximate native obstacle. Thus these are common-fluid norms,
not full-domain norms. No lookup failure or overlap is silently omitted.

## Remaining stationary-wall limitation

The actual excluded obstacle area is.790249997 versus pi/4=.785398163:
the excess is.004851853 (about0.618% of the circular area). Its maximum
radius is.5255534 rather than.5. Maximum geometric normal speed is.814124,
and instantaneous pressure-wall power is-1.4124124e-4. Convexity and no
penetration therefore do not prove an exact stationary circular boundary.

Over this interval, image-wall work plus contact kinetic-energy change is
-3.4178771e-5, although the ledger closes to roundoff. This physical wall
defect remains distinct from the corrected volume/pressure-work consistency.

## Separate offline solver-cost evidence

Previous LU factors were used only as a right preconditioner, with current
raw-face J matvecs and the existing1e-10 residual criterion. On six retained
matrices, measured total solve work decreases1.9324s to1.0577s. Five reused
solves require1--2 iterations. Maximum residual is5.44e-11. One stored-W
comparison differs3.47e-18; the other five are bit-exact. These are offline
linear algebra measurements, not an end-to-end physical speed qualification.

Ten additional synthetic cases pass: bounded-eight-iteration nonconvergence
falls back to current numerical factorization; changed dimension/pattern,
current-versus-old operators, zero RHS, invalid input rejection and recovery
are exercised. A deliberately injected old phase33 error successfully clears
the cache and uses a fresh factorization. This injection is not evidence of
an observed native solver failure. Native row-ID/anchor eligibility and
time-evolution qualification are still missing; this candidate is not adopted.

In the prototype logs, same_pattern=0 on a successful preconditioned solve
means the pattern comparison was bypassed, not that the pattern changed.
Its residual-check cost is included in solve, not certificate. Retain these
definitions when comparing timings. The earlier two-hop structural envelope
remains rejected on cost; it is a different algorithm.

## Evidence locations

Under Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/
cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/
production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/
entropy_trace_v1/prevention_t265_v1/continuation_t300025_v1:

- endpoint/audit.json and process_result.json: execution and restart state.
- assessment/endpoint.json and endpoint_geometry.json: heat/volume and polygons.
- morphology_v1/comparison.json and density_comparison.png: same-time fields.
- upstream_reference_v1/projection.json and wake_reference_v1/projection.json.
- review_record_v1/summary.json: compact record and source hashes.

Linear tests are under entropy_trace_v1/volume_coordinate_v1/precision_stop_v1/
lu_preconditioner_20261109_v1, including fallback_cases_v1/completed.json.

The same endpoint-method continuation to3.50025 is a bounded next Cylinder
validation interval. Full2D, MPI/3D and matched-accuracy GIZMO runtime remain
open. The first astronomy draft remains the priority.
