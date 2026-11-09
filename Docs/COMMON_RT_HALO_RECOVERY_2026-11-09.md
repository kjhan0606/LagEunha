# RT halo capture and qualified recovery — 2026-11-09

The frozen common RT128x256 run stopped at accepted t=3.095874084472433, with its last ordinary checkpoint at3.0. Its original failure remains recorded. This investigation changes no hydro law, physical boundary condition, pressure/work operator or timestep threshold.

## Actual failed geometry captured

`common_rt128_halo_capture_20261109_v2` changes only exam_gpu.o in the frozen native archive. The capture reads synchronized device polygons BEFORE the failed geometry returns, all double-precision x/y/W and IDs, linked-cell offsets/indices, periodic chart owners and RT mirror sidecar. It never calls the old float-W/clamped diagnostic or uses stale FaceCSR. First valid base geometry and first rejected predictor are preserved. Original failure behavior reproduces exit99.

The selected physical cell is ID15104 near the periodic x boundary, not a physical y wall. At predictor t=3.095936584472433, the conservative radius is.020667624048664478 and exceeds the left computational grid by3.75903891354501e-7. Its six real support faces are present; three use periodic images. There is no artificial support face.

Independent clipping and all294912 nearest periodic/reflected image checks find no omitted plane cutting the native polygon within the declared floating-point bound. Farther images have a positive distance/weight certificate. Relative area difference is7.82e-17 and vertex discrepancy1.25e-18. Absent-image minimum slack is9.15e-5. A separate check of the actual uploaded arrays gives maximum plane residual9.45e-21 and minimum slack4.31e-4 outside the actual linked-cell search window. All35857 uploaded slots occur once in CSR; owner weights, aliases, reflection and coordinate roundoff are checked. These are strong numerical certificates, not exact symbolic predicates or proof for future states.

Conclusion for this failure: an otherwise valid cell was rejected because the conservative search sphere exceeded the fixed allocated search grid. Halving dt approached the same capacity boundary rather than correcting a physical instability.

## Capacity-only controls and restart qualification

Existing GFS_KEPLER_HALO_CELLS controls both the computational grid and periodic/reflected padding. Counts3 and4 were tested on the same restored step-start input and original dt. Both complete the formerly rejected step with zero retries; all13 saved state fields are bitwise identical. Maximum audited volume defect2.9374e-12 passes the unchanged1e-11 contract. This does not prove that width3 suffices for all future times.

Before using the recovered input, `recovery_validation/report.json` verifies source rollback semantics, bitwise x/y/W agreement with native phase1 capture, exact original per-cell masses, positive U/V, native-volume agreement3.2833e-16 and gravity-inclusive total-energy drift-8.3463e-15 relative to the original IC. The diagnostic input is no longer treated as valid merely because a file exists. A subsequent one-step check with the ORIGINAL unmodified portable executable also matches the diagnostic build's complete saved state bitwise.

The original t3 checkpoint and failed trajectory remain unchanged. A detached finite continuation now starts from the qualified original-binary checkpoint t=3.0959365844724327, using halo3 and targets3.3,3.6,3.9,4.0. Same binarySHA11a0234acc882f6f2fc405b0053c2e5210262863be665f88cc5ae91b821ddf4e. Supervisor958973/native959091 were actually observed live after launch; inspect current /proc and status before future claims. Output bound15 checkpoints including copied inputs,342.10MB, with149TB free at launch. A new failure stops the batch; no automatic retuning/restart is enabled. This is the existing common control, NOT the rejected log/tangent candidate.

Live and terminal records are under `recovery_continuation_v1/{status.json,results.json,segment_*/}`. At the last observed accepted line it had reached3.1228115844723701 without retries. That observation is historical, not a perpetual live-state assertion.

## General correction still to implement

The existing generic power-search-box helper and pure-halo shortfall API can support measured capacity growth instead of a hand-selected constant. Extend grid and periodic/reflected padding together, re-evaluate the same stage before publishing the fluid update, and preserve failure on mixed invalid geometry, overflow or no growth. Reuse existing retry/rollback machinery. The present evidence qualifies recovery at halo3; it is not an implementation of that generic automatic capacity policy or MPI/3D verification.

## Reference limitation

The separate GIZMO diagnostic continuation reachedt4 with exit0 and its processes terminated. The strict mirror failure remains: maximum mirror position discrepancy.0097773, velocity.0753795, density relative.0168172 and internal energy relative.0117099. It is not a qualified reflecting-wall reference, even though its core conserved quantities pass. Preserve it as a diagnostic, and do not use its duration as the matched-accuracy GIZMO performance certificate.


## Automated capacity correction completed

The subsequent same-stage common search-capacity implementation, exact RT replay and periodic/cylinder/Noh regressions are documented in `COMMON_SEARCH_CAPACITY_2026-11-09.md`. This supersedes the implementation TODO above for the qualified single-rank periodic/RT driver only; the original ongoing RT remains on its frozen executable.
