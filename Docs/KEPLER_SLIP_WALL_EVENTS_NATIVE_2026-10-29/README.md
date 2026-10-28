# Native RK4 polygon-wall events: frozen v1, 2026-10-29

Follow-up: [v2](../KEPLER_SLIP_WALL_EVENTS_V2_2026-10-29/README.md) now passes native gapped-ID parity/rollback/collision tests and adds explicit endpoint checks. The preserved-ID Q1 IC has passed native initial-state validation and its H200 checkpoint chain has started; see [the Q1 experiment](../KEPLER_Q1_BOUNDED_GLOBAL_2026-10-29/README.md). The v1 evidence and earlier remaining-work list below are preserved.

This opt-in implementation extends the user's native CUDA+MPI+OpenMP solver. It passed a local replay of the original late Q2 wall crossing. It has not been promoted to a production executable and does not establish Q1 accuracy or long-time Kepler stability.

## Source and reproduction

`BUILD_V1.json` identifies the frozen executable, source hashes, build host and modules. Apply `NATIVE_EVENTS_V1.patch` to the exact first-RHS diagnostic source pinned by `../KEPLER_SLIP_WALL_NATIVE_2026-10-29/RHS_PROBE_BUILD.json`, checking the base hashes first. Use GNU `patch -p1` or `git apply --unidiff-zero`. Copy the two `Exam/` files in this directory alongside the previously published C-compatible `kepler_slip_wall.h`. The reconstructed files were verified byte-for-byte against the frozen candidate. The patch retains the legacy source encoding; do not reformat the base files.

The frozen binary is `/gpfs/kjhan/LagForce/Kepler_SlipWall_Integration_20261029/candidate_v1/eunha2.exe`, SHA256 `32fa3fe90c4c3a41126df1802c81ea382d6d8af2d4dbb3470c5c14f6bd0b9f35`. Large binaries and snapshots remain on GPFS. `runners/` preserves the effective test scripts with their original absolute paths; `RUNNER_MANIFEST.json` records their hashes. These scripts are execution receipts and must be staged into fresh directories before reuse.

Enable with `GFS_KEPLER_WALL_EVENTS=1`. V1 requires one MPI rank, contiguous IDs, W0, fixed masses and the tested native smooth pressure RHS. It rejects incompatible evolution modes. It introduces no generator steering, viscosity, intrinsic spin or thermal repair.

## Numerical behavior

The classical RK4 cubic continuous extension proposes the first outward polygon-plane crossing, including an out-and-back trajectory. Actual RK4 trials refine the event. A small, previously tested outward RHS extension supports trial predictors; no predictor is clamped or reflected. Only accepted contact receives the elastic normal impulse. The driver advances by the actual accepted prefix and recomputes the subsequent RHS/CFL. Unresolved corners, retry exhaustion and insufficient time progress fail closed.

Rollback restores the full saved particle payload, counts, timestep and diagnostic ledgers and invalidates GPU geometry caches. Pressure-wall impulse and torque are accumulated separately from collision impulse and torque. Static-wall collision work is zero; kinetic roundoff is reported without internal-energy compensation. The cubic proposal is third order; these tests do not establish an overall fourth-order event scheme.

## Evidence

- `UNIT_RESULTS_V1.json`: tests execute the actual event include for analytic flight, acceleration, grazing, out-and-back, repeated impacts, retries and corner rejection. MPI/geometry/cache adapters in these unit tests are test doubles; hardware tests are separate.
- `TRANSACTION_COMPARISON.json`: actual GPU flag-off parity, no-hit parity, forced retry and forced failure all pass. Comparisons preserve all stored record bytes except the runtime next pointer, plus the complete restart footer.
- `CROSSING_TO484.json`: the original accepted predecessor at t=48.37534965270793 advances through the former failure to exact t=48.4. All 18 audit checks pass. One collision occurs for ID55283, facet543, at t=48.37702401924956; no floors or energy injections are used.
- `HALF_DT.json`, `QUARTER_DT.json`, `CONVERGENCE.json`: all three native runs pass. Their collision times span 3.91e-13. Relative restart energy errors decrease from 1.48e-11 to 7.72e-13 to 3.65e-14. Endpoint field differences and collision impulses are reported without fitted pass tolerances. Thermal/pressure fields remain more sensitive than positions and velocities; this is a local convergence result.
- `WALL_TRACTION.json`: independent quadrature of actual wall faces at all four RK stages agrees with the fluid-force wall ledger. Absolute impulse/torque differences are at most 1.31e-29, below the declared tolerances. Each stage includes all 1024 wall facets. Captures leave the accepted restart body/footer unchanged. This check uses a tiny non-collision step and does not establish long-interval traction convergence.

The audit scripts are the three project-level `audit_kepler_wall_events_20261029.py`, `audit_kepler_wall_event_convergence_20261029.py` and `audit_kepler_wall_traction_20261029.py` files. They use the previously published native snapshot readers and comparison helpers.

## Remaining work

Before a clean Q1 bounded experiment, support preserved, noncontiguous global IDs and enforce an explicit final-domain invariant on every acceptance path. Validate both with analytic tests and the actual late GPU fixture. Preserve all canonical Q1 disk IDs, positions, velocities, masses and thermal state; initialize the changed exterior boundary only once and independently validate disk volumes, total polygon area and pressure. Then use the unchanged Q1 t=3 radial gate before continuing the same accepted chain to t=6.5 and 10.

The failed periodic Q1 t=3 checkpoint and older production binaries remain preserved. Boundary-radius/shape dependence, runtime overhead, GIZMO t=120/600 comparisons and the interpretation of the approximately 1000-orbit target remain unresolved. V1 backs up and scans the full particle set per step, so its cost has not been demonstrated to meet the requested production speed.
