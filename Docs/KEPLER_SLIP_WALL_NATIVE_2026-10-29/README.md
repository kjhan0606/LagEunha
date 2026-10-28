# Native slip-wall impulse primitive, 2026-10-29

Later update: the opt-in native RK4 event controller and its local GPU replay, rollback, timestep and independent wall-traction checks are now recorded in [the frozen v1 integration bundle](../KEPLER_SLIP_WALL_EVENTS_NATIVE_2026-10-29/README.md). The remaining-integration list below preserves the earlier planning state; Q1 and long-time validation remain pending.

This is a new opt-in component for the user's existing CUDA+MPI+OpenMP bounded Kepler solver. It is not installed into a production executable. The published fixed-polygon source overlay and the failed periodic t=3 checkpoint remain unchanged.

`Exam/kepler_slip_wall.h` is shared host/device C++ arithmetic. It accepts an already located, committed impact at a stationary boundary and returns reflected velocity, wall impulse, wall torque impulse about the supplied center, and kinetic roundoff. It neither moves a generator nor changes its mass, thermal energy or weights. It adds no steering, viscosity or energy floor. Pressure traction/work remains a separate continuous RHS ledger.

The plane convention matches the existing `fixed_polygon_postStage_blend`: outward unit normal at angle `(side+.5)*2*pi/N`, plane offset `R*cos(pi/N)`. Use this primitive with the existing polygon geometry; its wall torque must be retained. The circle primitive is a mathematical reference and can only be used after the tessellation clips to a true circular boundary. Mixing radial reflection with polygon clipping is excluded. Simultaneous corner impacts require a separate event policy and are not resolved by this primitive.

Host tests pass, and 4096 actual RTX 5000 Ada kernel cases pass host/device tolerance checks, normal/tangential impulse identities, kinetic conservation and circular center-angular-momentum conservation. A separate planar test verifies fluid angular-momentum change plus wall torque closes. Invalid/off-boundary inputs leave the output unchanged. `VALIDATION.json` records source hashes, compiler, GPU UUID and measured residuals. No production speed or Kepler stability is inferred from a tiny primitive kernel.

## Remaining integration bundle

1. Keep the current post-stage wall crossing abort until a real event controller is tested. Do not replace it with predictor clamping or per-stage reflection.
2. Add accepted-step event location and rollback covering particle/RK state, GPU caches and diagnostic ledgers. Commit the boundary impulse once and evolve the remaining physical interval. Reject unresolved simultaneous corner impacts. Only inspect wall-near cells/faces for events.
3. Accumulate pressure traction and collision impulse/torque separately; include wall impulse/torque in conservation reports, with zero work for a static elastic wall. Do not repair kinetic roundoff through internal energy.
4. Run uniform-pressure bounded-force/work and a short reproduction from the original Q2 wall-crossing predecessor. A successful Q2 replay validates the collision path and does not substitute for the Q1 test.
5. Build a clean boundary IC preserving Q1 disk IDs, positions, masses and thermodynamic state, with an adequate exterior buffer. Validate initial disk volumes; do not delete exterior sites from the contaminated t=3 state and call that a restart.
6. Test native Q1 t=3,6.5,10 using unchanged radial and zero-injection gates. Then establish boundary-radius/polygon convergence for t=120,600 and the interpreted approximately 1000 inner-orbit target. Closed-wall success is not automatically GIZMO boundary equivalence; meaningful wall dependence requires a genuine free/outflow branch.

Astra recommended this bundle after the saved-state exterior-support attribution. Full long-time stability remains unverified.

## First-RHS extension diagnostic added

`WALL_RHS_PROBE.patch` reconstructs the limited native diagnostic from the tested fast-CPU/selective-wall source. `RHS_PROBE_BUILD.json` pins all source, linked archive and executable hashes. It permits a bounded outward generator excursion (at most 1e-6) only with `GFS_WALL_RHS_PROBE_ONLY=1`, computes the first native K1 RHS for all owned particles and exits before any predictor or accepted evolution. This is not the RK4 event controller or a production boundary switch. Normal flag-off output matched the original entire restart body (excluding the next pointer) and full footer exactly in a tiny accepted replay from the late Q2 fixture.

The preimpact fixture is now captured at exact t=48.37534965270793; the unchanged native solver reproduces the original ID55283/side543 failure immediately after it. Inside/contact/outside/max-band probes pass independent volume/force/work and native selected-face support/reciprocity checks; band-exceeded and empty-cell fixtures reject. See the project-level probe result receipts. This validates a local trial RHS extension, not collision accuracy, rollback correctness or long-time stability.

`EVENT_CONTROLLER_DESIGN.md` records the next integration points and unimplemented tests. The C-compatible header now passes C, C++ and actual Ada CUDA compilation/runtime checks; `VALIDATION.json` preserves the initial C++-only receipt alongside the updated header validation.
