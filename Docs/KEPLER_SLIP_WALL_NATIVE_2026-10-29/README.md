# Native slip-wall impulse primitive, 2026-10-29

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
