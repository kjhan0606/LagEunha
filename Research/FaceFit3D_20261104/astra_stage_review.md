# Astra stage evaluation and next batch

The existing astra_dimensional_hydro_plan agent evaluated this stage once.
The review accepted the face-geometric weighting and strict geometry approach
within the measured scope. The review required bounded endpoint clamps, an
explicit double-precision contract, and a small 3D shock timestep refinement.
The final implementation and archived fixtures/actual runs satisfy those
requirements. The angular weight is a proxy, not an exact solid angle.

The next batch is Gresho, fixed-width smooth KH, a new-fit Kepler t0 RHS
audit, and removal of duplicate 2D geometry generation. Velocity-dependent
CFL, fit, viscosity coefficients, and active signs must still be refreshed.
Reuse requires unchanged positions, W, ownership/image mapping, and explicit
invalidation after restart or reordering. Contact/shear/shock regression and
standalone profiling precede a matched GIZMO speed claim. Sedov/Noh follow,
then RT/cylinder with separately defined force and boundary contracts.

Full physics validation, new-method long Kepler, and the speed goal remain
incomplete. Root performed the implementation; Astra performed evaluation
and planning only.
