# Energy requirement for the proposed cell-average transport model

The manufactured simple-shear diagnostic uses actual canonical t0 faces,
uniform density, zero current weights and nonzero relative weight rates chosen
to cancel net cell-volume rates. It does not evolve or reinterpret the existing
Kepler states. Data: `canonical_affine_energy_transport_2026-10-18.json`.

For v(z)=vf+gt*z and relative normal speed delta0+delta1*z on a face of length A,
the exact kinetic energy transport is

    A delta0 |vf|²/2
    + A³/12 [delta1 vf.gt + delta0 |gt|²/2].

Two-point Gaussian quadrature independently reproduces this cubic-integrand
formula within1.00e-19 on faces of the tested interior cells.

The actual average kinetic-energy derivative is

    Kmean_dot = (p/m).pdot - |p/m|² mdot/2.

Subtracting this alone from total-energy transport misclassifies the difference
Kdot-Kmean_dot as thermal energy. Here its volume-normalized RMS is1.42e-5 in
the inner region and1.48e-5 in the core, with both positive and negative values.
It represents internal kinetic redistribution; it is not guaranteed positive
entropy production and cannot simply be put into an AV heating term.

This supports the earlier 2026-10-09 moment energy decomposition and the finite
remap evidence. It does not establish that original Kepler heating has this
cause: the original point-state model and this proposed cell-average model
have different state interpretations. No finite-time, positivity, shock,
gravity, parallel-backend or long-Kepler result follows from this audit.

The next coupled candidate must evolve enough kinetic information to recover
U without this misclassification, while meeting the original-model reduction.
The current experiments do not yet supply that complete nonlinear closure.

## Existing runs observed 2026-10-17 17:26:34 UTC

Actual live process checks found baseline PID2286374 at t15.40140,
constrained-HG PID2295268 at t15.35438 and from-t0 timestep-cap PID2466162 at
t10.80151. All accepted records had finite=1 and retries=0. These are execution
status facts, not physical-accuracy passes; prior t10/t14 errors remain.
