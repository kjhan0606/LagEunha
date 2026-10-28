# Centroid-aware generator velocity: narrow actual-geometry check

This follows the formulation already recorded in
`staggered_affine_moment_closure_2026-10-09.md`; it is not a new velocity law.
New evidence uses the current canonical initial common-face geometry.

For affine v(x)=a+Gx, mass m, first spatial mass moment q, and integrated
momentum p, let c=q/m. The reconstructed point velocity is

    v(X) = p/m + G(X-c).

Using this as generator velocity is pointwise reconstruction, not a force or
centroid steering term. q evolves as qdot=p+integral_face x rho(s-v.n), including
the volume advection term p. A face-only q equation would be incomplete.

For manufactured steady simple shear v=(y,0), G²=0. With zero net mass flux,
the transport part of generator acceleration is

    [pdot - G qdot_transport]/m.

The omitted volume term and generator-motion term vanish after multiplication
by G for this particular field. Exact affine line integrals on actual t0 faces
give zero recovered acceleration in the reported floating-point calculation,
although pdot/m has inner RMS8.7513e-4. Independent two-point Gaussian
integration matches the analytic momentum face integral within7.28e-20.

Evidence: `canonical_affine_generator_recovery_2026-10-18.json`.
Script: `audit_affine_generator_recovery.py`.

## Limits

- G is prescribed exactly, not estimated from evolved cell averages.
- Uniform density and simple shear are manufactured on the saved geometry;
  this does not convert the canonical initial mass/momentum state.
- The cancellation is an algebraic consistency check. No time integration,
  shock, gravity, positivity, angular-momentum, MPI, or long-time test passes
  on its basis.
- Momentum-only closure is insufficient for energy: unresolved velocity
  variance must not be mistaken for thermal energy, as earlier remap checks
  already demonstrated.
- General time-dependent/nonlinear fields require Gdot and consistent moment
  evolution, not the simplified G²=0 cancellation used here.
- The original fixed-mass Voronoi reduction remains a separate gate; no
  candidate has been adopted in production.

Next relevant closure check is estimation of G from cell-integrated data at
actual mass centroids, rather than assuming the analytic shear gradient.

## Gradient estimation follow-up completed

`canonical_centroid_gradient_recovery_2026-10-18.json` fits gradients from
manufactured affine cell averages on the actual initial geometry. The fitting
algorithm receives average differences, centroid differences and face weights;
it is not supplied the exact gradient. Simple shear, rigid rotation and a
general trace-free affine field were checked for r<1.9.

The mass-centroid fit has maximum gradient error <=1.17e-14 and reconstructed
generator velocity error <=8.96e-16. Treating those same averages as samples at
generator positions gives maximum gradient error .0466–.0623. This comparison
concerns a proposed cell-average model, not a demonstrated gradient bug in the
original point-state code.

An initial centroid audit rejected common-intersection faces as an exact cell
polygon: their reconstructed area differed by up to5.98e-10 relative. Original
closed polygon edges reproduce stored area to4.45e-16 and were used instead.
Uniform density makes their geometric centroids equal to mass centroids here.
Nonuniform density, evolved gradients, limiters, energy and nonlinear-time
closure remain unverified. No restart fields were reinterpreted or overwritten.
