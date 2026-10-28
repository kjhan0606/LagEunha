# Relaxed weights with full residual transport: conditional candidate

At uniform incompressible density, write Vdot=Bv+L Wdot. Instead of enforcing
Vdot=0 indefinitely, consider the diagnostic relation

    L Wdot = -Bv - L W/tau.

For frozen geometry and forcing, with outer weights pinned and u=L^-1 Bv,

    W(t) = -tau u + (W(0)+tau u) exp(-t/tau).

This has bounded weights in the frozen model. The remaining volume rate is
-L W/tau. Full ALE mass balance at uniform density requires

    mdot = rho (Bv+L Wdot) = -rho L W/tau.

This is distinct from weight-swept-only transport rho L Wdot, which leaves the
generator-induced density source unchanged. At frozen stationary W=-tau u,
Wdot=0 but mdot=rho Bv remains nonzero.

## Actual-geometry tangent calculation

Using full-precision t0 common faces, the r<1.9 Dirichlet solve has relative
residual5.24e-16. For tau=.01,.1,1, the predicted maximum adjacent weight
difference/d² is .0041063,.041063,.410630. These are linear tangent predictions,
not rebuilt nonlinear-cell positivity or long-time results. The manufactured
uniform-density mass identity agrees within4.21e-19 in volume-rate units.
Data: `canonical_relaxed_weight_tangent_2026-10-18.json`.

## Unresolved consistency conditions

- A finite frozen-geometry stationary weight is not a nonlinear bound when
  L, topology, forcing, density and pressure evolve. No tau is selected for
  production and these values are not stability-tuned parameters.
- At a density interface, volume balance is not mass balance. Use actual
  reconstructed face states and preserve shared antisymmetric m,p,E transport.
- Zero net mass does not eliminate face circulation or its momentum/energy
  fluxes. Avoid donor mixing that converts shear kinetic energy to heat.
- Persistently equal weights must recover the original fixed-mass model under
  the approved reduction. A generic full ALE implementation does not prove
  that reduction. At W=0 with Bv nonzero this candidate immediately develops
  relative weight rates; equality at one instant is not the reduction case.
- Stationary nonzero relative weights with continuing mass transfer require
  explicit interpretation of the user's weight-linked relaxation condition.
  Do not silently strengthen or weaken that condition.
- Work, gravity, angular momentum, positivity, restart state and all RK stages
  must be consistently defined before a coupled evolution.

The candidate remains a mathematical possibility for further assessment, not
an adopted change or an established cure. Original runs and code are retained.

## Reduction and face-moment follow-up

At fixed finite tau, mdot=-rho LW/tau tends to zero with relative W. Sending
tau to zero along stationary W=-tau L^-1 Bv is a different, singular parameter
path: mdot remains rho Bv. The latter observation does not by itself violate
the approved condition for persistently equal weights at fixed parameters.
Do not silently replace that approved condition by a stronger requirement.
The complete original-solver reduction is still unproved.

`canonical_affine_weight_transport_2026-10-18.json` checks unit simple shear
v=(y,0), uniform rho=1, on actual t0 faces. At W=0 choose LWdot=-Bv in r<1.9.
The net mass-rate RMS is below1.6e-14, but cell-integrated momentum transport
per mass RMS is .00087513 in the inner region. A face midpoint evaluation
misses a contribution with RMS .00054453. These are diagnostic transport rates,
not measured physical acceleration errors or Kepler growth rates.

For a straight face coordinate z in [-A/2,A/2], let relative normal speed be
delta0+delta1*z and continuous affine velocity be vf+gt*z. Exact integrals are

    mass flux = rho A delta0,
    momentum flux = rho [A delta0 vf + A^3 delta1 gt/12].

Thus equal net cell mass does not remove momentum exchange, even in perfectly
smooth shear. The midpoint omits the second term. In this simple shear,
(v.grad)v=0: nonzero change of cell-averaged velocity can arise from moving
cell centroids, without acceleration of a material fluid element. Generator
velocity therefore cannot automatically be identified with p/m for this
cell-integrated candidate. This is a specific first-moment/face-quadrature
requirement; it does not justify adding expensive internal gravity quadrature
to the original solver or prove that such a change cures pairing.
