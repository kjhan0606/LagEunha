# Exact moving-contact controls across all certified segments

For each certified donor pair, manufacture equal pressure and normal velocity
using arithmetic averages of the t2 values, retaining donor densities and
tangential velocities. These are controlled contacts, not a Riemann solution
for the actual unequal states. Assume gamma5/3 and constant primitives within
each segment. The earlier fitted density is not integrated by this control.

With contact normal speed u_n and geometric face speed w_n(s), select the
upwind contact state by sign(u_n-w_n). The laboratory-frame ALE flux is

    mass = rho*(u_n-w_n),
    momentum = rho*v*(u_n-w_n)+P*n,
    total energy = E*(u_n-w_n)+P*u_n.

Split the segment where the affine relative speed changes sign before
integration. There are2540/2524 such crossings among5371/5365 segments for
weight fractions .25/1. A single midpoint upwind state would miss many of them.

Both cases pass:

- Analytic integration vs split two-point Gauss: max3.80e-19 absolute.
- Independently reversed normal, side ordering and parameterization: max1.31e-17.
- Uniform boost(.7,-.4), transforming face speed and flux consistently: max4.00e-17.
- Comoving contact w_n=u_n: mass flux zero to rounding, momentum P*n and
  energy P*u_n retained; max error2.07e-17.

The summed absolute contact pressure-work contribution is about.001764 in
these diagnostic flux units. It is not a global energy change. Omitting
pressure work when mass flux vanishes would fail even this simple control.

## Reduction and remaining physics

Zero mass flux follows from a comoving material contact, not merely from zero
weights. The geometric face speed and contact speed generally differ. Therefore
this control does not prove the approved original fixed-mass reduction.
Do not force u_n=w_n or suppress real advective flux to claim that reduction.

Next use the existing solver's general interface implementation as an object
of comparison against these exact-contact identities, before introducing a
new solver. General jumps, shocks/rarefactions, residual velocity structure,
full moment flux and work, original-model compatibility, backend/MPI parity
and long physical evolution remain unverified.

No production state or original Kepler trajectory was changed.
Evidence: `canonical_moving_contact_flux_2026-10-18.json`.
Script: `audit_moving_contact_flux.py`.
Remote: checkpoint-remap directory, `moving_contact_flux.json`.
