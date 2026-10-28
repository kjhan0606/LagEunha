# Weight-only swept transport does not cancel the initial shear source

## Narrow statement and derivation

At a fixed W=0 geometry let Bv denote the geometric volume rate induced by
generator motion, and let L Wdot denote the additional rate due to relative
power-weight motion. For the convention |x-X_i|²-W_i, each directed face gives
L_ij = -A_ij/(2d_ij), with positive diagonal sum. Then

    Vdot = Bv + L Wdot.

Consider a locally uniform density and a mass-relaxation prescription that
transports only the mass swept by changing weights:

    mdot = rho L Wdot.

Since rho=m/V, differentiating yields

    rhodot/rho = mdot/m - Vdot/V = -Bv/V.

Thus this specific physically motivated remapping term cancels its own
geometric volume change and leaves the generator-induced density source
unchanged at that state, regardless of the chosen weight rate. This is not a
proof that every weight-coupled method fails. Later geometry changes,
nonuniform reconstruction, physical transport, and other closures are outside
this instantaneous uniform-density statement.

An illustrative interpolation mdot=alpha rho L Wdot gives

    rhodot/rho = -[Bv + (1-alpha)L Wdot]/V.

For alpha<1, cancelling this source requires L Wdot=-Bv/(1-alpha). As alpha
approaches 1, the required volume correction increases rather than decreases.
At alpha=1 cancellation is impossible unless Bv=0. Arbitrarily scaling swept
mass by alpha is not a justified physical flux law and is not adopted here.

## Numerical check and limits

The script `audit_weight_only_density_response.py` uses saved full-precision
t0 common faces and original circular velocities with a manufactured uniform
density, not the actual nonuniform masses. Four amplitudes of a smooth weight
rate check the algebra on the actual geometry. No physical state is changed.
Results: `canonical_weight_only_density_response_2026-10-18.json`.

This explains why successful component remapping/conservation tests do not
establish a shear-stability cure. They remain useful building blocks, but
merely coupling mass changes to weight-swept volumes does not remove the
measured initial source.

## Consequences for the next coupled candidate

1. Identify the additional mechanism that modifies the generator-induced
   material/geometric mismatch; specify its shared-face mass, momentum and
   energy balance before evolving a long run.
2. Retain the approved reduction when weight differences remain zero. A
   generic ALE flux that persists at fixed equal weights does not meet it.
3. Do not recover the reduction by suppressing physical face transport without
   a model for the omitted transport. Do not choose alpha by a stability fit.
4. Existing volume-following DAE and weight-only remapping cannot simply be
   combined and assumed to avoid secular weights. The identity above is a
   necessary local check before any such claim.

No production correction is enabled by this audit; long-time physical
stability remains unachieved.
