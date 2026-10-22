# Inner Kepler refinement: benefit and limits

2026-10-21. User asks whether more Laguerre generators toward the center can
reduce shear-related error. Conclusion: useful and already measurably helpful,
but insufficient by itself; preserve the original physics and assess local
thermal/work error as well as spatial reconstruction error.

## Shear scale for the actual softened orbit

The native smooth IC uses GM=1 and

    Omega(r) = (r^2 + epsilon^2)^(-3/4),
    |r dOmega/dr| = (3/2) r^2 (r^2 + epsilon^2)^(-7/4).

These follow directly by differentiating the native velocity formula. The
second expression measures differential rotation, with its factor-of-two
strain-tensor convention left explicit. Across radial cell size h, the
leading shear velocity difference in a local rotating frame is approximately
h |r dOmega/dr|. Reducing h reduces this difference and curved-flow
reconstruction error; it does not remove the continuum shear itself.

For r much greater than epsilon, this is (3/2) Omega h, and relative to orbital
speed its scale is (3/2) h/r. Constant h/r is therefore a reasonable geometric
resolution criterion outside the softened core. A cold-flow criterion should
also examine h |r dOmega/dr| / c_s, with c_s from the actual density and P.
Constant h/r and constant shear Mach number are different requirements.

For r much smaller than epsilon, Omega approaches a constant and the shear
goes to zero proportional to r^2. The shear rate itself peaks at
r = (2/sqrt(3)) epsilon, about .05427 for epsilon=.047. It does not increase
without bound toward the origin in this test. Refinement should therefore
follow measured error and the softening transition, rather than an unlimited
r->0 particle concentration.

## Existing controlled evidence

Previous same-epsilon initial-field refinement changed65536 to71080 cells
(+8.46%) and minimum spacing .09375 to .01171875 (8x finer). Total initial
integrated mass remained75.97420601059443. Initial r<.5 sum of positive native
internal-energy rates fell1.81448e-4 to1.86462e-7, about973x. This includes
reversible pressure work; it is not an entropy-production convergence order.

The refined native RK calculation still failed near t=.102. Subsequent strict
SSPRK3 reached t=.135091 without committed negative U, but four focus cells
lost thermal energy to U/U_adiabatic=.01947-.02108. See `RESULTS.md` and
`SSPRK3_RESULTS.md`. Refinement has not passed t=.2, t6-10 or paper duration.
The signed face ledger attributes that captured excessive cooling to
pressure-contrast work, with negligible work/geometric velocity mismatch.
Increasing resolution alone does not repair this discrete thermal contract.

Pure Euler P2/P4/P6 projection independently improves away from four central
cells but fails their cold pressure. The literal original rho formula has a
separate r+1e-6 feature, unrelated to gravity epsilon. See
`DG_PROJECTION_RESULTS.md`. This is a projection obstacle, not proof of the
historical instability's cause. Original IC remains unchanged.

## Implementation consequences

Increase generator count while preserving the same continuous rho,v,P,
potential, total mass and conserved initial integrals. Smaller inner cells
need smaller initial cell masses at the same rho; those masses may still be
fixed individually in time. Keeping one identical mass for every cell while
arbitrarily increasing their count changes the density/mass distribution.
Changing Laguerre weights alone does not create additional degrees of freedom.

Use gradual resolution transitions and record geometry quality, pressure and
angular/energy ledgers by radius. Refinement can reduce the allowable global
timestep and makes unresolved cell modes more expensive; the different
historical TCF/ring branch had a mode growing as c_s/h, so finer is not a
universal stability cure. Do not rename/repeat that old method.

Keep the already refined original layout for the next operator comparison.
First repair the conservative Euler spatial/thermal representation using
`ASTRA_ENRICHED_EULER_PLAN.md`, then compare coarse and refined continuations
through the historical failure interval with the same validated operator.
Initial projection or a shorter run is not the requested long-time result.
