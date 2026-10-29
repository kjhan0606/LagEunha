# Geometric area rate versus material divergence in late Q2

The saved bounded Q2 W0 states at t=30 and t=47.50031653561342 were
reconstructed independently. For each cell, the instantaneous Voronoi area
rate follows the exact generator-driven face speeds. With fixed cell mass,
this imposes `d ln rho/dt = -d ln A/dt`. As a diagnostic comparison, the
saved velocity-gradient trace estimates material divergence under a local
affine reconstruction. The stored gradient may belong to a nearby RK stage;
it is not an independent exact continuum measurement.

An additional counterfactual keeps the **actual saved generator positions**
but substitutes the exact softened circular Kepler velocity at every site.
The analytic continuum velocity has zero divergence; nevertheless, the
Voronoi areas can change as those generators move. The exact-velocity face
area rate was checked against independently rebuilt polygons at generator
positions `x +/- 1e-5 v_K`; its relative RMS difference on the inner cohort
was 2.13e-6 at t30 and 3.77e-7 at t47.5.

Mass-weighted mean absolute rates normalized by local circular frequency:

| Inner disk, 3,931 cells | t=30 | t≈47.5 |
| --- | ---: | ---: |
| `|A_dot/A|` using saved generator velocities | 0.007868 | 0.019636 |
| `|div u|` from stored reconstructed gradient | 0.000292 | 0.000534 |
| `|A_dot/A - div u|` | 0.007705 | 0.019438 |
| `|A_dot/A|` using exact Kepler generator velocities | 0.007868 | 0.019618 |

The near equality of the first and last rows means the *magnitude at these
saved configurations* is not mainly removed by resetting the generator
velocities to the analytic orbit values. The evolving **positions and
geometry** still matter: this does not establish that the initial undisturbed
mesh has the same late rate. Absolute rates are snapshots, not integrated
errors or entropy generation, and their small difference is not a bound on
cellwise differences. Saved-state geometry matches native volumes to better
than 6e-8 relative. The earlier t0 diagnostic independently found the same
class of material/geometric face-speed mismatch.

Together with the separate adiabatic-pressure decomposition, these results
make the fixed-mass area-to-density coupling a concrete target for a
conservative closure. Simply replacing geometric work by `div u`, or adding
the computed difference as heat, would break the native force/work relation
without a matching derivation. The user-approved W=0 fixed mass and no
steering conditions still apply. Variable weights/masses and independently
transported subcell moments require their own coupled equations and tests;
no particular correction or long-time stability is established here.

Complete radial values and hashes are in the [receipt](KEPLER_Q2_VOLUME_DIVERGENCE_2026-10-29.json),
with reproducible [reader](audit_kepler_q2_volume_divergence_20261029.py).
