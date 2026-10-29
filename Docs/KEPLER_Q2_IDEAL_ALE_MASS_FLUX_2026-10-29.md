# Mass exchange required by a uniform-density Kepler reference

For the actual saved bounded Q2 generator positions, take the ideal softened
Kepler velocity field and hypothetical spatially uniform density `rho=1`.
Using the current W0 Voronoi face velocity `w_n` and exact analytic fluid
velocity `u`, evaluate the usual material flux

`m_dot_i = - sum_faces integral rho * (u_n - w_n) ds`.

This is a **counterfactual**: the user's W=0 reduction requires fixed cell
masses, so applying this full flux to the W0 native solver would violate that
condition. It is used only to quantify how much exchange a uniform-density
Euler reference requires on the same instantaneous geometry. Opposite
signs are assigned to each interior face's two owners. Exact Kepler velocity
has zero continuum divergence; hence uniform density stays unchanged if the
exchange equals the geometric area change. Four-point face quadrature gives
the following mass-weighted rates for the 3,931 initial inner-disk cells,
normalized by each cell's mass and local orbital frequency:
the normalization uses the hypothetical mass `rho*A`, while the averages
are weighted by the stored native masses to keep the Q2 cohort comparable.

| Measure | t=30 | t≈47.5 |
| --- | ---: | ---: |
| Mean absolute net mass rate | 0.007868 | 0.019618 |
| Mean sum of absolute face exchange rates | 0.041902 | 0.061810 |
| Cell p99 absolute net mass rate | 0.024027 | 0.111282 |
| Mean absolute residual density rate after exchange | 2.64e-13 | 2.69e-13 |

The last row is a spatial-quadrature consistency check for this analytic
uniform-density field, not evidence of an evolved stable solution. The
two- versus four-point difference in the mass-weighted mean normalized net
rate is 2.31e-9 and 7.45e-9. The sum of absolute net rates is only about
19% and 32% of the sum of gross face rates, respectively: substantial
facewise exchange cancels inside each cell. These are instantaneous rates,
not accumulated mass changes or a timestep limit.

This quantifies the tradeoff behind the user's preferred fixed-mass limit.
The previously tested weight-speed-only exchange transports the mass swept
by the **change of weight**, but does not by itself cancel the baseline
Voronoi generator/material slip. A viable variable-mass extension would
need an explicitly weight-linked mechanism for the remaining exchange,
shared momentum and energy fluxes, spin/angular momentum accounting,
positivity and a defined reduction to exactly fixed masses when W=0. No such
coupled mechanism is claimed or installed here. Full Euler ALE exchange at
W0, pressure-only modification, or generator steering are not selected.

Input hashes, reconstruction checks, radial metrics and the full counterfactual
are recorded in the [receipt](KEPLER_Q2_IDEAL_ALE_MASS_FLUX_2026-10-29.json)
and its [reader](audit_kepler_q2_ideal_ale_mass_flux_20261029.py).
