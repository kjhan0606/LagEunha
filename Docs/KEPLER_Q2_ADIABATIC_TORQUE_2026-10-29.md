# Adiabatic versus entropy contributions to late Q2 pressure forces

For the saved bounded Q2 W0 states, reconstruct the existing variational
pressure force on the actual saved geometry twice. The first reconstruction
uses the stored pressure. The counterfactual keeps the same density, sites
and neighbor polygons but uses each site's initial entropy coefficient:

`P_ad = P_initial * (rho_actual / rho_initial)^gamma`,

with the native parameter gamma represented as float32 `1.666667`. The
difference of the two force reconstructions isolates the effect of the
saved pressure departure from that adiabatic relation **at that state**.
Both actual and counterfactual pressures are used for every needed neighbor,
not merely for the inner cells. This is not a separate simulation.

For the 3,931 original inner-disk cells (`0.5<=r0<0.75`), pressure
acceleration is summarized by mass-weighted vector RMS:

| Time | Actual acceleration RMS | Adiabatic counterfactual RMS | RMS difference / actual RMS |
| --- | ---: | ---: | ---: |
| 30 | 3.7285e-7 | 3.7305e-7 | 1.18% |
| 47.5003 | 1.04317e-5 | 1.04324e-5 | 0.0516% |

The signed pressure-torque projection onto the existing `0.2<=m/N<0.3`
angular-momentum-error band is `0.191843` actual versus `0.192003`
adiabatic at t30, and `0.204999` versus `0.205017` at t47.5 (per unit time).
Those are instantaneous projections, **not measured growth rates**. The
mass-weighted RMS of `log(P/P_ad)` is 1.73e-4 and 1.97e-4, respectively;
its signed mass means are negative, -1.06e-4 and -1.24e-4. Initial pressure,
initial/current density and final pressure are stored float fields, so small
entropy terms include their representation and integration error.

This evidence says the late pressure feedback in this particular aligned Q2
state is overwhelmingly reproduced by **density-dependent adiabatic pressure**.
It does not support describing the measured late force primarily as entropy
production or net irreversible numerical heating. The reversible density
variation itself may be physical, mesh-induced, or a coupled error; the
saved-state comparison does not determine its origin. The geometry, mode
selection and pressure have evolved together, and no conclusion follows
about how the trajectory would behave if entropy were constrained during an
independent run. Inner core, other ICs, W!=0, spin, t600 and physical heating
rates remain untested by this diagnostic.

All endpoint hashes, five bands, values and geometry-volume checks are in
the [receipt](KEPLER_Q2_ADIABATIC_TORQUE_2026-10-29.json); the
[reader](audit_kepler_q2_adiabatic_torque_20261029.py) reconstructs the
same W0 pressure formula already checked against three native t47.5 pairs.
