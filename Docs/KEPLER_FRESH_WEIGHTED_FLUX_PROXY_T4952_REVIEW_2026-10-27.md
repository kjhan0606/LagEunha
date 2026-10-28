# Moving-weight Kepler face mass-flux proxy at t≈4.952

The [hash-checked, accepted one-step audit](KEPLER_FRESH_WEIGHTED_FLUX_PROXY_T4952_2026-10-27.json)
uses the existing fresh-weighted native snapshots at `4.951043599→4.952043599`.
It reconstructs each disk-disk power face at the generator midpoint and
estimates `rho_face * face_length * (v_Kepler−w_face)·normal`, with arithmetic
neighbor density, analytic softened Kepler velocity at the face midpoint,
and secant weight rate. W=0 and moving-power comparisons use **the same
generator positions, velocities, densities and fixed cell masses** to isolate
the geometric effect of weights. This is a frozen-state proxy; it is not a
native Riemann flux or a second W=0 trajectory.

| One-step disk-disk proxy | W=0 same state | Moving power |
|---|---:|---:|
| Gross absolute face transfer / disk mass | `1.63369e-5` | `1.48449e-5` |
| Median `|net cell transfer|/cell mass`, full disk | `4.89924e-6` | `1.06635e-7` |
| p99 `|net cell transfer|/cell mass`, full disk | `5.48395e-5` | `3.26798e-6` |
| p99, inner `0.5≤r0<0.65` | `1.05423e-4` | `3.86078e-5` |
| Maximum, inner `0.5≤r0<0.65` | `1.49142e-4` | `7.05343e-5` |

Moving weights reduce the full-disk p99 **net-cell** proxy by **16.8×**,
but gross face exchange by only **9.13%**. The distinction agrees with the
[independent face-speed audit](KEPLER_FRESH_WEIGHTED_MATERIAL_FACE_T4952_REVIEW_2026-10-27.md):
volume compensation can cancel face fluxes around a cell without making
each face material. The apparent weak p99 gain of **2.73×** in the innermost
disk bin is caused largely by omitted disk-core faces, as checked below.

The first calculation reconstructed only disk-disk faces. Its inner-bin net
flux omits faces at the core boundary, so it must not be read as an actual
cell mass source. The code holds masses fixed; the proxy measures the gap
between that representation and hypothetical material transport through
geometric faces.

The [split by the first three disk rings](KEPLER_FRESH_WEIGHTED_FLUX_PROXY_INNER_RINGS_T4952_2026-10-27.json)
shows why the core-interface omission matters. The weighted p99 net-cell
proxy is `5.06501e-5` in `0.5≤r0<0.52` (294 cells), but only `9.38356e-7`
in `0.52≤r0<0.55` and `8.55733e-7` in `0.55≤r0<0.65`.
The first-ring proxy is therefore **not evidence of a true mass leak**.
An independent [native orbital angular-momentum audit](KEPLER_Q2_FRESH_WEIGHTED_ANGULAR_MOMENTUM_INNER_RINGS_T65_2026-10-27.json)
also finds its largest weighted inner-ring error on that first ring, so
the disk/core native face ledger, pressure torque and weight-mask interface
are the next targeted checks.

## Boundary-face correction

The [all-finite-face rerun](KEPLER_FRESH_WEIGHTED_FLUX_PROXY_ALL_FACES_T4952_2026-10-27.json)
adds **575 disk-core** and **2,311 disk-exterior** faces to the same power
diagram. It reproduces every previous disk-disk gross rate and then sums
the missing boundary contributions into each disk cell's net rate.

| One-step all-face proxy | W=0 same state | Moving power |
|---|---:|---:|
| Gross absolute face transfer / disk mass | `1.65675e-5` | `1.50582e-5` |
| p99 `|net cell transfer|/cell mass`, full disk | `5.59352e-5` | `5.64640e-7` |
| p99, first ring `0.5≤r0<0.52` | `1.26555e-4` | `2.28529e-6` |

The weighted full-disk p99 is **99.1×** below same-state W=0, while gross
exchange is only **9.11%** lower. On the first ring, the weighted p99 is
**55.4×** below W=0; the previous disk-disk-only p99 overestimated its
net-face proxy by **22.2×**. Thus the alleged large first-ring *mass-flux*
residual was an incomplete-face accounting artifact. The first-ring
*orbital angular-momentum error* in the native t=6.5 snapshot remains real
and requires separate force/torque attribution.

Both calculations still use arithmetic face density, analytic Kepler fluid
velocity, and secant weight rate. They are not measured native Riemann fluxes
or stage-exact mass conservation tests, and small net-cell flux can coexist
with substantial opposite face flows. This result narrows the next check to
the native disk/core face force and angular-momentum ledger rather than a
putative missing core mass flux.
