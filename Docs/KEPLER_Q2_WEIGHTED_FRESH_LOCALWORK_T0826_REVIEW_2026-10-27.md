# Fresh weighted Kepler: local work hidden by global projection

An isolated one-step continuation used the 21/21 passing fresh weighted
`t=0.825` endpoint SHA-256
`b929fc51043675f0bdb09eb1858dcf374cd0dbc8f9bac42e7dffc6bd12d7a6dc`
and the same repaired v5 CUDA binary SHA-256
`feae0a48f737a7eab215f2e31be1d7938f45f894376dcdfda177f7e09f57d26b`.
It ran on the idle lageunha Ada to exact `t=0.826` in 3.625 seconds;
its [native endpoint audit](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_capacity64_v5_fresh_localwork_t0826_ada_v1/audit_t0826.json)
passes **21/21** checks, with all four weighted RK4 stages active,
independent power-area error `5.95e-8`, fixed masses, and zero recorded
injection. The endpoint SHA-256 is
`84a6e2bd8e4ea400fdc13e729326552453c84530192359d1ac905e1c22727f57`.
The isolated [stage manifest](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_capacity64_v5_fresh_localwork_t0826_ada_v1/stage_manifest.json)
and runner hashes preserve the parent state and inputs; this diagnostic
did not alter the still-running longer weighted chain.

The [independent local-work reconstruction](KEPLER_Q2_WEIGHTED_FRESH_LOCALWORK_T0826_2026-10-27.json)
uses the midpoint positions, weights, pressure, and a secant
`Wdot=(W_end-W_start)/0.001`. Each shared power face contributes
`g_ij(Wdot_i-Wdot_j)` to cell i's weight-driven volume rate, where
`g_ij=L_ij/(2d_ij)`. Its pressure power is `P_i Vdot_i^(W)`.
The cell sum and face sum agree within `1.33e-25` absolute power.

| One-step midpoint estimate | Value |
|---|---:|
| 99th percentile of absolute cell weight work / cell internal energy | `3.706e-5` |
| Maximum of that ratio | `1.041e-4` |
| RMS relative volume change from weights alone | `1.403e-5` |
| RMS observed total relative volume change | `8.626e-7` |
| Absolute face work / total internal energy | `7.735e-11` |
| Native maximum projected global stage work | `3.659e-29` |

Weights counteract a larger geometric volume change; the observed total
volume change is about 16× smaller in RMS than the weight component alone.
The nonzero local pressure-weight work largely cancels in the domain sum,
while the native GPU code projects only that global sum toward zero. This
is a **term requiring local accounting**, not by itself a measured violation
of the first law. The midpoint power area differs from the mean endpoint
volume by at most `2.90e-6` relative; the secant rate and midpoint geometry
are not the exact four RK4 stage states. A stage-level pairwise work and
face-transport ledger is still required before promoting the weighted
operator to a physically conservative method or extrapolating its excellent
early density result to the GIZMO longest-time test.
