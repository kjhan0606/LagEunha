# Kepler q2 finite-weight volume-feedback screen (2026-10-27)

The [read-only probe](KEPLER_Q2_WEIGHTED_TARGET_FEEDBACK_T7426_2026-10-27.json)
reconstructed the actual finite-weight power geometry at the accepted native
weighted `t=7.426299367` checkpoint. It kept the site velocity and one native
step length (`dt=0.00193983`) fixed, solved the same 16-step Jacobi volume
equation, and projected every candidate weight rate onto zero global
pressure-weight work. It then rebuilt power areas after one Euler geometry
step. This is a screen, not a native RK4 continuation.

| IC | Current disk mean absolute density-equivalent error | Zero-volume-rate next step | Initial-volume feedback tau=1 next step | Feedback tau=0.2 next step |
|---|---:|---:|---:|---:|
| aligned | 5.70310e-5 | 5.71056e-5 | 5.70770e-5 | 5.69633e-5 |
| independent 5% jitter | 1.09471e-3 | 1.09541e-3 | 1.09541e-3 | 1.09546e-3 |

The stronger feedback does **not** reduce the jittered error even in this
single-step screen. The signed pressure-weight power remains numerically zero
for all candidates, so the lack of improvement is a geometry outcome rather
than a failure of the global projection. For the aligned IC, the tau=0.2
change is only 0.12% of the current error in one step.

The initial-volume target is specific to a circular Kepler equilibrium. A
general hydrodynamic flow may genuinely compress, so such feedback is not a
general fixed-mass Laguerre law. It is not being attached to the native mode.
The next native comparison instead varies only the linear-solver iteration
count on the same `t=7.426` state. The [Jacobi residual screen](KEPLER_Q2_WEIGHTED_JACOBI_RESIDUAL_T7426_2026-10-27.json)
predicts jittered one-step area-change p99 of `7.25e-6` for 16 iterations and
`3.06e-6` for 64, at unchanged global pressure-weight power. Whether that
reduces the evolved `t=10` density error requires the native branch result.
