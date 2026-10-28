# Native fixed-mass Laguerre q2 Kepler at t=10 (2026-10-27)

The isolated v2 CUDA binary, SHA-256 `5214b5b3a68f3a309e8650de918ca906aaac7d8f41bfcfd9f07ef8f9c4a55342`,
continued the two accepted finite-weight `t=6.501` checkpoints without
reinitialization. Slurm 437384 (aligned) reached exact `t=10` and passed
19/19 independent [endpoint checks](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_weighted_t10_a10_v1/aligned/audit_t10.json).
Slurm 437385 (independent 5% site jitter) also reached exact `t=10`, but
its [endpoint audit](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_weighted_t10_a10_v1/jitter/audit_t10.json)
passed 18/19: only the chosen disk mean density-change `<1e-3` gate failed.
The failed Slurm exit records this physics gate, not a simulation crash.

The [hash-verified same-time comparison](KEPLER_Q2_WEIGHTED_NATIVE_V2_VS_W0_T10_2026-10-27.json)
uses the identical `t=0` IC within each pair and the original W=0 native
continuation as control. The [common-pixel density image](KEPLER_Q2_JITTER_W0_WEIGHTED_DENSITY_T10_CLEAR_2026-10-27.png)
shows the independent-jitter pair.

| IC and solver | Disk mean absolute `Σ/Σ₀−1` | Mean absolute entropy-log change | Mean orbit-position error |
|---|---:|---:|---:|
| aligned W=0 | 1.96891e-4 | 1.19396e-5 | 3.57456e-7 |
| aligned weighted | 1.45954e-4 | 9.06615e-6 | 3.35180e-7 |
| independent jitter W=0 | 2.37446e-3 | 2.03302e-5 | 1.62485e-5 |
| independent jitter weighted | 1.71685e-3 | 1.53116e-5 | 1.59889e-5 |

The weighted jitter density measure improves by a factor of 1.383 (27.7%)
but still exceeds `1e-3`. Its orbit-position error improves by only 1.6%.
The original W=0 jitter run failed a *different* gate, its native-minus-
ballistic Voronoi density residual; that exact-orbit W=0 area target is not
used as a weighted-cell gate.

Both weighted endpoints have exact fixed masses, positive states and power
areas, no floor/dual-energy/entropy-floor injection, whole-energy relative
error near `8.1e-9`, independent disk power-area relative error below
`6e-8`, and pressure-weight power below `1.2e-26` in every recorded RK4
stage. The weight ranges at `t=10` are `1.68e-6` aligned and `3.79e-6`
jitter. Their 7,216 stages each spent mean `2.91–2.92 ms` computing the
rate; the complete `t=6.501→10` jobs took about 48 minutes on one A10 each.

This confirms a modest short-time density benefit and a working bounded
native weight path. It does not establish local moving-face mass/energy
consistency or the long GIZMO `t=600` result. The initial one-step 15–22×
inner-disk density advantage did not persist.

The 64-iteration Jacobi continuation from the same audited `t=7.426` state
has now finished both ICs. [Aligned](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_weighted_jacobi64_from7426_to10_v1/aligned/audit_t10.json)
passed 19/19; [jitter](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_weighted_jacobi64_from7426_to10_v1/jitter/audit_t10.json)
reached exact `t=10` and passed 18/19. Its sole failure remains the chosen
disk-density gate: `1.710128e-3` versus `1.716845e-3` with 16 iterations,
an improvement of only **0.39%**. The aligned density measure changes from
`1.45954e-4` to `1.44529e-4` (**0.98%** improvement). The 64-iteration
jitter energy, power-area, zero-injection, stage-work and fixed-mass checks
all pass; Slurm **437394 FAILED 1:0** represents the audit gate, not a
runtime failure. Increasing Jacobi iterations alone does not account for
the remaining independent-jitter drift.

The separate ring-coherent IC control tests whether within-ring disorder is
the larger cause. If these controls do not clear the density gate, the
[fixed-mass compatibility limit](KEPLER_WEIGHTED_FIXED_MASS_COMPATIBILITY_2026-10-27.md)
points to a conservative moving-face flux or another explicit weight-work
reservoir as the next physical closure, rather than an arbitrary
initial-volume feedback.
