# Bounded q2 aligned W=0 Kepler reached t=30, with growing inner error

Slurm **437348** completed `0:0` in `1:44:42` on one A10. It continued
the accepted `t=20` aligned state with the bitwise-checked selective-wall
CUDA binary to exact `t=30`. The
[native endpoint receipt](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to30_a10_v1/audit_t30.json)
passes **17/17** checks, including exact fixed masses and zero weights,
positive state, EOS and polygon-area consistency, no injected energy,
energy and orbit gates, and no wall crossing. Its endpoint SHA-256 is
`12ea9ff96ec70c4816f834618229c972d0a8e47814bcdc2a38a20a8e89aad50e`.

| Disk mass-weighted diagnostic | t=20 | t=30 |
|---|---:|---:|
| Mean absolute density change | `2.28372e-4` | `3.66723e-4` |
| Native-minus-exact-orbit Voronoi density residual | `3.46573e-5` | `2.10251e-4` |
| Mean exact circular-orbit position error | `1.35457e-6` | `4.17803e-6` |
| Mean absolute entropy-log change | `2.08748e-5` | `3.14251e-5` |
| Exterior mean absolute density change | `0.38155` | `0.48564` |
| Relative total-energy change from t=0 | `1.61036e-8` | `2.39930e-8` |

The disk density residual grew **6.07×** during `20→30` although the
global density gate still passes. The
[ID-matched radial receipt](KEPLER_Q2_ALIGNED_RADIAL_T20_T30_2026-10-27.json)
shows that the original `0.50≤r₀<0.65` annulus increased in mean density
change from `6.10027e-4` to `1.67480e-3` and in mean orbit error from
`7.03938e-6` to `2.42402e-5`. The outer disk annulus `1.75≤r₀<2`
remains much quieter in density (`1.06262e-4→1.86527e-4`).

The independent frozen-state pressure-force reconstructions at
[t=20](KEPLER_Q2_ALIGNED_PRESSURE_FORCE_T20_2026-10-27.json) and
[t=30](KEPLER_Q2_ALIGNED_PRESSURE_FORCE_T30_2026-10-27.json) find disk
mass-weighted pressure acceleration RMS `2.9100e-8→1.2030e-7`
(**4.13×**). The `0.50≤r₀<0.75` band rises
`7.7870e-8→3.7285e-7` (**4.79×**); the t=30 band force is primarily
tangential (`3.5229e-7` RMS). Independently rebuilt Voronoi areas match
saved native volumes within `5.38e-8` relative at t=30. The native final
RK4 ledger has `n_pair=138`, up from `20` at t=20, while floor, dual-energy
and entropy-floor injection counters remain zero. The pressure-force
readouts are snapshots, not an integrated force-history attribution.

This aligned IC remains acceptable under the present finite-time global
gates at `t=30`, but the inner disk's accelerating density, orbit and
tangential pressure-force errors do not support an asymptotically stable
claim. The distinct nonaligned 5% jitter IC already failed its stricter
native-minus-ballistic density gate at `t=10`. The active fresh-weighted
from-`t=0` and half-timestep tests should be assessed before allocating
another long aligned-only continuation. No result here establishes the
GIZMO benchmark's longest time.
