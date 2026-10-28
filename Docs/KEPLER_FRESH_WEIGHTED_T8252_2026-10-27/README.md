# Fresh weighted Kepler comparison at t≈8.252

The native W=0 Voronoi run and fresh weighted Laguerre run start from the same
initial snapshot (`8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b`).
Their saved endpoint times are `8.251666223756754` and
`8.2516665027481`, a difference of `2.78991e-7`. The independent comparison
receipt passes the common-initial-state, time-alignment, entropy-lineage and
native branch-gate checks. Both use native CUDA evolution and fixed individual
cell masses; the weighted branch evolves Laguerre weights.

| Disk metric | W=0 Voronoi | Evolving weights | Weighted / W=0 |
|---|---:|---:|---:|
| Mass-weighted mean absolute density change | `1.96247e-3` | `4.05826e-5` | `0.02068` |
| Mean circular-orbit position error | `9.84606e-6` | `8.00112e-8` | `0.00813` |
| Mean absolute log-entropy change | `1.77206e-5` | `8.01735e-6` | `0.45243` |
| Relative total-energy error | `6.66183e-9` | `6.65947e-9` | — |

At this matched checkpoint, the weighted run has about 48 times less density
change and 123 times less orbital-position error. The absolute entropy-change
metric is lower by about 55%; the signed entropy changes differ in sign between
the branches, so these values alone do not establish irreversible heating or
its absence. Both branch energy errors are about `6.66e-9`.

This supports the hypothesis that controlling cell areas with evolving
Laguerre weights suppresses the fixed-mass Voronoi area error under Keplerian
shear, while site steering remains off. It does **not** prove that the
weighted moving faces obey a local mass, momentum and energy conservation law.
The branches use different binaries and trajectory/restart histories; their
global energy ledgers do not certify the face-level exchange. No conclusion
about stability through `t=20`, `t=60`, the 384² periodic run, or the GIZMO
`t=600` duration follows from this checkpoint.

Across the preceding matched late checkpoint at `t≈7.4253`, the weighted
disk density-change metric was `2.62810e-5` versus `1.83783e-3` for W=0.
At `t≈8.2517` these are `4.05826e-5` and `1.96247e-3`. The weighted error
rose between these saved times, while remaining about 1.4–2.1% of the W=0
value; it is small, not identically zero. The separate W=0 `t=10` ballistic
residual is `1.19668e-3`, above its `1e-3` gate. The active fresh weighted
run's exact `t=10` audit is still needed to determine whether this weight
prescription clears that wall.

This is **not a rotation run**. Its launcher explicitly unsets
`GFS_LAGUERRE_ROTATION`. In the current CUDA/CPU force code that option only
adds a face-normal-velocity correction for an off-midpoint Laguerre face; it
does not evolve a cell spin, spin angular momentum, or rotational energy.
Consequently this comparison says nothing about whether intrinsic cell
rotation suppresses shear-driven heating.

The comparison and branch receipts are preserved here. `density.png` is the
native surface-density field and per-cell density-change rendering on common
color scales, rasterized using the actual Voronoi/power-cell winners. The
rendering records the endpoint time difference explicitly.
