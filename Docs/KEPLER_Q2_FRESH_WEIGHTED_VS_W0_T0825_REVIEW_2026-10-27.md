# Fresh Laguerre weights suppress the earliest jitter IC error

The independently audited [fresh weighted t=0.825 receipt](KEPLER_WEIGHTED_FRESH_CAPACITY_FAILURE_2026-10-27/recovery_t0825_audit.json)
and accepted [W=0 t=0.825 receipt](KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T0825_2026-10-27.json)
share the exact same independent-site-jitter `t=0` physical snapshot SHA-256
`8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b`.
Both are native CUDA evolutions to exact `t=0.825`. The weighted run
preserves masses exactly and passes 21/21 audited checks, including an
independent power-area reconstruction, EOS, positive state, no energy
injection, and pressure-force activation. Its repaired 64-corner binary
SHA-256 is `feae0a48f737a7eab215f2e31be1d7938f45f894376dcdfda177f7e09f57d26b`.

| Disk metric at t=0.825 | W=0 Voronoi | Evolving Laguerre weights | W=0 / weighted |
|---|---:|---:|---:|
| Mean absolute density change | `1.619961e-3` | `1.203383e-5` | `134.6` |
| Mean circular-orbit position error | `3.427501e-8` | `2.209657e-10` | `155.1` |
| Mean absolute entropy-log change | `4.170007e-6` | `9.163168e-7` | `4.55` |
| Disk pressure standard deviation / mean | `3.464745e-3` | `4.181074e-5` | `82.9` |

The [same-scale native density fields](KEPLER_Q2_FRESH_WEIGHTED_VS_W0_DENSITY_T0825_2026-10-27.png)
show surface density above and `Σ/Σ₀−1` below. The weighted panel uses
certified power-cell winners for rasterization; both panels use the same
display points and color scales. Tabulated measures come from original
cells, not image pixels.

The last row was recalculated directly from the two hash-verified native
endpoint snapshots using the initial-radius disk mask. The W=0
[exact-orbit Voronoi reference](KEPLER_Q2_BOUNDED_JITTER005_BALLISTIC_T0825_T65_2026-10-27.json)
already predicts `1.62104e-3` density change at this time; the native W=0
minus exact-orbit residual is only `2.15725e-6`. This pins the original
early density contrast mainly to fixed-mass Voronoi cell-area changes under
shear. The fresh weighted calculation reduces that contrast while leaving
site steering off and masses fixed.

This is strong **early-time evidence** for using Laguerre area control to
break the geometry→pressure→orbit feedback identified by the
[pressure-onset series](KEPLER_Q2_INITIAL_PRESSURE_ONSET_REVIEW_2026-10-27.md).
It is not a t>6 or t=600 result. The weighted run's stage-global projected
weight work is small, but a local moving-face mass, momentum and energy
conservation law remains to be verified. An earlier
[local-work and face-slip audit](KEPLER_Q2_NATIVE_WEIGHTED_T6501_2026-10-27/README.md)
of a late activated weighted step found nonzero cell-scale work hidden by
the global sum and only about 10% reduction in face-normal slip. That audit
used a different start state and is a caution about interpretation, not a
measurement of the fresh `t=0.825` branch. A new
[fresh-branch one-step local-work estimate](KEPLER_Q2_WEIGHTED_FRESH_LOCALWORK_T0826_REVIEW_2026-10-27.md)
likewise finds nonzero cell-scale work despite nearly zero projected global
work; it is a midpoint estimate, not a stage-exact conservation verdict.
The active fresh weighted t=0.825→6.5
continuation and subsequent exact-time comparison are required before
claiming long-time repair.
