# Periodic Q1 seam and close-pair postmortem

This note records a read-only analysis of the earlier 384² periodic Q1 trajectory. It supplements the same-time comparison in [RESULT_T3.md](RESULT_T3.md); it does not modify or supersede the bounded-domain t=3 or t=10 results.

## Saved states and method

The compared full native snapshots are:

- t=0: `/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/khout.000000.dat`, SHA-256 `edb1f9d30700a95093a4f442afd3461238dd4e93ad9b075302557c15dcd693c9`.
- t=3: `/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/halo40_from_t2654918_a100_20261028/to3p0/khout.027642.dat`, SHA-256 `cfa09e4f0f8b074bed39b35fe15d2e48fb14d1f3bb7c86a3968aaf2539c02bae`.

Both contain 147,456 sites in a 6×6 periodic box centered at `(3,3)`. Using the native snapshot reader and `scipy.spatial.cKDTree(..., boxsize=6)`, we measured each site's nearest periodic generator distance `d_nn` and normalized it by `sqrt(V)`.

| Metric | t=0 | t=3 |
| --- | ---: | ---: |
| Minimum `d_nn` | `1.495723e-2` | `1.180661e-4` |
| Minimum `d_nn/sqrt(V)` | `0.7213` | `0.01049` |
| Sites with `d_nn/sqrt(V) < 0.05` | 0 | 125 |

At t=3, 91 of the 125 close sites are at `2 <= R < 2.75`; 10 are inside `R < 2`; the remaining 24 are at `R >= 2.75`. The closest pair is near `R=2.25`, with separation `1.18e-4`, relative speed `1.19`, and densities `0.0635` and `0.0214`. A group near `R=2.55` has separations around `1.39e-4` and relative speeds around `0.99`.

A separate pair crosses the periodic x seam: positions approximately `(5.99064, 3.20672)` and `(0.00936, 3.20657)` are `2.30e-4` apart under the periodic metric, at `R≈3.0`, with relative speed `1.15`. This is consistent with a large velocity difference across the seam for circular motion about the box center.

## Interpretation and limits

The t=0 to t=3 comparison demonstrates that very close generator pairs form during the periodic evolution. Combined with the prior active-source limiter replay, which found near-coincident exterior pairs setting the face-CFL minimum around t≈1.5, this is a plausible explanation for the tiny global steps and a performance cost concentrated outside the scientific disk. These saved-state distances alone do **not** prove which face set the live minimum at t=3.0047. The exact per-call limiter replay remains the causal check.

Hopkins (2015) specifies the softened circular velocity and density tail outside the main disk, and notes that mesh methods use a `Σ>=0.01` floor and smoothed edges. Its effective mesh setup is `320²` over `[-2.5,2.5]²`; this Q1 branch instead uses `384²` sites in a 6×6 periodic box. The cited §4.2.4 passage does not state the Kepler test's boundary condition, so this comparison identifies a setup difference and a seam risk without asserting that the paper uses the same periodic boundary. See [Hopkins (2015), §4.2.4](https://arxiv.org/abs/1409.7395).

The bounded-domain continuation remains the stronger physical branch: it preserves the original disk state while replacing the periodic exterior with a finite polygon. Its disk gates through t=10 pass, while the core remains inaccurate; see [RESULT_T10.md](RESULT_T10.md). This postmortem supports retaining the periodic branch as a diagnosed comparison rather than treating it as the final GIZMO reproduction.
