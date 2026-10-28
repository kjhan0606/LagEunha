# Matched bounded q2 Kepler at the historical t>6 gate

The native CUDA+MPI+OpenMP aligned polar q2 job **437317** and independent
5% site-angle-jitter q2 job **437327** both completed exact `t=6.5` with
Slurm `COMPLETED 0:0`. Their [aligned](KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json)
and [jittered](KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json)
states each pass all **17** native physics gates: finite positive EOS,
fixed mass, `W=0`, polygon-area closure, zero RK4 floor/dual-energy
injection, energy/orbit limits, and no generator crossing of the closed
1,024-gon. The dependent CPU-only comparison **437345** also completed
`0:0`; its [paired receipt](KEPLER_Q2_BOUNDED_ALIGNED_JITTER_PAIR_T65_2026-10-27.json)
uses the two complete accepted endpoints.

| Initial `0.5≤r<2` disk metric | Aligned | 5% site jitter |
|---|---:|---:|
| Mean absolute density change | 1.86493e-4 | 1.69828e-3 |
| Exact circular-site Voronoi prediction | 1.85727e-4 | 1.47879e-3 |
| Native-minus-ballistic mean absolute density residual | 2.88254e-6 | 3.44253e-4 |
| Mean absolute log-entropy change | 7.47903e-6 | 1.49389e-5 |
| RMS radial/circular speed | 1.75195e-8 | 1.22209e-6 |
| Whole-domain relative total-energy change | +5.26279e-9 | +5.24739e-9 |

The jittered native-minus-ballistic residual is **119.4 times** the
aligned residual across the disk. In the initial `0.5≤r<0.75` band it is
`5.73423e-4` versus `1.15159e-5`, a factor **49.8**. The inner jitter
residual is 8.5% above the prior three-point `t²` extrapolation; the
short-time near-quadratic growth continues but its coefficient is rising.
The [native density image](KEPLER_BOUNDED_ALIGNED_JITTER_NATIVE_DENSITY_T65_2026-10-27.png)
uses a shared color scale and shows an intact disk in both ICs with a
much stronger fine-scale density pattern under independent site jitter.

The [jitter sign-resolved entropy audit](KEPLER_Q2_BOUNDED_JITTER005_ENTROPY_SIGN_T65_2026-10-27.json)
finds disk positive and negative mass-weighted log-entropy components
`7.30077e-6` and `7.63813e-6`; the signed mean is `−3.37355e-7`.
In the inner band the signed mean is `+1.47936e-7`, compared with
absolute mean `6.99605e-6`. These nearly balanced errors do not establish
positive **mean** irreversible heating. They also do not rule out local
dissipation because no per-face irreversible-work ledger was integrated.

The exact-orbit adiabatic pressure comparison remains strong: the
[aligned](KEPLER_Q2_BOUNDED_BALLISTIC_PRESSURE_T65_2026-10-27.json)
and [jittered](KEPLER_Q2_BOUNDED_JITTER005_BALLISTIC_PRESSURE_T65_2026-10-27.json)
disk pressure correlations between native and geometry-only predictions
are `0.99907` and `0.98078`. Jittered native mean absolute relative
pressure change is `2.83063e-3`; circular-site area evolution alone
predicts `2.46465e-3`. The unpredicted pressure component grew from
`3.97992e-4` at `t=5.775` to `5.74957e-4` at `t=6.5`.

The [jitter wall-margin audit](KEPLER_Q2_BOUNDED_JITTER005_WALL_MARGIN_T65_2026-10-27.json)
finds minimum generator-wall gap `0.101359`, effectively identical to
the [aligned audit](KEPLER_Q2_BOUNDED_WALL_MARGIN_SERIES_TO65_2026-10-27.json).
This excludes a generator-wall collision through `t=6.5`, not an
unmeasured pressure wave or wall torque.

**Interpretation and next gate.** This establishes survival beyond the
old `t≈6` failure time for **both** a specially aligned and a slightly
disordered bounded IC. It does not establish generality, an open-boundary
GIZMO-equivalent problem, or stability to the paper's `t=600` image.
The current W0 CUDA source applies `ΔP` times a face-geometry coefficient
to each cell's momentum; the geometric pressure pattern is therefore a
plausible source of the growing orbit error. A selected-face force dump
at an accepted late restart should compare that measured pressure
acceleration with the subsequent gravity-subtracted velocity change
before changing Laguerre weights or introducing cell spin. Existing
aligned **437328** and jittered **437344** continuations are gated to
exact `t=10` using their own accepted original-binary checkpoints.

The subsequent [original-binary selected-face probe](KEPLER_Q2_BOUNDED_FACEFORCE_CAUSAL_T6501_2026-10-27.md)
completed this one-step force check: on 64 inner-disk cells the reconstructed
W0 pressure acceleration matches the observed gravity-subtracted velocity
change within 1.8% relative vector L2 in the aligned case and 0.8% in the
jittered case. Exact-orbit cell-area pressure explains much, but not all,
of the jittered pressure force. The result identifies a concrete operator
target while leaving long-time and boundary validation open.
