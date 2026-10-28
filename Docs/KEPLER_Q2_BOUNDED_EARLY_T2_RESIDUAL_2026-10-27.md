# Early native-minus-ballistic error grows approximately as t²

The accepted [aligned](KEPLER_Q2_BOUNDED_RADIAL_TAILS_T330024_2026-10-27.json)
and [independent 5% site-jitter](KEPLER_Q2_BOUNDED_JITTER005_RADIAL_TAILS_T330024_2026-10-27.json)
bounded q2 radial-tail receipts give the mean absolute native-minus-exact-orbit
Voronoi density residual in the initial `0.5≤r<0.75` band. The exact-orbit
reference includes each IC's own initial site positions and fixed masses;
therefore the residual excludes the density variation caused solely by
differential rotation of those sites.

| Accepted time | Aligned residual | Aligned residual / t² | Jitter residual | Jitter residual / t² |
|---:|---:|---:|---:|---:|
| 2.475808 | 2.04663e-6 | 3.33892e-7 | 6.94625e-5 | 1.13323e-5 |
| 3.300241 | 3.37538e-6 | 3.09907e-7 | 1.34349e-4 | 1.23351e-5 |
| 4.126612 | 4.78350e-6 | 2.80905e-7 | 2.11438e-4 | 1.24164e-5 |
| 4.951043 | 6.76070e-6 | 2.75802e-7 | 3.08653e-4 | 1.25915e-5 |
| 5.775473 | 8.83845e-6 | 2.64972e-7 | 4.34672e-4 | 1.30313e-5 |
| 6.500000 | 1.15159e-5 | 2.72565e-7 | 5.73423e-4 | 1.35721e-5 |

The later three jitter points have a through-origin `c t²` coefficient
`c=1.25110e-5` and ordinary three-point `R²=0.99933`. Their inner-band
RMS **position** departure from exact circular orbits also follows `t²`
(`c=2.11464e-7`, `R²=0.99963`). The aligned position departure has
`c=1.85877e-8`, `R²=0.99914`. The aligned density fit from `t=3.30–4.95`
predicts `1.18839e-5` at `t=6.5`; the observed value is `1.15159e-5`.
Before the jitter endpoint completed, the same short-time trend predicted
an inner-band `t=6.5` residual of about `5.29e-4`. This was a **near-term
discriminant**, not a long-time forecast.

The independent [matched `t=5.775473` native pair](KEPLER_Q2_BOUNDED_ALIGNED_JITTER_PAIR_T577547_2026-10-27.json)
subsequently confirmed an inner-band jitter residual of `4.34672e-4`.
Its `residual/t²=1.30313e-5` is about 4% above the earlier three-point
coefficient, so the approximate short-time law still holds but is not
an exact scaling relation. Its position-error RMS divided by `t²` is
`2.19357e-7`, similarly close to the earlier `2.11464e-7` coefficient.
The completed [matched `t=6.5` pair](KEPLER_Q2_BOUNDED_ALIGNED_JITTER_PAIR_T65_2026-10-27.json)
has inner jitter residual `5.73423e-4`, 8.5% above the early three-point
prediction `5.28592e-4`. Position-error RMS divided by `t²` rises to
`2.26435e-7`. Thus the short-time near-quadratic trend persists through
the historical failure interval, while its coefficient is drifting upward;
it cannot be projected to hundreds of orbits.

This pattern is consistent with a persistent small force or trajectory
error accumulating into a position error proportional to `t²`, with a much
larger coefficient under site disorder. It does **not** identify which
pressure/mesh/gravity term is responsible: a time-dependent shearing error
could also imitate this short trend. The jitter's inner-band signed
log-entropy at `t=4.951` is only `+1.709e-7`, so these early density and
position residuals do not themselves establish mean irreversible heating.
The matched `t=6.5` endpoint now confirms continued growth; the
[sign-resolved entropy audit](KEPLER_Q2_BOUNDED_JITTER005_ENTROPY_SIGN_T65_2026-10-27.json)
still shows nearly balanced positive and negative components. Do not
extrapolate this fit to `t=600`.

An independent [first-step velocity-defect audit](KEPLER_Q2_BOUNDED_INITIAL_VELOCITY_DEFECT_2026-10-27.json)
compares the accepted `t=.001` native velocity with each site's exactly
rotated `t=0` circular velocity. In the inner band, its mass-weighted RMS
velocity defect divided by `.001` is `6.84e-10` aligned and `1.86e-9`
jittered. Treating those first-step values as **constant** acceleration
errors would produce only `1.14e-8` and `3.11e-8` RMS position offsets by
`t=5.775`; the actual offsets are `6.35e-7` and `7.32e-6`. This simple
constant-initial-bias model underpredicts the later errors by 56 and 235
times. It is not a force decomposition, but it argues against attributing
the whole `t²` pattern to the finite-cell IC mismatch present at `t=0`.
An evolving mesh/force error or dynamical amplification remains to be
measured directly.
