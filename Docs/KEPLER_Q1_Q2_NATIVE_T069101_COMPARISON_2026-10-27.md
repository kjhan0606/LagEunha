# Matched native polar Kepler states at t=0.6910082452

The accepted q1 and q2 halo-8 step-start snapshots precede their later
periodic-halo aborts and share the same physical time. Their independent
[q1](KEPLER_Q1_FAILURE_START_RADIAL_T069101_2026-10-27.json) and
[q2](KEPLER_Q2_FAILURE_START_RADIAL_T069101_2026-10-27.json) radial audits
each pass complete IDs, fixed mass, W=0, positive state, and density
consistency. Source JSON SHA-256 values are `9d498c755b26841a97a6d113cca77140ac2d2a309d3f2dce04ed72123579cce3`
and `79b11bc4d1430dcae7ca166425d878322307d587d3018293685e163f1c920554`,
respectively. The figures below combine the three disjoint initial-radius
bins `0.5–0.75`, `0.75–1`, and `1–2`, weighting each bin mean by its recorded
initial mass. They are native gas results, not exact-orbit predictions.

| Initial disk, 0.5 <= r < 2 | q1 | q2 | q2/q1 |
|---|---:|---:|---:|
| Fixed initial mass | 11.86552375 | 11.90497938 | 1.0033 |
| Mean absolute relative density change | 5.51643e-4 | 1.13661e-4 | 0.206 |
| Mean absolute relative pressure change | 9.17421e-4 | 2.31622e-4 | 0.252 |
| Mean absolute entropy-log change | 1.34969e-4 | 1.07981e-4 | 0.800 |
| Signed mean entropy-log change | -6.07079e-6 | -1.40253e-5 | — |
| Mean absolute radial speed | 2.07913e-8 | 4.28591e-9 | 0.206 |
| Mean exact-orbit position error | 7.46717e-9 | 1.73206e-9 | 0.232 |

The density benefit is broad across the disk: q2/q1 is 0.404 in the inner
`0.5–0.75` bin, 0.226 in `0.75–1`, and 0.160 in `1–2`. The entropy trend
does not scale with it. In `0.75–1`, q2's absolute entropy-log error is
`1.30393e-4`, **larger** than q1's `8.93788e-5`, while q2's density error
is much smaller. Both disk-wide signed entropy changes are negative at this
time; neither is evidence of net irreversible heating. This comparison
supports separating geometric density preservation from entropy/force error
when evaluating the q2 IC and any proposed cell-spin repair.

The shared outer Cartesian region `r>=2.75` has the same reported density
change `0.0676667` and absolute radial speed `0.0602113` for both ICs at
this time. The severe periodic-seam disturbance is therefore present even
when q2 improves the disk. The two layouts differ in disk site distribution
and initial disk mass, so these ratios measure an IC comparison, not an
isolated change to one face operator. Native q2 `t=.825`, q1/q2 `t=2.221`
and especially `t>6` with the same physical settings are needed to decide
whether the density advantage delays late degradation. The queued
site-jitter comparison tests how much of q2's advantage relies on the
regularity of its rings; its early endpoint alone cannot establish a
long-time heating result.
