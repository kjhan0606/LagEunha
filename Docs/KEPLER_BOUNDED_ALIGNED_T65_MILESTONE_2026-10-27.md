# Bounded aligned q2 passes the historical t>6 gate

The original accepted near32 CUDA+MPI+OpenMP binary completed A10 Slurm
job **437317** (`COMPLETED`, exit `0:0`, elapsed `01:20:48`) at **exact
`t=6.5`**. Its [native receipt](KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json)
hashes the complete 55,351-site checkpoint and passes **17 of 17** gates,
including finite positive state, fixed mass, `W=0`, polygon area, no
generator crossing of the wall, zero energy injection, and orbit and
conservation checks. This exceeds the historical `t≈6` breakdown time in
the **aligned, closed-wall control**. The original binary's dependent `t=10`
continuation is Slurm **437328**.

For the 49,038 initial disk sites (`0.5≤r<2`), the mass-weighted absolute
density change is `1.86493e-4`. Exact circular-orbit Voronoi geometry
predicts `1.85727e-4`; the cellwise native-minus-ballistic absolute mean is
`2.88254e-6`. The absolute mean logarithmic entropy change is `7.47903e-6`,
and the signed mean is `−4.93911e-6`. RMS radial speed divided by circular
speed is `1.75195e-8`; whole-domain total-energy relative change is
`+5.26279e-9`. The [radial-tail receipt](KEPLER_Q2_BOUNDED_RADIAL_TAILS_T65_2026-10-27.json)
places the largest inner-band (`0.5≤r<0.75`) native-minus-ballistic mean at
`1.15159e-5`, with 99th percentile `1.39557e-4`.

A separate [sign-resolved entropy audit](KEPLER_Q2_BOUNDED_ENTROPY_SIGN_T65_2026-10-27.json)
reproduces the native entropy absolute mean, then splits it into positive
`1.26996e-6` and negative `6.20907e-6` mass-weighted components across the
disk. About 28.5% of disk mass has positive log-entropy change. In the
inner `0.5≤r<0.75` band, the positive and negative components are
`5.29399e-6` and `3.06321e-5`, respectively. Thus the aligned run has
local positive deviations but no measured **positive mean** entropy growth
at `t=6.5`. These signed errors alone are not a face-work or irreversible
heating ledger.

The [wall-margin series](KEPLER_Q2_BOUNDED_WALL_MARGIN_SERIES_TO65_2026-10-27.json)
shows the closest generator still `0.101359` from the polygon at `t=6.5`.
Its outward wall-normal speed is `0.00242144`; this demonstrates no
generator-wall collision so far, **not** absence of wall pressure work or
future boundary influence. No wall impulse/torque time integral has been
verified. The closed polygon, intentionally aligned polar IC, and finite
cell-average initialization prevent this result from establishing the
generality or `t=600` stability of the GIZMO Kepler benchmark.

The independent 5% site-angle-jitter control **437327** was still running
when this milestone was written. The same-time `t=6.5` aligned/jitter
comparison and its radial entropy distribution are required before
attributing any failure or survival to reduced shear heating. Its queued
postprocessor **437345** depends on both endpoint gates.
