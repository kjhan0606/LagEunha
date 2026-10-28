# Equal-time bounded q2 Kepler comparison at t=5.775473

The 5% independent site-angle-jitter native CUDA+MPI+OpenMP run saved a
complete checkpoint at `t=5.775473491594761`, within `2.63e-12` of the
aligned run's accepted `t=5.775473491592131`. The
[jitter native receipt](KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T577547_2026-10-27.json)
passes all 17 physics checks; the already accepted
[aligned receipt](KEPLER_Q2_BOUNDED_NATIVE_T577547_2026-10-27.json)
passes the same. Both have fixed mass, `W=0`, positive EOS, zero RK4
injection and no generator crossing of the fixed 1,024-gon.

The [paired audit](KEPLER_Q2_BOUNDED_ALIGNED_JITTER_PAIR_T577547_2026-10-27.json)
compares two real native checkpoints at the same physical time. For
the initial `0.5≤r<2` disk, the mass-weighted absolute density change is
`1.75064e-4` aligned versus `1.63149e-3` jittered. The exact circular-site
Voronoi prediction is `1.74687e-4` versus `1.48716e-3`; much of the
raw density contrast is geometric. The **native-minus-ballistic** density
residual is `2.40041e-6` versus `2.37918e-4`, a factor **99.1**. In the
inner `0.5≤r<0.75` band, the same residual is `8.83845e-6` versus
`4.34672e-4`, a factor **49.2**. This extension follows the previously
observed approximate `t²` inner-band growth.

The [native density-field figure](KEPLER_BOUNDED_ALIGNED_JITTER_NATIVE_DENSITY_T577547_2026-10-27.png)
uses the same logarithmic surface-density and signed relative-change scales
for both ICs. The disk remains intact in both top panels; the lower panels
resolve the larger, fine-scale jitter density pattern. The figure is from
the native gas states, while the exact-orbit comparison above is needed to
separate geometric variation from physical-state departure.

The sign-resolved [aligned](KEPLER_Q2_BOUNDED_ENTROPY_SIGN_T577547_2026-10-27.json)
and [jittered](KEPLER_Q2_BOUNDED_JITTER005_ENTROPY_SIGN_T577547_2026-10-27.json)
entropy audits use the same accepted snapshots and reproduce the native
absolute entropy metrics. In the jittered disk, the mass-weighted positive
and negative log-entropy components are `6.76329e-6` and `7.03121e-6`;
the signed mean is only `−2.67919e-7` and 49.9% of disk mass has a
positive change. The jittered inner band has a signed mean `+2.30076e-7`
with positive and negative components `3.40216e-6` and `3.17208e-6`.
Thus the much larger density and orbit departure does **not** yet
coincide with a demonstrated positive mean entropy production. These
components are state differences, not an integrated face-work ledger.

An additional adiabatic pressure comparison uses each initial pressure
and the **exact circular-site Voronoi area** to predict
`P_ballistic/P₀=(A₀/A_ballistic)^γ`. Its [aligned](KEPLER_Q2_BOUNDED_BALLISTIC_PRESSURE_T577547_2026-10-27.json)
and [jittered](KEPLER_Q2_BOUNDED_JITTER005_BALLISTIC_PRESSURE_T577547_2026-10-27.json)
receipts find disk mean absolute native relative pressure changes of
`2.91850e-4` and `2.71941e-3`. Geometry alone predicts `2.91146e-4`
and `2.47862e-3`; mass-weighted native/predicted pressure correlations
are `0.99930` and `0.98872`. In the jittered inner band, native and
predicted mean absolute pressure changes are `3.03727e-3` and
`2.55981e-3`. Thus moving fixed-mass Voronoi volume changes account for
most of the pressure pattern even on exact circular generator orbits.
The unmeasured next link is whether the native W0 pressure **force** from
that pattern drives the observed orbit residual; pressure correlation
alone does not prove force causality.

The matched `t=6.5` endpoint and dependent CPU-only postprocessor remain
the decisive next short-time gate. This one seed and one deliberately
aligned control do not establish general solver accuracy or `t=600`
stability; the noncircular closed wall remains a separate long-time
boundary uncertainty.
