# Matched bounded Kepler IC sensitivity at t=4.126612

The aligned q2 and independent 5% site-angle-jitter q2 native CUDA
calculations reached matching accepted times `4.126612443667704` and
`4.126612443668104`. Both use 55,351 fixed-mass W=0 generators, the same
1,024-sided closed polygon, pressure `1e-6`, softened central gravity,
no steering, no AV and infinite Reynolds number. Jitter changes polar
site angles and one-time cell areas/masses, while retaining the target
continuum density and circular-speed profile. Both
[aligned](KEPLER_Q2_BOUNDED_NATIVE_T412661_2026-10-27.json) and
[jittered](KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T412661_2026-10-27.json)
17-check native receipts pass positivity, fixed mass, W=0, polygon area,
zero injection, energy, orbit, and no-wall-crossing gates.

| Disk r=.5–2, mass-weighted diagnostic | Aligned | Jitter 5% |
|---|---:|---:|
| Absolute density change | 1.606e-4 | 1.532e-3 |
| Exact circular-site Voronoi prediction | 1.605e-4 | 1.493e-3 |
| Native minus ballistic cellwise absolute density residual | 1.472e-6 | 8.408e-5 |
| Absolute log-entropy change | 6.525e-6 | 1.099e-5 |
| RMS radial/circular speed | 1.736e-8 | 8.672e-7 |

The approximately 9.5-fold density contrast is largely present in the
exact-orbit geometry prediction. The **native-minus-ballistic residual** is
about 57 times larger with jitter, however, and is growing in the inner
disk. The independent [aligned radial-tail receipt](KEPLER_Q2_BOUNDED_RADIAL_TAILS_T412661_2026-10-27.json)
and [jitter radial-tail receipt](KEPLER_Q2_BOUNDED_JITTER005_RADIAL_TAILS_T412661_2026-10-27.json)
put the `0.5<=r<0.75` mean absolute residual at `4.78e-6` versus
`2.11e-4`; jitter's 99th percentile is `7.61e-4` and maximum
`1.26e-3`. At `t=3.30024` the corresponding jitter mean was
`1.34e-4`, so this nonballistic error has increased. The inner-band
mass-weighted **signed** log-entropy change in the jittered calculation
is only `+1.05e-7`; the aligned value is `-1.63e-5`. Thus these data
show increasing orbit/partition departure under site disorder, but do
not yet establish positive mean shear heating or identify a particular
force term.

The outer `2.35<=r<2.5` absolute pressure change is `0.32789564` in
both calculations to the displayed precision, and neither has a generator
crossing `r=.5`, `2`, or `2.35`. This localizes the **difference between
these two ICs** primarily to the perturbed disk at this time; it does not
prove that the fixed wall has no influence at later times.

The running aligned **437317** and jittered **437327** calculations must
both pass their `t=6.5` native endpoint gates before interpreting survival
past the historical `t>6` failure. If the residual keeps growing without
matching entropy production, compare pressure/mesh force terms and orbit
errors before introducing a spin or viscous-stress operator. A successful
aligned run alone remains a special polar-layout result, not a general
Kepler benchmark pass.

The [jitter t=10 runner](kepler_q2_bounded_jitter005_to10_a10_20261027.slurm)
and [guarded prestager](stage_kepler_bounded_jitter005_to10_20261027.py)
are queued as Slurm **437344** with verified `afterok:437327`. They copy
the jittered `t=6.5` checkpoint only after its complete 17-check native
receipt passes, retain the exact original near32 binary, and audit the
`t=10` endpoint. The aligned control has its separate dependent `t=10`
job **437328**. Neither `t=10` result exists yet.

The [paired endpoint auditor](audit_kepler_bounded_aligned_jitter_pair_20261027.py)
was smoke-tested on these `t=4.126612` accepted receipts; it reproduced
the 57.14-fold whole-disk nonballistic residual ratio and 44.20-fold inner
band ratio. A CPU-only [postprocessor](kepler_bounded_aligned_jitter_post_t65_cpu_20261027.slurm)
is queued as Slurm **437345**, with both `afterok:437317` and
`afterok:437327` verified in the scheduler. It requests two CPU threads,
16 GB, no GPU, and produces same-time radial tails, a paired JSON report,
and a common-scale density-field image in the isolated
`/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/aligned_jitter_post_t65_cpu_v2`
directory. The first unsubmitted staging directory `v1` is preserved:
its runner specified a `normal` partition absent from the `syntax` Slurm
controller. The corrected `v2` uses the available `a10` partition without
requesting a GPU. No `t=6.5` comparison result exists yet.

## Matched continuation at t=4.951043

The next equal-time [aligned](KEPLER_Q2_BOUNDED_NATIVE_T495104_2026-10-27.json)
and [jittered](KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T495104_2026-10-27.json)
native checkpoints both pass all 17 physics checks. Whole-disk mean
absolute density changes are `1.76155e-4` and `1.56104e-3`, while exact
circular-site geometry predicts `1.75887e-4` and `1.47901e-3`.
The cellwise native-minus-ballistic density residual has risen to
`1.93825e-6` aligned and `1.46732e-4` jittered: the latter is now
75.7 times the aligned residual and 1.75 times its own `t=4.1266`
value. The jittered disk radial/circular speed RMS is `1.031e-6`, still
small in absolute terms but 59 times the aligned value.

The [aligned radial tails](KEPLER_Q2_BOUNDED_RADIAL_TAILS_T495104_2026-10-27.json)
and [jittered radial tails](KEPLER_Q2_BOUNDED_JITTER005_RADIAL_TAILS_T495104_2026-10-27.json)
place the largest difference in `0.5<=r<0.75`: mean absolute
nonballistic residual `6.76e-6` versus `3.09e-4`, with a jittered 99th
percentile `1.12e-3` and maximum `1.69e-3`. The jittered band's signed
mass-weighted log-entropy change remains only `+1.71e-7` (aligned
`-1.95e-5`). Thus the disorder-sensitive departure is growing before
any demonstrated positive **mean** entropy production. This does not
identify which part of the pressure/mesh coupling causes it, or predict
whether it becomes irreversible after `t=6`.

## Matched native result at t=6.5

Both the unchanged-binary aligned job **437317** and 5% independent
site-jitter job **437327** reached exact `t=6.5`, exited `0:0`, and passed
all 17 native fixed-mass, W=0, positivity, geometry, conservation and
orbit checks. The latter's [copied, snapshot-rehashed receipt](KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json)
has endpoint SHA-256
`d02e15ebd01193f1da002d87d689f83525b0150541a9f78dedbc451702c5a2ff`.
The dependent original-binary jitter `t=10` job **437344** started from
this accepted state. This is a continuation, not a restart from the IC.

The [automated paired receipt](KEPLER_Q2_BOUNDED_ALIGNED_JITTER_COMPARISON_T65_2026-10-27.json)
and common-scale [density-field image](KEPLER_Q2_BOUNDED_ALIGNED_JITTER_DENSITY_T65_2026-10-27.png)
show disk mass-weighted mean absolute density drift `1.86493e-4` aligned
versus `1.69828e-3` jittered (9.11-fold). Exact circular-site Voronoi
geometry alone contributes `1.85727e-4` versus `1.47879e-3` (7.96-fold).
The **cellwise native-minus-ballistic** density residual is `2.88254e-6`
versus `3.44253e-4`, or **119.43-fold**. Thus initial geometric disorder
explains much of the absolute density contrast, but the evolved jittered
state also departs much more strongly from its own exact-orbit geometry.
The two absolute density images still look similar at the common broad
color scale; their relative-density panels expose the different fine
structure.

The [jitter radial-tail receipt](KEPLER_Q2_BOUNDED_JITTER005_RADIAL_TAILS_T65_2026-10-27.json)
places the largest nonballistic residual in `0.5<=r<0.75`: `5.73423e-4`
mean absolute, versus `1.15159e-5` aligned, with jittered 99th percentile
`2.08640e-3`. The jittered disk signed mass-weighted log-entropy change
is `-3.37355e-7` and the inner-band value `+1.47936e-7`; both are near
zero. The whole-disk energy errors are essentially equal at about
`5.25e-9`. These facts identify a growing disorder-sensitive **local**
orbit/partition error before demonstrated positive mean heating. They do
not distinguish pressure-force error from trajectory-induced Voronoi
sensitivity or establish a need for spin. Both `t=10` controls and a
face-force/trajectory diagnostic are needed before changing the operator.
