# Bounded Kepler radial localization at accepted t=1.651374

The A10 fixed-polygon continuation, Slurm 437317, produced a complete native
checkpoint at `t=1.651374332206546` (`khout.000427.dat`, SHA-256
`b3993e3569307fceab64a5c53cb3a46fe69edac6abaf979fd73010a5337d649f`).
The [read-only radial-series auditor](audit_kepler_bounded_radial_series_20261027.py)
compared it and the accepted `t=.825` checkpoint with the same initial native
snapshot. The [machine-readable receipt](KEPLER_BOUNDED_RADIAL_T0825_T165137_2026-10-27_v2.json)
records all initial-generator radial bins, complete file hashes, whole-domain
energy and angular momentum, and generator-radius crossings. Both endpoints
retain exactly fixed mass, W=0, positive state, and polygon area tiling.

| Quantity, initial-generator bin | t=.825 | t=1.651374 |
|---|---:|---:|
| Disk r=.5–2, mass-weighted absolute entropy log | 2.52e-6 | 4.08e-6 |
| Disk r=.5–2, rms radial/circular speed | 6.76e-9 | 1.20e-8 |
| Disk r=.5–2, relative angular momentum | -1.88e-13 | -4.95e-13 |
| Inner disk r=.5–.75, absolute entropy log | 5.87e-6 | 1.05e-5 |
| Outer r=2.35–2.5, absolute density change | 8.48% | 13.98% |
| Outer r=2.35–2.5, absolute pressure change | 14.28% | 23.99% |
| Outer r=2.35–2.5, absolute entropy log | 3.76e-6 | 4.71e-6 |
| Whole-domain relative total-energy change | -6.67e-10 | -1.34e-9 |
| Whole-domain angular momentum change | -1.16e-10 | -2.33e-10 |

No generators crossed r=.5, 2, or 2.35 in either snapshot. The much larger
outer density and pressure changes coexist with very small outer entropy
change; they are primarily reversible compression/geometry at these times,
not evidence by themselves of shock heating. The disk orbit and entropy
metrics remain small through this checkpoint. The angular momentum change
of the whole domain is a useful constraint on wall torque plus numerical
error, but cannot independently separate those two contributions. Radial
bins follow **initial generators**, not face-integrated fluxes. These early
results do not establish what reaches the disk near or after `t=6`.

The next required gate is the already-running job's `t=6.5` endpoint. At its
saved intermediate checkpoints, the same auditor can locate the first radial
bin in which entropy, pressure and radial motion begin to grow, without
restarting the calculation. A separate wall-face impulse/torque ledger or
boundary-control run is still needed before attributing late disk behavior
solely to interior shear.

## Later accepted checkpoint at t=2.475808

The same A10 job subsequently saved a complete checkpoint at
`t=2.4758079598035847`. The expanded
[radial-series receipt](KEPLER_BOUNDED_RADIAL_T0825_T165137_T247581_2026-10-27.json)
finds disk mass-weighted absolute entropy log `5.04e-6`, disk rms
radial/circular speed `1.49e-8`, and inner-disk `r=.5–.75` entropy log
`1.48e-5`. In the outer `r=2.35–2.5` bin, absolute pressure change grew
to `29.20%` while entropy log remains `6.45e-6`. Whole-domain energy
relative change is `-2.01e-9`. The spatial contrast persists; no causal
claim about the wall reaching the disk follows from these snapshots alone.

## Later accepted checkpoint at t=3.300241

The [four-time radial receipt](KEPLER_BOUNDED_RADIAL_T0825_T330024_2026-10-27.json)
adds the complete native checkpoint `t=3.300240722632501`. Disk
mass-weighted absolute entropy log is `5.93e-6`, radial/circular speed
RMS is `1.64e-8`, and the inner disk `r=.5–.75` entropy log is
`1.91e-5`. The outer `r=2.35–2.5` pressure change is `31.95%`, but its
entropy log remains `7.57e-6`. No generator crosses `r=.5`, `2`, or
`2.35`. The bounded disk therefore remains accurate by these early
measures; the strong outer pressure evolution still needs localization
at later times.

## Later accepted checkpoint at t=4.126612

The running A10 job saved `khout.001703.dat` at
`t=4.126612443667704`. The read-only
[radial receipt](KEPLER_BOUNDED_RADIAL_T412661_2026-10-27.json)
confirms fixed mass, W=0, positive state and polygon area tiling. Relative
whole-domain energy change is `-3.35e-9`, and relative angular-momentum
change is `-3.82e-11`. No generator crossed `r=.5`, `2`, or `2.35`.
Disk `r=.5–2` mass-weighted absolute entropy log is `6.53e-6` and
radial/circular speed RMS is `1.74e-8`; inner disk `r=.5–.75` entropy log
is `2.33e-5`. Outer `r=2.35–2.5` absolute pressure change reached
`32.79%`, while entropy log is `7.98e-6`. These are still intermediate
snapshots before the required `t>6` gate; the outer change cannot yet be
assigned a causal role in disk heating.

## Later accepted checkpoint at t=4.951043

The [read-only radial receipt](KEPLER_BOUNDED_RADIAL_T495104_2026-10-27.json)
adds the complete `khout.002128.dat` checkpoint at
`t=4.951043426752876`. Fixed mass, W=0, positive state and area tiling
persist; relative whole-domain energy and angular-momentum changes are
`-4.01e-9` and `-4.59e-11`. Disk `r=.5–2` absolute entropy log is
`6.80e-6` and radial/circular speed RMS is `1.74e-8`. The inner disk
`r=.5–.75` absolute entropy log is `2.76e-5`. Outer `r=2.35–2.5`
absolute pressure change is `31.71%` while entropy log is `9.89e-6`.
No generator crosses `r=.5`, `2`, or `2.35`. This is still before the
historical `t>6` gate and is not a long-time benchmark pass.

## Last saved checkpoint before the historical wall, t=5.775473

The running original near32 binary saved a complete accepted
`khout.002553.dat` at `t=5.775473491592131`. Its
[read-only radial receipt](KEPLER_BOUNDED_RADIAL_T577547_2026-10-27.json)
retains exact fixed mass, W=0, positive state and polygon area tiling.
Relative whole-domain energy and angular-momentum changes are
`-4.68e-9` and `-5.35e-11`; no generator crosses `r=.5`, `2`, or `2.35`.
Disk `r=.5–2` absolute entropy log is `6.95e-6` and radial/circular speed
RMS is `1.73e-8`, essentially unchanged from `t=4.951`. Inner disk
`r=.5–.75` absolute entropy log is `3.19e-5` with signed change
`-2.26e-5`. The outer `r=2.35–2.5` pressure-change magnitude has fallen
to `29.08%` while its entropy log is `1.12e-5` and radial/circular speed
RMS has risen to `1.33e-3`. The disk remains quiet by these diagnostics,
but only the pending exact `t=6.5` endpoint can establish passage through
the previous `t>6` failure interval.
