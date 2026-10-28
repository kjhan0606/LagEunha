# Pressure-force onset in the three native Kepler ICs

The [frozen-state series receipt](KEPLER_Q2_INITIAL_PRESSURE_SERIES_2026-10-27/summary.json)
reconstructs the native W=0 variational face pressure force from five saved
states in each accepted `t=0→6.5` CUDA run. The same auditor previously
matched the native GPU force capture on 64 disk cells at `t=6.5` to
`≤1.01e-11` absolute acceleration and rechecks independent Voronoi areas
against saved native volumes. The source audit and snapshot hashes are in
the series receipt. A separate reconstruction of each hash-verified `t=0`
snapshot gives exactly zero pressure acceleration and pressure variation
for all three ICs: [aligned](KEPLER_Q2_INITIAL_PRESSURE_FORCE_ALIGNED_2026-10-27.json),
[ring-coherent](KEPLER_Q2_INITIAL_PRESSURE_FORCE_RING_2026-10-27.json),
and [independent jitter](KEPLER_Q2_INITIAL_PRESSURE_FORCE_JITTER_2026-10-27.json).

| Saved time, approximately | Aligned pressure acceleration RMS | Ring-coherent RMS | Independent jitter RMS | Jitter / aligned |
|---:|---:|---:|---:|---:|
| `0` | `0` | `0` | `0` | undefined |
| `0.825` | `1.4320e-8` | `1.4915e-8` | `2.2212e-7` | `15.51` |
| `2.476` | `1.6008e-8` | `1.6596e-8` | `1.9832e-7` | `12.39` |
| `4.127` | `1.5741e-8` | `1.6454e-8` | `2.1434e-7` | `13.62` |
| `5.775` | `1.6662e-8` | `1.7180e-8` | `2.3557e-7` | `14.14` |
| `6.5` | `1.7461e-8` | `1.8197e-8` | `2.5069e-7` | `14.36` |

The independent-jitter force is about 95% tangential in RMS magnitude at
all five saved times. At `t≈0.825`, its pressure standard deviation/mean
is `3.46e-3`, compared with `3.96e-4` aligned. The separately accepted
[native t=0.825 audit](KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T0825_2026-10-27.json)
measures disk mean absolute density change `1.61996e-3` and
native-minus-exact-orbit density residual only `2.15725e-6`. Its
[exact-orbit Voronoi reference](KEPLER_Q2_BOUNDED_JITTER005_BALLISTIC_T0825_T65_2026-10-27.json)
already predicts `1.62104e-3` density change at the same time. Thus the
large early pressure contrast follows a density pattern that is nearly
entirely present under prescribed circular shear and fixed-mass Voronoi
areas, before a large site-orbit error has developed.

This sequence strengthens a specific mechanism: irregular within-ring
spacing creates kinematic cell-area and density variation; the EOS turns
it into pressure contrast, and the native variational force acquires a
large tangential component. The later
[inherited-orbit control](KEPLER_Q2_JITTER_INHERITED_ORBIT_T10_REVIEW_2026-10-27.md)
shows that pressure-off evolution from the already perturbed `t=6.5`
state carries nearly all of its `t=10` residual forward. The subsequent
[same-t0 pressure-force-off control to t=0.825](KEPLER_Q2_EARLY_PRESSURE_CAUSAL_TRIAD_T0825_2026-10-27.md)
directly separates the earliest response: turning the force off leaves the
`1.62e-3` geometric density change in place but reduces mean orbit error
from `3.43e-8` to `4.09e-13`. The initial density disturbance is therefore
kinematic, while pressure force creates nearly all measured orbital
deviation by t=0.825 in this diagnostic comparison. This does **not** yet
quantify the force contribution to the t=6.5 error or establish a
physical pressure-free solution. The five intermediate force snapshots
have no separate full native physics gate. The fresh weighted-from-`t=0`
and half-timestep continuation remain the relevant live tests before
choosing another long branch.
