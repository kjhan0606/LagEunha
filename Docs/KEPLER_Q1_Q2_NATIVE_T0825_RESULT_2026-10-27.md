# Native polar q1/q2 Kepler at exact t=0.825 (2026-10-27)

## Acceptance and provenance

The full q2 halo-12 recovery **436995** finished `COMPLETED`, exit `0:0`,
at the exact accepted time `t=.825`; its endpoint is
`/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/halo12_after_failure_to0825/khout.010856.dat`,
SHA-256 `a8b3a3dfde83f13b9621c945d5cb52816f1cd48e51ee02be933b523f20f`.
The q1/q2 [strict comparison](/gpfs/kjhan/LagForce/Kepler_Q2_BallisticResidual_T0825_20261027/q1_q2_comparison.json)
uses independently audited complete 147,456-cell snapshots at the same time.
Both pass all 15 native gates, including exact initial mass and W=0, positive
state, periodic area sum near 36, `rho=m/V`, zero cumulative floor/dual/entropy
injection, and energy/orbit/density tolerances. The q2 [strict receipt](/gpfs/kjhan/LagForce/Kepler_Q2_BallisticResidual_T0825_20261027/native_strict.json),
[radial receipt](/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/halo12_after_failure_to0825/audit_radial_t0825.json),
and [source/input hashes](/gpfs/kjhan/LagForce/Kepler_Q2_BallisticResidual_T0825_20261027/inputs.sha256)
record the checks. Independent CPU Slurm analysis **1223641** finished and
produced the [q2 native/ballistic receipt](KEPLER_Q2_NATIVE_BALLISTIC_T0825_2026-10-27.json)
and [q2 density field](/gpfs/kjhan/LagForce/Kepler_Q2_BallisticResidual_T0825_20261027/q2_native_density.png).
The receipt copied into Docs is byte-identical to that CPU job's result; no
second interpretation or solver run is represented by the copy.

## Same-time disk measures

For cells initially in `0.5<=r<2`, using initial-mass-weighted means:

| Measure | q1 | q2 | Interpretation |
|---|---:|---:|---|
| Native absolute density change | `5.23520e-4` | `1.16521e-4` | q2 4.49 times smaller |
| Exact-orbit geometry density change | `5.23767e-4` | `1.16545e-4` | essentially the same q1/q2 ratio |
| Native versus exact-orbit cellwise density residual | `4.80819e-7` | `8.76967e-8` | small in both |
| Absolute log-entropy change | `1.53318e-4` | `1.29736e-4` | q2 15.4% smaller |
| Signed log-entropy change | `-8.45393e-6` | `-1.89861e-5` | neither has positive mean entropy change |
| Position error from exact orbit | `1.11215e-8` | `2.50253e-9` | q2 4.44 times smaller |
| Disk angular-momentum relative drift | `-9.31e-13` | `-1.86e-13` | tiny in both |

The independent [q2 periodic-image certificate](KEPLER_Q2_T0825_PERIODIC_IMAGE_SAFETY_2026-10-27.json)
finds minimum disk `edge_distance/(2*circumradius)=41.0794` at exact
`t=.825`; the nonperiodic exact-orbit polygons used for this central-disk
comparison are unaffected by omitted periodic images. All native-disk
density changes at this early time are closely predicted by exact-orbit
Voronoi geometry. This identifies q2's **aligned IC geometry** as the main
reason for its smaller early density error; it does not identify a late
entropy source or show that the general solver improved. The native W0 force
is the variational pressure operator, not an HLLC pressure flux.

The [q1 density field](/gpfs/kjhan/LagForce/Kepler_Q2_BallisticResidual_T0825_20261027/q1_native_density.png)
and q2 density field above use the same plotted density and relative-change
color limits. The q2 early result passes its acceptance gate, but the
historical failure interval `t>6`, outer periodic seam, inner core, 5%
site-jitter control, and GIZMO `t=600`/longest-orbit target remain open.
The existing q2 continuation **437076** to `t=2.221` remains pending H200
quota, with **437077** to `t=6.5` dependent on its valid completion.
