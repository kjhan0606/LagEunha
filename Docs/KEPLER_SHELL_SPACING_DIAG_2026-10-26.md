# Kepler shell-spacing hypothesis and C010 failure (2026-10-26)

## What the accepted C010 snapshots show

C010 begins with Cartesian generator positions, so it has no designed circular
shells. In the first disk annulus, `2 <= r0 < 2.25` (368 cells), the initial
Cartesian nearest-neighbor distance is exactly `H=24/256=0.09375`. To test the
shell idea operationally, classify a neighbor as approximately on the same
shell when `|r0,j-r0,i| < H/4 = 0.0234375`; all others are cross-band neighbors.
The classification is an analysis convention, not a feature of the IC.

| Quantity, first disk annulus | t=0 | t=6.500079 |
|---|---:|---:|
| Median nearest same-band distance | 0.132583 | 0.144404 |
| Median nearest cross-band distance | 0.093750 | 0.064731 |
| Fraction with cross-band neighbor closer than same-band neighbor | 63.0% | 78.3% |
| All-neighbor median minimum distance | 0.093750 | 0.062971 |

At t=6.5, 29.3% of cells in this annulus have some neighbor at less than
`H/2`; the closest separation is `0.026336`. Across the 368 cells, Spearman
correlation between nearest-neighbor distance and absolute density change is
`-0.4006`. The closest 40 cells have mean absolute density change `0.4233`,
versus `0.1187` for the rest. Repeating the same/cross-band classification with
radial thresholds 0.015, 0.020, and 0.030 preserves the trend: the fraction
whose cross-band neighbor is closer rises from 73.9% to 88.0%, 67.4% to
81.5%, and 52.2% to 66.3%, respectively.

These observations support cross-radius crowding as a geometry mechanism.
They do not prove nearest-neighbor distance alone causes a volume or pressure
change. The independent exact-orbit Voronoi calculation predicts the native
first-annulus mass-weighted absolute density change at t=6.5 to `7.68e-5`
residual, so generator geometry is already established as the dominant source
of that disk signal. The accepted C010 density change there is `0.1469215`.

## Controlled geometry-only polar layouts

The [control script](kepler_ring_spacing_ballistic_control_20261026.py)
constructs staggered polar rings over `1.25<r<4.5`, chooses radial spacing
`a=H sqrt(q)` and target azimuthal spacing `b=H/sqrt(q)` for `q=a/b`, and
assigns fixed cell masses `m_i=(2/r_i) A_i(0)`. This preserves nominal local
area `ab=H^2`. Every generator then follows its exact softened circular
Kepler orbit; SciPy independently recomputes unweighted Voronoi cell areas.
The reported statistic is the mass-weighted `|A_i(0)/A_i(t)-1|` for cells with
`2<=r0<2.25`. The constant-count variant has equal particle count on all
rings, while the variable-count variant approximately keeps azimuthal
spacing constant across radii.

| Layout | q=0.5, t=6.5 | q=1, t=6.5 | q=2, t=6.5 | q=2, t=17.779 |
|---|---:|---:|---:|---:|
| Polar rings, constant count | 0.001162 | 0.000172 | 0.000001 | 0.000005 |
| Polar rings, variable count | 0.004426 | 0.001139 | 0.000952 | 0.000856 |
| Original C010 Cartesian IC, exact orbit | 0.146948 | — | — | — |

These are **ballistic geometry controls**, not native hydrodynamic runs.
Crucially, even the polar `q=0.5` layout, whose cross-ring spacing is smaller
than its along-ring spacing, has much less area noise than the Cartesian IC.
Within the more realistic variable-count polar layouts, moving from `q=1`
to `q=2` improves t=6.5 error only from `0.001139` to `0.000952` (about 16%).
Thus **spacing ordering alone is neither a sufficient explanation nor the
dominant measured difference in this control**; circular alignment and ring
regularity appear to matter greatly. The polar-ring symmetry changes along
with the spacing, and the `2<=r0<2.25` measurement contains 2–4 distinct
rings across variants. The test cannot establish stability of a gas run to
the GIZMO paper's final orbit. An affine-sheared infinite lattice can also
keep constant Voronoi areas, so shear amplitude by itself is insufficient;
radial curvature, discrete pattern, and finite transition regions matter.

The next causal test is an independent native fixed-mass, W=0 polar-ring IC
with the same physics, resolution budget, outer/inner ambient treatment,
complete-snapshot audits, and density images. Test both `q=1` and `q=2`, then
follow beyond t=10 and toward the paper's final orbit only if valid. Do not
substitute this geometry-only result for that test.

## Separate C010 whole-domain stop

Slurm job 436489 ended `FAILED` after the accepted t=6.5 checkpoint. The
native log states `RK4_INVALID_STATE` at t=6.5989249904980696, step 40724,
`k4_input`, cell ID 71058: `U=-5.3970e-10`, volume `2.8233e-6`, mass
`1.4051e-7`. Its initial radius was `0.4900`, so the failure is in the
low-density inner atmosphere, not the first disk annulus. Its radius had
moved to `0.3620`. The start-of-step nearest generator was only `4.7e-5`
away, versus `0.011719` at t=0; its volume was `1.7604e-6`, versus
`1.3733e-4` initially. At the failed RK stage, that pair was `2.9e-5` apart.
This is a severe central crowding and stage-positivity problem and must be
diagnosed separately from the disk's exact-orbit Voronoi density signal.
At the failing step start, the two almost coincident cells (IDs 71058 and
68323) are respectively 0.817 and 0.789 position units from their own exact
softened circular-orbit predictions. Their initial radii were 0.4900 and
0.4367. Thus this central collapse has substantial actual trajectory drift;
it cannot be attributed solely to changing Voronoi volume while all
generators remain on perfect orbits.
The conditional resume job 436713 correctly refused to restart a parent
whose state was `FAILED`; it was designed to resume only preemption/timeout.

Machine-readable measurements:
[neighbor data](KEPLER_NEIGHBOR_SPACING_T65_2026-10-26.json),
[polar layout data](KEPLER_RING_SPACING_BALLISTIC_CONTROL_2026-10-26.json),
and [native versus ballistic t=6.5 data](BALLISTIC_VORONOI_CONTROL_T65_2026-10-26.json).
