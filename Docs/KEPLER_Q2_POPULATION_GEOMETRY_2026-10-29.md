# Q2 ring-population geometry control

An exploratory control changes each of the 68 original disk ring populations
alternately by +1 and -1 while preserving its radius, uniform angular spacing,
one angular anchor, the total number of disk sites, and all non-disk sites.
The modified inner band has 3,932 sites instead of 3,931. This creates an
initial placement only; there is no generator steering during evolution.

For each IC, frozen Voronoi cell areas are reconstructed at t=0 and after
prescribed softened circular motion to t=6.5. Assuming initial uniform
density and fixed individual masses, their ratio gives the kinematic density
change. This is not a hydro evolution. Initial original areas agree with the
stored native areas to 3.03e-8 relative maximum; independent radial support
cutoffs 1.2 and 1.4 agree to 1.4e-12 relative maximum for both ICs and times.

| Inner band `0.5<=r0<0.75` | Original | Alternating population |
| --- | ---: | ---: |
| Initial area peak modes on 11 rings | all 13 | 14 on one, 15 on ten |
| t=6.5 mass-weighted mean absolute density change | 4.74052e-4 | 5.12006e-4 |
| t=6.5 cell p99 absolute density change | 1.41776e-3 | 1.66913e-3 |

The population change shifts the initial geometric pattern, but raises the
kinematic mean density variation by about 8% and p99 by about 18% in this
example. The result does not support adopting this IC as a cure. It provides
a controlled hypothesis for whether late modes follow the initial population
pattern, provided a future native run uses independently recomputed masses,
full-domain wall geometry, matched physics and conservation checks. No
conclusion about pressure forces, entropy, long-time stability, or equivalence
to the GIZMO IC follows from the kinematic reference alone. The altered
placement and its masses have not been installed or run in the native solver.

The [receipt](KEPLER_Q2_POPULATION_GEOMETRY_2026-10-29.json) records all
counts, per-ring peaks and source hashes; the [reader](audit_kepler_q2_population_geometry_20261029.py)
reconstructs the experimental positions without changing the stored snapshot.
