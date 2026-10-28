# q1 native disk versus exact-orbit geometry at t=0.825

An independent CPU Slurm job **1223481** (completed in 4 s, exit `0:0`)
applied the existing [early polar checkpoint auditor](audit_kepler_polar_checkpoint_20261026.py)
to the completed native q1 halo-12 `khout.010856.dat`. The staged script,
reader layout, initial and endpoint snapshots are SHA-256 recorded at
`/gpfs/kjhan/LagForce/Kepler_Q1_BallisticResidual_T0825_20261027/inputs.sha256`.
The complete [machine-readable receipt](KEPLER_Q1_NATIVE_BALLISTIC_T0825_2026-10-27.json)
passed finite, fixed-mass, W=0 and positive-state checks.

The control moves the *same initial generators* along their softened exact
circular orbits to `t=.825`, reconstructs an independent W=0 Voronoi disk,
and compares native cell density to the resulting fixed-mass density. The
initial `0.5<=r<2` disk contains 48,680 q1 cells.
An additional omitted-periodic-image certificate at this exact time gives
minimum `edge_distance/(2*circumradius)=46.35` across those polygons, so the
independent nonperiodic reconstruction is valid for the disk.

| Initial radial range | Native mean absolute density change | Exact-orbit geometry change | Native/geometry pointwise residual | Residual / native change |
|---|---:|---:|---:|---:|
| 0.5–0.75 | 9.12990e-4 | 9.13501e-4 | 9.14571e-7 | 0.100% |
| 0.75–1 | 8.49647e-4 | 8.50168e-4 | 6.73391e-7 | 0.079% |
| 1–2 | 4.35240e-4 | 4.35419e-4 | 4.07496e-7 | 0.094% |
| Entire 0.5–2 disk | 5.23520e-4 | 5.23767e-4 | 4.80819e-7 | 0.0918% |

All listed values except the last column are initial-mass-weighted mean
absolute relative quantities. The last column compares two such norms; it
is **not** a signed error budget or a proof that exactly 99.91% of the
variance has one cause. Native positions differ from exact circular orbits
by only `1.11215e-8` mass-weighted in the disk. At this early time, the
native density field is therefore very close to the exact-orbit tessellation
prediction. An early shear-induced density error is dominated by the
discrete mesh geometry, not visible loss of circular generator motion.

The [strict native gate](KEPLER_Q1_HALO12_T0825_STRICT_2026-10-27.json)
also gives disk absolute entropy-log deviation `1.53318e-4` but signed
mass-weighted entropy-log change `-8.45393e-6`. The absolute deviation does
not by itself prove positive net heating. The outer periodic-seam region
already has much larger density/pressure deviations, and this `t=.825`
decomposition cannot establish its later effect on the disk. The matching
[q2 native-versus-ballistic audit](KEPLER_Q1_Q2_NATIVE_T0825_RESULT_2026-10-27.md)
is now complete. Continue both physics and boundary diagnostics through
`t>6` before choosing a physical correction.
