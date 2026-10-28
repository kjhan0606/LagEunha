# Sparse q2 A10: provenance of the outer periodic-seam disturbance

The [read-only auditor](audit_kepler_seam_provenance_20261027.py) tags
native cells by their **initial** distance to a box face and initial
radius, then follows their generator centres and fixed masses by cell ID.
It checks complete finite positive W=0 snapshots, exact fixed mass and
periodic area 36 before measuring crossings. A centre crossing is not a
face-integrated fluid mass flux, and a force or pressure disturbance can
propagate through neighbors without any tagged centre crossing.

At the first common saved output, actual `t≈1.650242`, the independent
[wide](KEPLER_SPARSE_WIDE_SEAM_PROVENANCE_T165024_2026-10-27.json) and
[narrow](KEPLER_SPARSE_NARROW_SEAM_PROVENANCE_T165024_2026-10-27.json)
receipts agree closely:

| Initial-source diagnostic | Wide | Narrow |
|---|---:|---:|
| Cells initially within 0.2 of a box face | 152 | 152 |
| Those now inside `r=2.75` | 20 | 20 |
| Their fixed mass now inside `r=2.75` | `0.00493599` | `0.00493599` |
| Those now inside guard radius `r=2.35` | 0 | 0 |
| Those now inside disk edge `r=2` | 0 | 0 |
| Minimum current radius of initial `r>=2.75` cells | `2.69024136` | `2.69024129` |
| Radial gap: that minimum minus outermost original-disk site | `0.70150352` | `0.70150301` |
| Initially seam-tagged cells' absolute density change | `0.649392` | `0.649393` |
| Initially seam-tagged cells' absolute log-entropy change | `0.00171845` | `0.00171841` |
| Initial `2.35<=r<2.75` buffer's absolute density change | `0.090354` | `0.090226` |

The seam region is badly deformed already, while no initially seam-tagged
generator has entered the guard or disk and the outer-to-disk generator
gap remains about 0.70. The buffer's density change is substantial even
though its generator-centre radial speed is only about `5e-5` by the
[independent radial audit](KEPLER_SPARSE_A10_T10_PARALLEL_2026-10-27.md#first-common-saved-checkpoint).
This does **not** establish that the disk is dynamically isolated: Voronoi
faces and pressure work can transmit a disturbance before centres cross
these radii. It does show that simple advection of the tagged seam cells
is not yet the cause of the disk's early density residual.

Repeat the same provenance audit on later complete `t>6` and `t=10`
snapshots, together with disk entropy, native-versus-exact-orbit density
residual, torque and radial velocity. If the outer source crosses the
guard or the disk while those errors accelerate, run a matched isolated
boundary control before promoting a periodic sparse solution as a GIZMO
Kepler benchmark. An unchanged centre gap alone is insufficient to rule
out face-mediated coupling.
