# Kepler polar IC: material and Voronoi face-motion gap at t=0

The native W=0 Kepler solver fixes each cell's mass and computes its density
from its *geometric* Voronoi area. The active CUDA binary was verified to
contain `GFS_VARPRESS_W0`: its pressure force and energy work use the
geometric face area rate, with an adjoint force. Ordinary phase-1 AV5 faces
also use the geometric speed in their reversible work. **HLLC contact speed
is not substituted for that geometric speed in the active variational work
law.** This makes the thermodynamic update geometrically consistent in the
smooth W0 branch; it does not make a Voronoi cell a material volume under
Kepler shear. The calculation below compares the face's geometric normal
speed `w_g` with the **exact analytic Kepler material velocity** at its
centroid. It measures a kinematic source of fixed-mass density variation,
not a work-law defect, HLLC flux, entropy source, or late-time heating.

The [read-only auditor](audit_kepler_face_motion_gap_20261027.py) constructs
ordinary Voronoi faces from the full native q1 and q2 t=0 snapshots. For a
face between generators i and j, normal `n=(x_j-x_i)/d` and centroid `c`,
the exact W=0 bisector motion is

`w_g = n·(v_i+v_j)/2 - (v_j-v_i)·[c-(x_i+x_j)/2]/d`.

It compares this with `v_K(c)·n`, where `v_K` is the analytic softened
Kepler field. All reported faces touch an initial `0.5<=r<2` disk cell;
same-shell and cross-shell rows require *both* cells to be in the disk.

| Measure at t=0 | q1 | q2 | q2/q1 |
|---|---:|---:|---:|
| Same-shell `sum L |v_K·n-w_g|` | 0.172045 | 0.069341 | 0.403 |
| Cross-shell `sum L |v_K·n-w_g|` | 0.044410 | 0.036171 | 0.814 |
| Cross-shell face-length-weighted mean `|v_K·n-w_g|` | 5.553e-5 | 6.749e-5 | 1.215 |
| Disk-cell mass-weighted mean `|d ln A_g/dt|` | 0.0041387 | 0.0020369 | 0.492 |

The **same-shell** contribution is larger than the cross-shell contribution
in both layouts. q2 reduces the integrated cross-shell mismatch through
shorter total contact length even though its typical *per-length* mismatch
is larger. Thus a model attributing the geometric volume discrepancy only
to adjacent-shell face shear is incomplete. The signed face contributions
also cancel differently from the absolute face sums, so neither row alone
predicts the final density or entropy error.

The geometric per-cell area rate was independently checked by rebuilding
the full Voronoi diagram at exact-orbit times `t=±1e-4` and differencing
areas. The mass-weighted mean absolute rate discrepancy divided by cell area
is `5.65e-7` in q1 and `2.97e-7` in q2, much smaller than the geometric
rates above. Native t=0 stored areas and independently reconstructed areas
agree to a maximum relative difference of `3.03e-8` in disk cells. A few
topology-sensitive outliers have a larger finite-difference discrepancy,
so the check supports the aggregate rates rather than every cell's derivative.
The full numbers are in the [validated receipt](KEPLER_FACE_MOTION_GAP_T0_V2_2026-10-27.json).

**Next native attribution:** At matched q1/q2 checkpoints beyond t=6,
separate (i) geometric `dV/dt` and the native pressure-work/entropy residual,
(ii) face force and torque by same-shell/cross-shell and radial bin, and
(iii) the outer periodic seam from the physical disk. The optional HLLC
contact speed may be logged to quantify *material* flux through the moving
Voronoi face, but `P_*L(v_*-w_g)` is **not** the active variational scheme's
missing heat term. Only a measured native entropy/force residual can motivate
weight evolution, mass exchange, or a cell-spin operator.
