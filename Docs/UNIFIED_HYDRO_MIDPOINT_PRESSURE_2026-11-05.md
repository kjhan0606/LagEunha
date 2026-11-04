# Midpoint-offset divergence stops at the pressure-consistency gate

The previous local `psi_i G_i` construction was retired in `de66900`. A bundled
Astra plan selected one finite eligibility test for a different midpoint-offset
face divergence. The new test is complete and rejects the proposed operator
under the preregistered pressure gate. No production force, geometry, pressure,
or timestep changed. The complete user objective remains active.

## Definition and difference from the historical reconstruction

Use the existing common face-angular least-squares velocity gradient, original
native face areas and centroids, and periodic image displacements. Define

    r_ij = x_j(image)-x_i.
    o_ij = x_face-(x_i+x_j(image))/2.
    u_face = (v_i+v_j)/2 + (G_i+G_j) o_ij/2.
    (Dv)_i = sum_faces A n_i dot u_face.
    F = D^T P.
    Udot_i = -P_i (Dv)_i.

The force transpose includes the transpose of the gradient reconstruction.
Using an unrelated central pressure force would define a different operator.
The actual geometric volume rate remains Bv until a compatible weight law is
implemented. No D-only thermal update is installed in the fluid code.

The Nov2 reconstruction extrapolated each cell velocity separately to the
face. With the same fit, the old expression has an additional term

    u_face_old-u_face = (G_i-G_j) r_ij/4.

Only the formula difference is replayed with the current common fit. The old
Kepler campaign is not rerun. The historical fit used different angular weights,
so the matched-fit formula control is not presented as an exact reproduction
of the historical numerical experiment.

At W=0 in 1D, the face offset is zero and D=B algebraically for every velocity
field. Existing native interval geometry gives matrix differences below
3.2e-15. The old extra term would change the 1D operator by about 0.269 in the
same matrix-entry norm. The proposed change is therefore a distinct candidate
with an unchanged 1D reduction, despite sharing the reconstructed-Gauss family.

## Completed algebra before the pressure gate

The six original native geometries are 1D 16/32, 2D 3²/6², and 3D 3³/6³.
The following checks pass without a new tessellation or coefficient solve:

- Direct face action agrees with the sparse D matrix.
- Independent face/gradient transpose agrees with D^T.
- The old-minus-new term agrees with the explicit gradient-difference flux.
- Translation, constant pressure, force/work adjoint, and net-volume weight-RHS
  compatibility hold to roundoff.
- Affine dilation, shear, rotation, and general strain use exact image
  continuation and an explicit periodic cut-work ledger.

The maximum basic algebra defect is 5.7e-16. The largest affine/cut defect is
1.93e-14. Compatibility of `L Wdot=(D-B)v` is only a necessary zero-sum condition.
No weight solve, geometry positivity, or finite-time integration is qualified.

## Pressure representation and fixed-pattern refinement

Pressure is tested as both a generator sample and an exact cell average of a
linear field. The primary gate was fixed before the run to the cell-average
thermodynamic interpretation of `P=(gamma-1)U/V`. The nodal test separately
measures the loss of the original B patch property. Neither representation is
selected afterward to obtain a pass.

Exact linear cell averages use the original polygon/polyhedron centroid.
The centroid follows the first-moment divergence theorem on actual native
faces. Each original primitive pattern is tiled and rescaled analytically.
The rescaled coordinates and volumes agree with the archived doubled native
geometry to 1.2e-16. This uses exact geometric scaling, not a substitute mesh.

One complete central primitive pattern is sampled at both resolutions, with
both reconstruction layers excluded from the periodic seam. Sample types
and their relative geometry are identical. The primary gate stops when an
error above 1e-9 retains at least 99% of its amplitude after factor-two
refinement. This is a detector of the measured consistency floor, not a chosen
allowable physical error for the final solver.

For a unit pressure gradient along x, the force-density L1 errors are:

| Dimension | Coarse/fine global side | B, nodal P | D, nodal P | B, cell-average P | D, cell-average P |
|---|---|---:|---:|---:|---:|
| 1 | 48 / 96 | roundoff | roundoff | 0.0172832 | 0.0172832 |
| 2 | 9 / 18 | roundoff | 0.0410618 | 0.0294977 | 0.0415866 |
| 3 | 9 / 18 | roundoff | 0.0484122 | 0.0187311 | 0.0398783 |

Every nonzero entry retains the same amplitude under the matching refinement
to reported precision. The y/z primary tests also fail to decrease. All six
dimension/axis primary comparisons fail. D improves affine volume transport
but loses the original nodal pressure-gradient patch in 2D/3D and worsens the
cell-average pressure result. The first pressure gate therefore rejects D.
Bloch, N128 stiffness, weight realization, and candidate evolution are not run.
This is a result for the selected operator and repeated irregular patterns,
not a theorem excluding every weak Euler discretization.

## Independent 1D cause of the inherited cell-average defect

The 1D comparison exposes an inherited issue because D=B. For an unweighted
Voronoi interval, its volume and actual centroid are

    V_i = (x_(i+1)-x_(i-1))/2.
    c_i = (x_(i-1)+2 x_i+x_(i+1))/4.

For a linear pressure field with gradient b, the exact cell average is
`P_i=a+b c_i`. The original operator gives

    F_i = (P_(i-1)-P_(i+1))/2.
    (F_i+V_i b)/V_i
      = b [(c_(i-1)-x_(i-1))-(c_(i+1)-x_(i+1))]/(2 V_i).

The centroid offsets and V scale by the same refinement factor. Their
force-density ratio therefore does not decrease on the repeated irregular
pattern. Nodal pressure instead gives the original exact linear patch.

An independent C calculation evaluates the intervals, closed-form defect,
and existing `gfs_variational_face`/`gfs_variational_pressure` helpers. At 48/96
cells the cell-average L1 errors are 0.01728315099810 at both resolutions.
The closed-form and native-helper defects are below 2.2e-14. The independent
C result agrees with the sparse operator result to better than 1e-12.

The earlier stationary-contact initialization cause remains valid. Constant
pressure has no nodal/average difference, and the 33 contact controls concerned
the initial mass/volume/thermal projection. The new probe concerns a pressure
gradient and does not replace the earlier causal result. A point-quadrature
model and a literal cell-average model require separate definitions; the
pressure-location result is not a license to relabel the stored state silently.

## Evidence and remaining scope

Run root is `/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/midpoint_divergence_v1`.
All diagnostic processes are terminal with exit zero. The scientific gate
fails explicitly in `results/result.json`. The independent C proof is in
`pressure_location_1d`. Compact evidence is in `Research/MidpointPressure_20261105`.
Original native geometry, code, manuscript results, and trajectories are
preserved. No PDE step, new global solve, or multidimensional mesh construction
was performed. Self-consistent state representation, an affordable common
1D/2D/3D Euler method, the remaining 2D suite, and matched GIZMO speed remain open.
