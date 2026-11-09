# Compatible reconstructed pressure work: useful identity, incomplete closure

Follow-up to KH_WEIGHT_TRANSPORT_SOLVE_2026-11-06.md. No native time evolution
or production force changed. Uses the same32x64/64x128 native exact row-shear
geometries at t=1, not another long KH simulation.

## Definition and code reuse

The existing dimension-generic C/CUDA-compatible gfs_velocity_fit.h supplies
the small least-squares inverses. fit.c calls its actual init/add_face/solve
functions; no fit fallback or tunable regularization is introduced. In2D
the fit normal matrix uses area/distance^3 times r r^T. Note that the header's
current production use is viscosity fitting, not primitive extrapolation;
the proposed use here is a diagnostic, not an already-qualified production use.
This weighting differs from the older inverse-distance-squared offline screen
in WEIGHT_VOLUME_PROJECTION_SCREEN_2026-11-06.md.

Extrapolate velocity from each side to the face centroid and take the symmetric
average normal velocity T v. With signed face-area incidence D, define C=D T.
If the material target rate is C v, weights would need to satisfy

    L Wdot = (C-B) v.

Pairing momentum force F=C^T P with internal energy rate -P_i (Cv)_i gives
exact instantaneous pressure-work cancellation. The force is local once the
fit is known, although finding compatible geometric weights is still global.
This is not the retired psi-times-strain model.

## Checks and limitation

Uniform-pressure force is below1.19e-17. Interior affine shear, rotation and
expansion reproduce their expected divergence to7.96e-14; periodic boundary
wraps are excluded from affine fields that are not globally periodic.
Pressure-work cancellation is below6.32e-15 with nonuniform diagnostic pressure.
These are consistency/adjoint checks, not proof of physical accuracy.

For the actual nonlinear smooth shear U(y), native B v/V RMS is below2.11e-14.
The reconstructed target C v/V RMS is0.001861 at32 and0.005791 at64. It would
generate immediate relative pressure RMS rates0.003102 and0.009652 per time.
This is a spatial truncation residual in a steady continuous field. Nonzero
truncation alone does not invalidate a consistent method, and these two
shear phases do not establish its asymptotic order. It does mean that affine
exactness and energy conservation must not be advertised as preserving this
nonlinear shear exactly.

More fundamentally, on a regular square lattice with x-only variations,
the same symmetric face rule has divergence symbol

    i sin(theta) [1 + (1-cos(theta))/2] / h.

Its two-cell alternating mode theta=pi is still zero. A direct64-column
alternating-velocity calculation gives exactly zero face divergence. By the
same translationally invariant adjoint, the matching alternating pressure
mode remains unseen by C^T. Raising the reconstruction order in this symmetric
central rule therefore does not remove the checkerboard null. It cannot be
accepted as a complete explanation or cure of the observed short-wave issue.

## Decision

Retain the compatible force/work construction as a required design condition,
but do not port this symmetric reconstruction alone as the solution. The
common method must also provide an acoustic response to alternating pressure
while avoiding false dissipation of resolved shear. The actual1D material
Riemann method supplies an acoustic interface response; replacing it by a
symmetric average in2D would not meet the user's common-method requirement.
How to couple that response to constrained Laguerre motion and energy remains
unfinished. No sensor, viscosity switch, or new physical parameter was added.

Artifacts: Research/MaterialVolume_20261105/kh_reconstructed_work_20261106_v1/:
fit.c, audit.py, native fit inverses/condition reports, results.json, analysis.log.
No3D work started. Full2D accuracy and <=2x GIZMO runtime remain unverified.
