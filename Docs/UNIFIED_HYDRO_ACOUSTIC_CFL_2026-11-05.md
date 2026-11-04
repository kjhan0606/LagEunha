# Complete pressure Hessian, KH endpoints, and a bounded timestep experiment

This stage continues `078a226` in the same LagEunha checkout. The goal remains
contact causality, consistent 1D/2D/3D hydrodynamics, full 2D validation, simple
algorithms, and matched GIZMO wall time <=2x. It is not complete. All calculations
use the original native program; Python performs orchestration/analysis only.
Builds and trajectories are under `/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1`.

## KH completion and direct field audit

All three periodic KH cases now reached exactly t=2 with exit zero. Their
per-ID masses are bit-identical to birth, W and weight momentum remain zero,
U/V are positive, domain-volume sums agree to roundoff, and material entropy
increments are positive. Density is always m/V. The same common pressure,
face-angular velocity fit, central BE viscosity, and global KDK were used.

| Resolution | Relative energy change from birth | Mode amplitude at t=2 | Transverse kinetic energy | Material overlap width |
|---|---:|---:|---:|---:|
| 64 by 128 | 1.67269e-4 | 0.17206 | 0.02930 | 0.51950 |
| 128 by 256 | 7.19548e-4 | 0.27187 | 0.07673 | 0.64134 |
| 256 by 512 | 2.46529e-4 | 0.30621 | 0.10121 | 0.66976 |

The 64 by 128 case continued its accepted t=0.05 checkpoint. Higher resolutions
started at birth. Their physical layer width and perturbation are identical.
The overlap width measures material interleaving, not molecular diffusion.
These are evolution diagnostics, not independent-reference error norms. The
energy errors are not monotone in resolution under the original CFL.

The density figure in the evidence archive samples exact periodic W=0 Voronoi
regions with a nearest-site query and assigns each region its actual m/V.
It does not interpolate density between particles. All three panels use the
same [0.9,2.2] scale; no cell density lies outside that colour range.

## Complete fixed-entropy pressure derivative

For fixed masses and material entropy,

    H = B^T diag(gamma P/V) B - sum_i P_i Hessian(V_i).

The native diagnostic displaces generator positions along a mass-normalized
direction, rebuilds the original geometry, applies
`U_i = U_i0 (V_i0/V_i)^(gamma-1)`, refreshes that state on the GPU, and differences
the resulting pressure acceleration. Six displacement sizes are used. It takes
no accepted evolution steps and supplies no production force correction.

Along the dominant positive-acoustic direction, the full directional curvature
gives native dt times frequency of 0.547 at initial 64², 1.741 at 64² t=0.5,
2.0367 at 64² t=1, and 2.1134 at 128² t=0.5. Force derivatives converge as the
perturbation shrinks; quadratic energy differences agree before cancellation
dominates. At t=0.5, independent-direction symmetry residuals are 1.17e-10 and
7.69e-11 for 64²/128². Uniform-translation derivatives are about 1e-11 of the
acoustic derivative norm. These are directional checks, not a full spectrum.
Negative eigenvalues have not been bounded or ruled out.

The values above 2 identify a missing fast-pressure-mode restriction in the old
cell-size CFL. The threshold is the fixed positive quadratic KDK criterion;
it is not a stability theorem for the nonlinear, moving, dissipative system.

The first diagnostic version omitted the post-geometry thermal-state refresh;
its force curvature disagreed with the energy curvature. That version is invalid
for the full Hessian and is preserved as failed diagnostic evidence. The corrected
probe uses the same refresh as the production final pressure kick. No previous
evolution trajectory is reclassified using the invalid probe.

## Cheap dimension-common bound, still an opt-in candidate

`Exam/gfs_acoustic_cfl.h` implements a rotation-invariant block bound:

    T_ij = sqrt(gamma P_i/V_i) norm(B_ij) / sqrt(m_j),
    lambda_positive <= max_row_sum(T) max_column_sum(T).

The candidate ceiling is `2*courant/sqrt(bound)`. It requires only a face sweep,
not an eigensolve or additional mesh construction. Dense acoustic matrices from
actual project geometry in 1D/2D/3D verify the bound against independent NumPy
eigenspectra. Bound/exact maximum ratios are 1.315, 3.236, and 5.212 for the
small tested geometries; rotation changes the bound by at most 3.33e-16.

The full pressure-weighted geometric Hessian is not covered by this bound.
The opt-in native adapter (`GFS_UNIFIED_ACOUSTIC_CFL=1`) currently supports only
periodic one-rank 2D. The common algebra is tested in all three dimensions;
operational 1D/3D timestep integration is not yet implemented or qualified.
The default timestep and all pressure/viscosity laws remain unchanged.

The first adapter attempted to read a host face CSR before downloading it and
aborted before any accepted step. The corrected version explicitly downloads
and validates the device CSR. Both attempts are retained with separate immutable
builds and run directories; existing runs were not overwritten or restarted.

## Gresho timestep comparison and spatial error separation

The first four rows below start from the identical accepted t=0.05 state and
reach exactly t=1. The final row checks the candidate from the original IC.

| Timestep | Steps | Relative energy change from start | Vector velocity L1 | Segment viscous heat |
|---|---:|---:|---:|---:|
| Original | 515 | 1.14040e-3 | 0.07656 | 0.03832 |
| Original / 2 | 1000 | 1.49544e-5 | 0.07171 | 0.02843 |
| Original / 4 | 1996 | 4.07731e-6 | 0.07318 | 0.02898 |
| Acoustic bound | 1721 | 4.92564e-6 | 0.07298 | 0.02832 |
| Acoustic bound, birth | 1742 | 5.20645e-6 | 0.07756 | 0.02945 |

The candidate controls energy error while retaining roughly one third loss of
initial kinetic energy. Its velocity error is comparable to the uniform small
step controls. Bound evaluation costs about 1.2% of the measured native phase,
but the step count is 3.34 times the original continuation. Timings were taken
under different concurrent loads and do not constitute a matched speed result.
The candidate is not promoted to the default and does not solve vortex accuracy.

On identical saved geometry/thermal states, replacing only generator velocities
by the exact stationary Gresho field leaves pressure acceleration identical and
changes instantaneous heat rates as follows:

| State | Evolved velocity | Exact velocity | Ratio |
|---|---:|---:|---:|
| 64² t=0.5 | 0.09079 | 0.005587 | 16.3 |
| 64² t=1 | 0.05374 | 0.005945 | 9.0 |
| 128² t=0.5 | 0.12856 | 0.001867 | 68.9 |

Exact affine core rotation gives roundoff-level direct heat. Direct closure
truncation remains around the piecewise-profile transitions, while evolved
noncircular motion substantially amplifies heat in smooth regions. This does
not by itself assign all of that motion to a single cause.

## Next work and evidence

The next causal test should distinguish geometric Voronoi volume work from
physical divergence, and pressure-field noise from force discretization, on the
same saved states. Timestep reduction alone is insufficient. Do not add another
shear limiter/restoring force or launch new long Kepler trajectories from this
result. After selecting a consistent treatment, qualify its dimensional and
shock contracts, remaining 2D suite, and matched GIZMO cost.

`Research/AcousticCFL_20261105` holds the compact results, failed-attempt records,
build commands/manifests, source versions, and figures. Binaries and native
trajectories remain on GPFS with hashes. Builds are partial relinks; a clean
native build remains outstanding. The manuscript retains its earlier verified
baseline; these experiments have not replaced its published-method results.
