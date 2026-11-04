# Linear modes and known-wave response of the constrained-weight candidate

This bundled diagnostic continues the coefficient-error stage in `312bd4d`.
Production hydrodynamics is unchanged. No fluid was advanced, no additional
native coefficient reference was solved, and no limiter or solver tuning was
introduced. The complete dimensional accuracy and GIZMO throughput goals remain
open. The current candidate is retired after the uniform-pressure equilibrium
control documented in `UNIFIED_HYDRO_EQUILIBRIUM_RETIREMENT_2026-11-05.md`.

## Periodic symbols retain image phases

The original irregular 1D/2D/3D primitive geometries have 16, 9, and 27 cells.
The previous independent finite-matrix reference definition reconstructs the
small coefficient fields. The 2D/3D fields agree with the archived doubled
pattern after the expected factor-of-four weight scaling to 1.0e-16.

For each periodic face, the Bloch phase uses the actual image shift
`r_ij-(x_j-x_i)`. The finite complex symbols obey

    C(k) = B(k) + L(k) R(k).

The uniform-state test sets rho=1 and gamma P=1, giving sound speed one.
Cell masses are the actual positive cell volumes and remain fixed. With
`H=diag(1/V)` and `M=diag(m)` the positive acoustic pressure block is

    A_P = H^(1/2) C M^(-1) C^dagger H^(1/2).

Complex symbols use the conjugate transpose. Nonnegative squared frequencies
follow the adjoint construction and do not alone establish physical accuracy.
The pressure-plane branch is selected by overlap with the volume-weighted
continuum pressure wave. No global maximum over all wavevectors is claimed.

Thirty-five symbol cases cover axial and diagonal directions at fixed wave
norms 0, .01, .02, .04, pi/4, pi/2, and pi. In 1D C=B exactly. Both 1D operators
retain the two previously known zero-wave pressure null modes. The 1D plane
branch is nearly degenerate with the alias branch, so branch overlap is recorded.
In 2D/3D the zero-wave symbols have one pressure null mode for both B and C.
No additional pressure null mode is detected at the sampled nonzero wavevectors.

At wave norm .01, the candidate axial sound speeds are 0.999998151 in 2D and
0.999998188 in 3D. The errors shrink as the wave norm decreases. The original
operator also has the correct long-wave acoustic speed. The candidate improves
the transverse, divergence-free response rather than correcting a wrong bulk
sound-speed limit.

## Fixed physical wavelength and repeated-pattern refinement

The primitive geometry is rescaled analytically to repeat over the same
physical box. Volumes scale as length^D and geometric coefficients follow
their exact length dimensions. This control uses original native cells and
does not substitute a different tessellation or perform additional evolution.
The archived doubled-pattern coefficient check validates the scaling.

False transverse compression is normalized by `|k| amplitude total_volume`.
The exact divergence is zero. The ratio C/B is not used when the B response
is nearly zero. Axial fixed-wavelength results are:

| Dimension | Global side | B false compression | C false compression |
|---|---:|---:|---:|
| 2 | 12 | 0.07942 | 0.01142 |
| 2 | 24 | 0.08063 | 0.00569 |
| 2 | 48 | 0.08093 | 0.00284 |
| 2 | 96 | 0.08101 | 0.00142 |
| 3 | 12 | 0.06624 | 0.01141 |
| 3 | 24 | 0.06641 | 0.00572 |
| 3 | 48 | 0.06639 | 0.00285 |
| 3 | 96 | 0.06637 | 0.00142 |

Diagonal controls also decrease approximately linearly for C. The measurements
support one necessary transverse-consistency property on the repeated irregular
patterns. The results do not prove convergence for an evolving deformed mesh.
The highest frequencies sampled on the mild primitive patterns are slightly
smaller for C than for B. A regular or mildly irregular pattern therefore does
not capture the practical stiffness found in the evolved native state.

## Known-wave probes on the evolved 128² geometry

The accepted original-operator Gresho checkpoint at t=.5 supplies geometry,
mass, and background pressure. Only the previously saved reference coefficients
are used. Axial and diagonal sinusoidal velocity/pressure probes have amplitude
one and integer wave components 4, 8, and 16. An additional axial mode 64 is
marked as a cell-scale diagnostic rather than an accuracy gate.

The divergence and pressure-gradient references are exact polygon integrals
of the specified Fourier wave. Straight-face integrals have closed forms from
the original face centroid, length, normal, and tangent. The divergence theorem
then gives the cell integral. Independent uniform-square checks agree with the
analytic product of sinc functions to 3.6e-16. This is an analysis calculation,
not runtime subcell integration in the hydrodynamic code.

Velocity is sampled at generators. Pressure is tested both as a generator
sample and as the exact cell average. Pressure representation is not changed
silently. Errors are normalized by `|k| amplitude total_volume` and use the
same exact integrated gradient for both pressure inputs.

| Axial wavelength | Quantity | B error L1 | C error L1 |
|---:|---|---:|---:|
| 1/4 | transverse false compression | 0.05246 | 0.00902 |
| 1/4 | longitudinal divergence | 0.05431 | 0.01158 |
| 1/4 | pressure force, generator samples | 0.01071 | 0.02502 |
| 1/4 | pressure force, exact cell averages | 0.04538 | 0.16160 |
| 1/8 | transverse false compression | 0.05320 | 0.01941 |
| 1/8 | longitudinal divergence | 0.06048 | 0.02806 |
| 1/8 | pressure force, generator samples | 0.02586 | 0.05676 |
| 1/8 | pressure force, exact cell averages | 0.05425 | 0.17180 |
| 1/16 | transverse false compression | 0.05746 | 0.04000 |
| 1/16 | longitudinal divergence | 0.08473 | 0.07249 |

The candidate improves smooth volume response while worsening pressure force
in the reported cases. Maximal local errors and shorter diagonal wavelengths
can worsen even when an L1 volume error improves. A larger response to the
previous evolved-state deviation is consequently neither a sufficient reason
to reject C as unstable nor evidence of a fully improved physical operator.
The pressure representation does not explain away the measured force error.

## Practical stiffness is strongly geometry dependent

For the same positive background pressure and actual cell masses, a bounded
sparse diagnostic estimates the largest eigenvalue of the positive block for
each operator. The cap is 256 matrix applications and 15 seconds per operator.
The comparison concerns the same positive acoustic block as the earlier audit.
The moving Gresho background and remaining geometric terms are not included.

| Operator | Highest estimated frequency | Ritz applications | Relative Ritz residual |
|---|---:|---:|---:|
| B | 3033.19 | 42 | 6.35e-10 |
| C | 150319.65 | 18 | 5.28e-14 |

The frequency ratio C/B is 49.5583. The corresponding inverse-frequency scale
is about 0.0202 of the B scale. This is a substantial explicit-timestep cost
warning, not a growth rate, a measured trajectory instability, a complete
KDK stability limit, or a measured fiftyfold PDE runtime.

The surprising candidate value was checked by a separate face reassembly.
The verifier computes each 2x2 fit inverse directly rather than reading the
captured native inverse. Thirty-two power iterations per operator stay inside
the same combined 256-application/15-second limits. The candidate Rayleigh value
agrees with the Ritz result to 1.4e-15 relatively, with relative residual
4.24e-14. The B power iteration is still 1.18% below the converged B Ritz value
and is explicitly not reported as a converged B verification.

The candidate mode has pressure participation 2.74 cells and kinetic
participation 2.90 cells. Most pressure norm lies near IDs 13617--13619 in the
shearing annulus around R=.35. The cell volumes are near the mean volume and the
minimum neighbor distance is about .326 nominal cell spacings. The large
frequency is therefore not explained solely by a vanishing cell volume.
The associated local weight-rate stencil coefficients reach about 6.6--9.1
nominal cell spacings. Localization is measured on a frozen state, not a
guarantee that a future candidate trajectory reaches the same state.

## Current interpretation and provenance

The candidate passes the sampled affine/long-wave/transverse necessary checks.
The evolved-geometry pressure force and acoustic stiffness remain practical
problems. A more accurate or faster coefficient solve would not remove the
reference-operator stiffness measured here. No coefficient optimization or new
long PDE trajectory is justified by the completed bundle. The earlier CPU
cost/tolerance failure is preserved.

The first analysis is retained at
`/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/material_volume_modes_v1`.
An adapter divided a two-component velocity vector by one mass per cell and
stopped on a shape error before a native spectral result. The corrected analysis
repeats the mass for both components and passes an explicit two-cell matrix
check with unequal masses. No PDE or native coefficient solve was repeated.

Authoritative results are in `material_volume_modes_v2/results`, with the
independent check in `material_volume_modes_v2/verification/result.json`.
Compact source, commands, checks, and results are archived under
`Research/ModeResponse_20261105`. The manuscript and accepted trajectories keep
their previously verified state. The next method decision must retain the
requirements of self-consistency, simplicity, no steering, preferred fixed
mass, and matched GIZMO runtime at most twice the reference.
