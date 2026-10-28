# Controlled strip shear: spatial projection mixing decreases with refinement

Analytic 2D unit-width strip cells, rho=1, v=(y,0), fixed generators, bounded
y domain[0,1]. Refine only y with N16,32,64,128. Alternate interior face shifts
of magnitude .15h, zero exterior shift; four increments out and four back.
Construct adjacent weight differences from delta_face=(Wi-Wj)/(2h), verified
exactly in this run. Weight amplitude scales as h^2: fixed deformation in cell
units, not a fixed absolute weight function across resolutions.

Two tests start with the same piecewise-constant sampled shear and zero
internal kinetic variance: constant projection and affine moment recovery.
A third exact-affine initialization is a separate consistency control with
its own resolved internal kinetic energy. It is not an energy-free conversion
of the point initial state. All integrals are analytic polynomial integrals.

| N | Point initialization, affine remap residual K | Total internal K |
|---|---:|---:|
| 16 | 1.49409e-4 | 2.54428e-4 |
| 32 | 3.81166e-5 | 6.53659e-5 |
| 64 | 9.62471e-6 | 1.65613e-5 |
| 128 | 2.41812e-6 | 4.16782e-6 |

Last-pair exponents are1.9929 for residual K and1.9905 for total internal K.
Constant projection also approaches h^2 internal-energy scaling, with larger
coefficient. In this controlled geometry, time-subdivision plateaus therefore
do not imply spatial nonconvergence.

The exact-affine control preserves generator velocities and produces zero
residual K to printed precision at all resolutions; its initial nonzero
affine K scales as h^2. Affine remapping preserves represented Lz, mass,
momentum and total kinetic moment to floating precision. Constant projection
does not preserve represented Lz, consistent with the earlier geometry audit.

## Boundary caveat

Whole-domain point velocity RMS scales as h^1.5 for both sampled-shear
projections, but the interior .25<y<.75 RMS is at most about8.01e-14. Symmetry
cancels the interior mean errors and boundary effects dominate the global norm.
Do not report h^1.5 as a general scheme order or a Kepler convergence result.
The v2 audit adds the explicit interior check; the initial result is retained.

## Consequence

The retained residual energy can be a spatial representation error that
decreases with h, rather than evidence of a new physical heat source. It must
still not be silently thermalized. This strip result does not establish the
same scaling on irregular Kepler cells with density jumps and curved flow.
Before production integration, test irregular geometry/velocity curvature
and retain the original cold point-state energy convention; the exact-affine
control cannot be substituted for the authorized Kepler initial condition.

No physical time, pressure, gravity, AV, coupled CUDA/OpenMP/MPI or long-time
stability is tested. Original simulations were not modified.

Evidence: `canonical_shear_remap_refinement_2026-10-18.json`.
Script: `audit_shear_remap_refinement.py`.
Remote: `/gpfs/kjhan/LagForce/Codex_shear_remap_refinement_v2_20261018`.
