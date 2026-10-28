# Smooth pressure moment rate and compatible internal-energy work

The experimental header `Voro/Laguerre/moment_pressure_rate.h` supplies two
functions in the existing17moment layout:

- `mpr_face`: pressure traction, position-weighted traction, velocity-weighted
  traction and total-energy pressure flux for affine traces on a straight face.
- `mpr_volume`: the integration-by-parts terms delta_ab integral(P) in mixed
  moments and (G+G^T) integral(P) in velocity second moments.

Adding only face terms is insufficient. Combined rates obey
Udot=Edot-tr(Pidot)/2=-integral(P div(v)) for smooth affine fields.

The actual t0 closed polygon geometry was tested with manufactured affine
pressure and velocity, including nonzero divergence. Independent three-point
face quadrature plus direct volume pressure-gradient integration checks all
moment components. On11,498cells and68,986faces, CPU1/OpenMP4/CUDA outputs are
bitwise equal with FMA disabled. Maximum assembled error divided by volume is
9.43e-16; the internal-pressure-work identity residual is4.71e-16.

The first diagnostic mixed closed polygon integrals with generator-direction
normals and failed its1e-10 check at2.07e-10. The corrected diagnostic uses
the actual polygon-edge normals throughout. Four nearly vanished edges of
length about4e-11 have normal differences .51–.54; for A/d>1e-6 the maximum is
4.26e-13. This does not establish a material Kepler error source. Original
production normals and forces were not modified.

Evidence: `canonical_pressure_rate_backend_2026-10-18.json` and
`canonical_pressure_normal_audit_2026-10-18.json`. Failed preparation directory
is preserved; completed outputs are in
`/gpfs/kjhan/LagForce/Codex_pressure_rate_backend_v2_20261018`.

This is not a discontinuous interface/Riemann pressure law. It does not yet
include gravity, evolving nonlinear reconstruction, positive-state guarantees,
MPI, coupled RK evolution, or the approved original fixed-mass reduction.
No production call site was added and no long-Kepler success is claimed.
