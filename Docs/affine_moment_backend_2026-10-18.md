# Existing moment layout reused for face transport

`Voro/Laguerre/affine_moment_flux.h` reuses `wm_atom` from
`weight_remap_moments.h`. Signed mass quadrature weights are face fluxes, not
positive physical material atoms. The existing17-component layout retains mass,
spatial first/second moments, momentum, mixed position/velocity moments,
velocity second moments and hydrodynamic energy. Potential flux is zero in
this experimental function; external gravity is not implemented here.

Two-point Gauss is exact for the stated affine traces through degree3. The
test reference uses independent three-point Gauss. On69,114inputs, CPU1,
OpenMP4 and CUDA agree bitwise with FMA disabled. Maximum scaled reference
error is8.41e-16, maximum absolute component error1.25e-14.
Evidence and hashes: `canonical_affine_moment_backend_2026-10-18.json`.

This supplies face transport only. In particular, spatial moments also require
volume sources (qdot includes p), and pressure/gravity supply additional
moment terms. Simply adding these face rates to all17state variables would
not be a complete physical evolution. No production RHS call site was added.

Remaining gates include consistent pressure/thermal closure, gravity work,
nonuniform reconstruction, moment realizability, integration and topology
changes, approved original-model reduction, MPI and historical long Kepler
validation. Neither positivity nor total evolved energy conservation follows
from component backend parity.
