# Common face transport implementation

`Exam/affine_face_transport.h` implements analytic integration of affine normal
mass flux, position and velocity along a straight 2D face. Six outputs are
mass, two spatial first moments, two momentum components and kinetic energy.
It is C-compatible and marked host/device for CUDA. No production call site
has been added.

On lageunha, `check_affine_face_backends.cu` evaluated 68,986 actual-geometry
inputs and 128 manufactured general affine inputs, including zero face length.
OpenMP runtime team sizes were checked as 1 and 4. CUDA allocation, execution,
synchronization and copies were checked. All three outputs agree bitwise with
FMA disabled. Independent two-point Gaussian quadrature gives maximum error
scaled by 1+|reference| of3.57e-16, maximum absolute component error1.78e-15.

Results and source/input/executable hashes are stored in
`canonical_affine_flux_backend_2026-10-18.json`. Remote artifact directory:
`/gpfs/kjhan/LagForce/Codex_affine_flux_backend_20261018`.

Compiler: gnu13/13.2.0 and cuda/13.0.2 modules, nvcc -O2 --fmad=false
-Xcompiler=-fopenmp,-ffp-contract=off. FMA-enabled production flags, MPI and
face accumulation reductions have not been checked by this test.

This is one reusable component, not an evolved solver. Thermal and pressure
fluxes, Riemann states, nonlinear density reconstruction, cell moments,
generator motion, weight evolution, boundary exchanges, positivity, checkpoint
semantics and the approved original Voronoi reduction remain separate work.
Do not infer long-time physical accuracy from this backend agreement.
