# Native periodic weighted K1 qualification — 2026-11-06

This is not an accepted time evolution or KH/Gresho accuracy pass. The 1D
material Riemann/PLM model is not claimed equivalent to this native variational
pressure/weight-inertia model. 2D remains the priority; no 3D work was launched.

## Actual execution

Syntax A10, one MPI process, two OpenMP threads. Shared /gpfs files are accessed
directly from LagEunha. No /scratch use. Artifact directory:
Research/MaterialVolume_20261105/native_periodic_rhs_20261106_v1

The original native reader roundtrip passed for frozen and dynamic KH/Gresho.
The new explicit GFS_WEIGHT_PERIODIC_RHS_ONLY=1 gate requires both existing
K1-only flags and the periodic poststage callback, and rejects wall events.
It does not enable periodic RK4 evolution or the different unified force model.
The production RK4 wall gate is unchanged.

C and CUDA now permit periodic images only for this explicit diagnostic path.
The immutable geometry-weight snapshot comparison remains enforced.
Existing CUDA archive lacked sm_86; rebuilt exam_gpu and gfs_unified_gpu for
A10. Rebuilt C/extraction objects and relinked in an isolated artifact directory.

## Reproduced defect and correction

Tessellation writes double-precision volume only for owned cells. The periodic
field refresh omitted weight_volume64. The weighted force's compression bound
reads neighbor V64, including periodic aliases. In v4 all owned volumes were
positive, but 60 of 256 Gresho cells and 92 of 512 KH cells had nonfinite bounds.

The existing device alias refresh now copies V64 and checks matching W, Wdot,
inertia and U64. Native periodic K1 requires that refresh to succeed; the
existing direct ID-to-index assumption fails closed rather than being bypassed.
General reordered-ID and multi-rank mapping remains to be implemented/qualified.

## v5 measured outcome

Both K1 probes exit zero and explicitly stop before any accepted evolution.
All exported cell fields and compression bounds are finite.
Gresho 16x16: volume sum 1; sum weight force 3.89e-17; augmented instantaneous
energy-rate sum -1.10e-17.
KH 16x32: volume sum 2; sum weight force 0; augmented instantaneous energy-rate
sum -1.13e-18.
Initial W and pi are zero, so these tests do not validate finite evolving W or
weight energy exchange. They do exercise active-weight kernels and their CFL.

## Remaining work

Nonzero weight/rate periodic states and reordering; energy/work and geometry
finite-difference validation; bounded periodic stage acceptance/rollback and
actual 2D evolution. The current rectangular alias checker uses one box length;
a general rectangular-domain check must use separate x/y lengths.
KH/RT endpoints, shock consistency, MPI parity, and the GIZMO <=2x runtime
target remain unverified. Do not call this a completed 2D implementation.

## Nonzero weight/rate extension

Artifact: Research/MaterialVolume_20261105/native_weight_exchange_20261106_v1.
These are manufactured column-weight fixtures using the square/rectangular
IC layouts, not Gresho/KH solutions: velocities reset to zero, W=.01*h^2*cos(2pi*x),
Wdot=sign*.01*h^2*sin(2pi*x), P=10+.2*sin(2pi*x); masses retained.
Exact weighted-column volumes initialize thermal energy. All six signed/zero
rate K1 probes complete on the native A10 weighted kernel.
Thermal power plus sum(Wdot*Q) cancels to <=2.72e-19. Reversing Wdot reverses
both rates; zero Wdot produces zero heat. Both grids have positive volumes
and volume sums agreeing with domain area to roundoff.
Analytic per-cell volume and -P*dV/dt comparisons are saved separately.

This verifies the instantaneous energy ledger includes weight kinetic energy.
It does not justify a new physical inertia model or prove the model matches the
1D material Riemann method. No accepted time evolution or KH/RT pass is claimed.
