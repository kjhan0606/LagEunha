# Static condensation screen of matched internal acoustic representation

Continuation of WEIGHT_METRIC_PRESSURE_CONSISTENCY_2026-11-06.md. Uniform1D
linear rest-state analysis only; no new hydrodynamic trajectory or3D test.
Let X be mean displacement and W the geometric weight perturbation. Units
h=rho=c=1. The matched kinetic matrix M and compression matrix A are retained.
Instead of evolving the fast internal coordinate independently, impose the
stationarity condition partial_W compression_energy=0 at each X:

    W = -A_WW^-1 A_WX X = i sin(theta)/(2+cos(theta)) X.

This is static condensation, not adding damping or dropping energy after an
accepted step. Lift both kinetic and potential energy onto z=(X,R X):
Mred=[1,R]^H M [1,R], Ared=[1,R]^H A [1,R]. The reduced velocity obeys the
same lift for this constant-coefficient linear problem. Retaining Mred is
necessary; retaining only Ared would be a different approximation.

## A cheap, well-conditioned relation in this special geometry

Restoring h, the real-space relation is

    2 W_i + (W_{i-1}+W_{i+1})/2 = h (X_{i+1}-X_{i-1})/2.

Its Fourier symbol is 2+cos(theta), so the coefficient condition number is
at most3, independent of grid size. It does not require the unbounded inverse
Laplacian appearing before symbolic cancellation. CG at relative tolerance
1e-12 takes21 iterations for random data at64,256,1024,4096 cells. Diagnostic
Python solve time is .00055--.00090 seconds; this is NOT a native PDE timing.
A cyclic tridiagonal solver or native implementation has not been built.

## Acoustics

With c=cos(theta),

    omega^2 = 20 (2-c-c^2)/(c^2+13c+16).

The derivative with respect to c is negative on[-1,1], so max omega=sqrt(10)
in normalized units. The full uncondensed matched model had high internal
frequency approaching sqrt(60); its removed branch is not time-integrated.
The original arbitrary-diagonal-inertia model had maximum2, while the
centered frozen-W0 volume model has maximum1 and a checkerboard null; do not
confuse these references or translate spectral ratios into GIZMO runtime.

Reduced low-wave speed ratio is1.000000000001 at1024-cell wavelength and
1.00000006464 at64-cell wavelength. At4-cell wavelength it is1.00658424 versus
1.00375412 for the full model's low branch. Static condensation is therefore
an approximation to the slow branch, not exact elimination of finite-frequency
dynamics. Both effective mass and stiffness are positive for nonzero modes.

## Adoption limits and next gate

This supplies a parameter-free, neighbor-coupled linear closure derived from
matched energies; it is distinct from the retired psi-times-strain model and
from pressure relaxation with geometric-mismatch heating. It does not require
all faces to move at the same speed.

Do not impose this relation on existing nonlinear checkpoints: it would
change their state/energy unless a consistent reduced-model initialization
and thermodynamic map are supplied. Fast internal transients after heating
or shocks are not captured by the reduction. Nonuniform entropy, large
material deformation, actual Laguerre2D geometry, weight admissibility,
compatible nonlinear forces, and shock entropy remain unresolved. The
constant-condition bound is proved ONLY for the uniform1D coefficient
operator; it must not be assumed for arbitrary2D/3D cells.

Next useful work is to derive the corresponding geometric/thermodynamic
closure in2D from existing moments and reject it if complexity or conditioning
loses the simplicity/runtime advantage. No production code or physical
parameter was changed. Full2D validation and <=2x GIZMO runtime remain open.

Artifact: Research/MaterialVolume_20261105/weight_static_condensation_20261106_v1;
audit.py/results.json retain generalized-eigenvalue and coefficient-solve
checks. The actual1D Riemann/PLM reference and successful Kepler code remain.
