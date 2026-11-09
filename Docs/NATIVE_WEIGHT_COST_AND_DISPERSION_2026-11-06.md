# Native weight cost and uniform-grid dispersion

## Measured native K1 comparison

Artifacts: Research/MaterialVolume_20261105/native_weight_cost_20261106_v1.
Same IC generator, same native A10 executable, one MPI rank, two OpenMP threads.
Only the frozen/dynamic feature differs. KH initial conditions on 16x32 and
32x64 periodic grids. All four K1 probes terminate normally without evolution.

| Nx | frozen candidate dt (printed) | dynamic candidate dt (printed) |
|---|---:|---:|
|16|0.006123|0.0007813|
|32|0.003062|0.0003906|

The extra weight bound reduces the initial candidate step by approximately
7.84 at both resolutions. This is not a measured trajectory runtime ratio and
not a GIZMO comparison. The dynamic implementation uses dt <= .25/omega_bound;
omega_bound is 320 and 640 respectively. Do not loosen this safety factor
merely to pass the speed target. No long dynamic run was launched.

## Independent linear uniform-grid analysis

For constant density, pressure, cell side h, and dimension d:
m=rho*h^d, volume=h^d. The dimensionally consistent inertia scaling matching
the existing 2D choice mu=rho is mu=rho*h^(d-2). This scaling is an analysis
assumption, not a production change or a derivation of physical subcell inertia.

The volume operator is Vdot=B v + L Wdot. B is the centered volume derivative;
L is half the graph Laplacian times h^(d-2). With the native adjoint forces and
Wdot=pi/mu, the positive linear acoustic block is

K = (gamma P / h^d) [B M^-1 B^T + L mu^-1 L^T].

For Cartesian Fourier angles theta_a its eigenvalue is

omega^2 = (c^2/h^2) [
  sum_a sin(theta_a)^2 + 4 (sum_a sin(theta_a/2)^2)^2 ].

In 1D this reduces exactly to 4(c/h)^2 sin(theta/2)^2.
The weighted term removes the centered-pressure checkerboard null.
In 2D the maximum is omega=4c/h; the low-wavenumber sound speed is c.
Matrix eigenvalues agree with the formula to 7.5e-16 relative for the saved
1D/2D small algebraic checks. These are not PDE runs; no 3D test was launched.

Thus a regular-grid dispersion calculation does not show the ~36x frequency
explosion of the already retired strain-correction model. It also does not
establish stable steps on deformed meshes, correct shear/entropy behavior, or
an equivalence with the 1D material Riemann/PLM solver.

## Decision and next requirement

The current bounded dynamic RK4 setting is not selected as a production
solution: its measured initial timestep penalty is substantial, and the
physical interpretation of auxiliary weight kinetic energy remains a separate
requirement. Preserve the force/work and periodic-field fixes independently.

Before production adoption, a common pressure/work/entropy model must agree
with the qualified 1D results and fit the runtime target. A KDK formulation
must include the same position and weight canonical pair and geometric
positivity checks, not simply switch the time-integrator flag or omit the
weight energy. Its admissible step must be established from actual evolution,
not inferred from the regular-grid eigenvalues alone.

2D KH/RT and shock suites, full timing against GIZMO, and subsequent 3D remain
incomplete. No claim of algorithm completion or candidate retirement of all
weighted methods is made.
