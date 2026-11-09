# Physical closure of evolving Laguerre weights

## Scope and history

The user asks for weights to represent material deformation, with fixed
mass, convex cells, no steering, few parameters and practical runtime.
Do not interpret the failure of one auxiliary-inertia model as rejecting all
dynamic weights. The successful Kepler W0 model remains the control.

The previously retired psi-times-fitted-strain construction remains retired:
its equilibrium acoustic frequency ratio lower bound is about36, it worsens
forces and its coefficient solve is too costly. A new volume-driven proposal
must not silently reuse that construction. The recent canonical inertia
candidate is also not selected: Gresho radial drift and KH geometry stop near
1.07313 remain, and the extra inertia lacks a physical velocity-field derivation.

## New modal audit: correct frequency is insufficient

This is an algebraic 1D reduction of the current 2D inertia candidate, NOT a
rerun or diagnosis of the successful 1D material Riemann/PLM solver.
Assume the previously documented dimensional scaling mu=rho*h^(d-2).
For a uniform state and Fourier angle theta, use mass-normalized physical
velocity and auxiliary weight velocity. After removing a harmless Fourier
phase, the pressure-coupling vector (in units c/h) is

    a = (sin(theta), 1-cos(theta)).

The positive kinetic-coordinate acoustic matrix is a a^T. Its nonzero
eigenvalue is 4 sin(theta/2)^2, but its eigenvector has auxiliary kinetic
energy fraction sin(theta/2)^2. This is a fraction of modal kinetic energy,
NOT of kinetic-plus-compressive total energy. For wavelengths16,8,4,2 cells,
the fractions are3.806%,14.645%,50%,100%. At the Nyquist limit the acoustic
mode is carried entirely by the auxiliary coordinate rather than particle
velocity. A second neutral eigenvector is orthogonal to a: physical and
weight contributions cancel in the linear volume rate. It can have nonzero
physical velocity despite zero pressure response. This does not prove that
the neutral mode is excited in KH or causes the geometry failure.

The low-wavenumber auxiliary fraction tends to zero quadratically; therefore
this result alone is not a proof of nonconvergence. It does show why matching
sound speed and conserving the augmented energy cannot establish that the
resolved variables represent the intended Euler fluid. The extra coordinate
must be constrained or associated with a derived physical internal motion.
Artifact: Research/MaterialVolume_20261105/weight_physical_modes_20261106_v1;
audit.py numerically verifies both eigenspaces and analytic fractions.

## Required power identity for a kinematic alternative

For fixed masses let actual geometric volume rate be

    Vdot = B v + L Wdot,
    Udot_i = -P_i Vdot_i  (smooth isentropic flow).

The native translational pressure force is B^T P. If its pressure work is
left unchanged, total fluid energy changes by -P^T L Wdot. Simply replacing
the independent weight acceleration with a convenient geometric Wdot thus
breaks conservation; adding that deficit as heat is not automatically physical.

If a new closure is Wdot=R v (linear in velocity at fixed state), define
C=B+LR. Then pressure force C^T P and Udot=-P C v form the matching power
pair. This algebra is dimension-independent. It is a necessary smooth-flow
identity, not a complete shock solver or a guarantee that all faces are
material boundaries. A pressure-dependent closure requires a separate power
and entropy derivation; do not apply the velocity-only formula to it.

## Design decision

Prefer a geometric constraint derived from the represented fluid velocity
and material deformation over an independently assigned weight inertia.
But do not implement another R solely because it reproduces affine strain:
require bounded acoustic norms, cheap construction, no hidden compressive
neutral motion, equilibrium tests, and preservation of the successful 1D
material pressure/work behavior. No extra damping coefficient or artificial
radial force is justified by the current evidence.

If independent weight inertia is retained, derive its kinetic metric from
an explicitly represented cell-internal velocity field (including any needed
cross terms); the standalone pi^2/(2mu) term is not such a derivation.
This is a research requirement, not authorization to add a large subcell model
or an expensive per-step volume integration. Offline quadrature diagnostics
must not become the production method by default.

Remaining implementation is unresolved. Continue with native C/CUDA geometry,
keep previous results, finish 2D before3D, and verify end-to-end runtime against
GIZMO. This note does not claim the objective is complete.
