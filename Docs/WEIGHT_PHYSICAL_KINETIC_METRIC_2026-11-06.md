# Weight rates and a physically represented internal velocity field

Reused existing Exam/gfs_moment.h, specifically gfs_moment_accumulate and
gfs_moment_primitive. This is the existing affine moment representation,
not a new spin implementation. Previous archived affine RK4 failures remain
as documented in ROTATION_IMPLEMENTATION_RECONCILIATION_2026-11-05.md.
No old branch was enabled or overwritten.

## A representable alternating-column weight mode

Take rectangular cells hx by hy, density rho, mass m=rho hx hy. Generators
are initially stationary, W=0, and neighboring columns have alternating
weight rates +q and -q. The left/right x-face velocities of a +q cell are
-q/hx and +q/hx, while y-face normal velocities are zero. A compatible affine
physical field is ux=2q*x/hx^2, uy=0, with zero mean velocity relative to the
cell center. Neighboring cells share the same interface normal velocity.

Its integrated kinetic energy is

    K_internal = m q^2/(6 hx^2).

The current auxiliary model uses .5*rho*q^2 in2D. For square cells this is
three times the affine-field kinetic energy. For hy/hx=.5,1,2 the ratios are
6,3,1.5. These are comparisons to this explicitly specified velocity field,
not universal correction factors or proof that every mass-lumped metric is
invalid. Do not divide the current inertia by3 and call the physical problem
solved: other modes and geometries have different mappings and cross terms.

The unchanged moment helper recovers the analytic energy for all three
aspect ratios, and retains the supplied thermal energy to <5.4e-24 absolute.
This is a state decomposition check, not time evolution or a stability result.
Two-point tensor quadrature is exact for the quadratic integrands here and
is used only in this small diagnostic, not proposed as a production quadrature.

## General1D boundary traces expose further requirements

For a cell of width h, prescribe mean physical velocity v and boundary
velocities wL,wR. A single affine profile can satisfy all three only if
v=(wL+wR)/2. If that fails, adding an independent affine moment alone cannot
repair the geometric mismatch.

A minimal polynomial illustrating the missing information is
u(x)=v+a*x+b*(x^2-h^2/12),
a=(wR-wL)/h, b=3*(wR+wL-2v)/h^2.
It has the required mean and traces, and exact kinetic energy

    K = .5*m*v^2 + m*(wR-wL)^2/24 + m*(wR+wL-2v)^2/40.

The last term is a non-affine internal contribution, not heat. The identity
was independently checked in1000 random cases by4-point Gaussian integration,
maximum relative error1.45e-15. This polynomial is a diagnostic representation,
not a new evolved state or an adopted closure. A kinetic reconstruction alone
does not supply pressure/entropy moments, full dynamics or a cheap solver.

## Implication for common dimensional design

The physical metric must follow from the chosen fluid-velocity representation
and the actual geometric boundary map. It generally couples bulk velocity,
internal motion and weight rates; a standalone diagonal weight inertia lacks
that derivation. Existing affine moments are a useful component but are not
already connected to a compatible Laguerre weight law. Higher-order internal
motion must not be silently counted as thermal energy, nor introduced without
a cost and closure justification. No new production physics/state was added.

This also qualifies the earlier modal audit: auxiliary-mode energy is not
necessarily unphysical IF it represents actual internal fluid motion. The
current code has not established that mapping. Correct acoustic frequency
alone cannot substitute for it.

Artifact: Research/MaterialVolume_20261105/weight_physical_metric_20261106_v1.
audit.c/results.csv reuse the C header; quadratic_identity.py/.json retain the
independent polynomial check. No PDE,3D or GIZMO performance claim follows.
