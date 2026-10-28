# Restoring transported spatial moments resolves first-remap covariance failure

Reintegrated the same uniform donor fields and first-remap overlap geometry
as the failed affine-path diagnostic, now retaining all 17 moment slots.
Fractions .25,.125,.0625,.03125 correspond to its first transitions. Slot16
potential energy is unused/zero in this diagnostic; no potential claim.

The shared reduced m,p,D,Pi,U states agree with the failed audit within
1.518e-18 maximum absolute error. Replacing the imposed uniform spatial moments
by the actual transported q,Q yields zero residual-covariance failures in all
four cases, versus 238,227,212,194 failures previously. Minimum residual
eigenvalue is -8.96e-19 (roundoff; not strict positivity); minimum S eigenvalue
exceeds 6.04e-10. Thermal energy stays positive, minimum 3.70952e-8.

Thus the first-remap failure was caused by the inconsistent spatial projection
in the tested candidate, not by an unavoidable failure of affine recovery.
No clipping, added energy, changed donor velocities or reduced thermal energy
was used to remove the failures.

## Uniform boost check

Transform p,D,Pi analytically under a uniform velocity boost b=(7,-3), using
the transported q. Recovered central C changes by at most 1.735e-18, G by
2.102e-11, R by 3.146e-17 absolute. Small S and raw-moment subtraction amplify
floating-point cancellation in G; this is not bitwise invariance. The test
confirms numerical consistency of the specified moment transformation at
these scales, not a full moving-frame simulation.

## Remaining gate

The transported q,Q have an explicit positive realization: the original
piecewise donor densities on overlap fragments. A compact reconstruction used
for the next remap must preserve these moments on the receiver polygon. The
present calculation does not provide that reconstruction or complete a full
weight path. PSD of the joint moments alone is not sufficient to justify an
arbitrary spatial density reconstruction. Retaining a realizable fragment
representation is a diagnostic reference; compression must be evaluated for
moment preservation and cost before production adoption.

Next compare a positive moment-preserving density representation against
the known overlap-fragment realization, before testing repeated affine remaps.
Keep U and residual kinetic covariance separate, and verify angular momentum
from D. No original point-model reduction, pressure/gravity evolution, MPI,
backend parity or long-time Kepler stability is established by this audit.

## Long original runs

At 2026-10-17 18:02:25 UTC the one-shot process checker confirmed live processes
and finite latest accepted states: baseline t15.611318, fixed-HG t15.464273,
from-t0 timestep-cap t11.507506. They were not restarted or modified.

Evidence: `canonical_remap_true_spatial_2026-10-18.json`.
Script: `audit_remap_true_spatial.py`.
Remote: `/gpfs/kjhan/LagForce/Codex_remap_true_spatial_20261018`.
