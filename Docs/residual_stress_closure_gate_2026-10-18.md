# Residual covariance cannot simply be a passive energy ledger

Read saved N4/N8 hybrid remap endpoint states. Compare residual covariance
stress scale R/V with thermal P=(gamma-1)U/V, using explicitly chosen gamma=5/3.
This is a diagnostic comparison, not a readback or change of production gamma.

For N8 inner r<.5, lambda_max(R)/[(gamma-1)U] has median .09563,
95th percentile .45602, maximum3.49604. Nine inner cells exceed one.
One additional edge cell exceeds one (edge maximum1.37141).
Energy-weighted residual anisotropy (lambda_max-lambda_min)/trace(R) is
.5395 inner and .8382 edge. N4 values are nearly identical.

These are volume-integrated tensor scales, not measured face forces. R/V
cannot be inserted into a face flux without a spatial reconstruction and
compatibility checks. Nor is this covariance automatically molecular thermal
pressure, turbulence, or an independently evolved physical gas component.

## Exact transport requirements

For a resolved underlying velocity u and moving face speed w, advective
momentum and kinetic-moment fluxes contain, respectively,

    integral_face rho*u_a*(u-w).n dA,
    integral_face rho*u_a*u_b*(u-w).n dA.

Replacing u by a mean/affine trace drops correlations unless the residual
field is accounted for. The second expression additionally requires third
velocity information, which is not determined by cell-integrated second
moments. Earlier donor-face diagnostics demonstrate that the omitted field
is not generally negligible. This remains true if its energy is stored in a
passive scalar without allowing it to affect transport.

The pressure-work terms must exchange kinetic and thermal energy consistently,
and subcell gravitational work must use the same state interpretation. Neither
an arbitrary relaxation of R into heat nor an extra isotropic pressure chosen
from trace(R) is justified by these conservation identities.

## Gates before physical integration

1. Define one represented subcell density/velocity/thermal field and its common
   face flux. Preserve m,q,Q,p,D,Pi with the positive realization.
2. State the approximation for residual spatial/third-velocity information;
   quantify its error, not just the conservation of its stored energy.
3. Verify pressure and gravity work with that same representation, including
   no spurious heat for the cold constant-pressure remap control.
4. Prove the approved persistently-equal-relative-weight fixed-mass reduction
   for full dynamics. Vanishing remap increments or equal final weights alone
   do not prove it. Do not erase historical mass or residual state by reset.
5. Only then integrate into the original RK4/CUDA/OpenMP/MPI solver and test
   through the historical t6--10 failure range using preserved checkpoints,
   followed by t20 and the longer orbital target.

These gates do not add user approval requirements; they are scientific
acceptance criteria under the existing authorization. They prevent a positive,
conservative remapper from being mistaken for a complete physical solver.
No production closure or new physical trajectory was adopted this turn.

Evidence: `canonical_residual_stress_scale_2026-10-18.json`.
Script: `audit_residual_stress_scale.py`.
Remote: existing hybrid directory, `residual_stress_scale.json`.
