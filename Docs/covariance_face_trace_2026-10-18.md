# Constant specific covariance with zero skew is insufficient in the saved face audit

Reuse the actual conditional t2 remap states and retained face/donor segments
from `remap_face_trace_2026-10-18.md`. Add a candidate with affine local mean u,
constant specific residual covariance Sigma=R/m and zero third central moment.
All velocities below are relative to the frozen receiver-mean frame; this is
not the evolving physical Laguerre face speed.

At a quadrature point, the candidate moments of the flux are

    mass = rho (u.n),
    momentum = rho [u(u.n)+Sigma n],
    kinetic = rho/2 [(|u|^2+tr Sigma)(u.n)+2u.Sigma n].

The missing kinetic term for a general velocity distribution is one half of
rho times the normal component of its central third contracted moment.
No gas-kinetic/Riemann model or physical thermal pressure is supplied here.

| Weight fraction | Inner relative K-flux L2 difference, affine only | With covariance |
|---|---:|---:|
| .25 | .978845 | .944414 |
| 1 | .882738 | .831494 |

Relative differences in edge kinetic flux remain .999968 and .999414.
Inner momentum-component differences improve from .945/.959 to .916/.927
at fraction .25, and from .811/.859 to .777/.816 at fraction1.
Mass-flux differences are unchanged (.850/.713 inner), because the covariance
term cannot change the local mean velocity. Thus missing third moments are
not the sole issue: the spatial mean trace itself differs from the donor field.
Order8/12 quadrature differences stay below6.64e-15 absolute.

## Decision and limits

This covariance closure modestly improves selected flux comparisons but does
not recover most of the missing donor trace in this diagnostic. Do not add it
to production as a completed physical closure. Do not infer that every
zero-third-moment approximation is inconsistent from this one trace test:
exact donor discontinuity preservation is not a universal convergence criterion.

The comparison uses the same donor density, one-sided receiver states, and
excludes ambiguous coincident donor boundaries. It is neither a common
conservative interface flux nor a complete closed-cell balance. Large relative
differences are not global energy losses or measured Kepler force errors.

Next separate spatial residual mean at faces from covariance/third-moment
effects before extending the state with more moments. The unchanged mass-flux
difference shows why simply adding a higher velocity moment cannot by itself
repair this candidate's face reconstruction. Preserve original production
trajectories and approved reduction conditions.

Evidence: `canonical_covariance_face_trace_2026-10-18.json`.
Script: `audit_covariance_face_trace.py`.
Remote: original checkpoint-remap directory, `covariance_face_trace_audit.json`.
