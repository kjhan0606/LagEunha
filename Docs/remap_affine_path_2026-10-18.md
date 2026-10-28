# Uniform-density affine projection fails realizability after first remap

Compared the same 0 -> .25 -> 0 weight path at 1,2,4,8 subdivisions per leg.
The candidate retains m,p,D,Pi,U and reconstructs affine velocity plus explicit
nonthermal residual covariance. It reimposes uniform density in each receiver,
therefore replacing the actual mass centroid and spatial covariance with the
geometric uniform-density values. No physical evolution or production change.

Every case fails at recovery after its first remap; none completes the path.
Negative residual eigenvalues are not clipped and no energy is added.

| Subdivisions | First fraction | Failing cells /1590 | Most negative residual eigenvalue |
|---|---:|---:|---:|
| 1 | .25 | 238 | -4.90584e-5 |
| 2 | .125 | 227 | -1.77826e-5 |
| 4 | .0625 | 212 | -5.70364e-6 |
| 8 | .03125 | 194 | -1.66043e-6 |

For the worst cell at fraction .25 (ID691), the imposed affine kinetic energy
is 2.45360e-5, while available internal kinetic energy is 7.68147e-9. Tolerance
for that residual eigenvalue is 4.00234e-17. This is not roundoff. Aggregate
m,p,D,Pi,U drift in each first remap is at most 7.11e-15, demonstrating that
aggregate conservation does not ensure a realizable reconstructed state.

## Inconsistent spatial mass projection

Let q_x be the true transported spatial first mass moment and c_true=q_x/m.
Replacing it by m*c_geom while holding D and p fixed yields

    C_imposed = D - m*c_geom*vbar^T
              = C_true + (q_x-m*c_geom)*vbar^T.

The extra term is not a genuine central position-velocity correlation. Under
a uniform velocity boost b, true D transforms as D+q_x*b^T and p as p+m*b;
the imposed C acquires (q_x-m*c_geom)*b^T. A correct central covariance is
boost invariant. Resetting spatial mass moments while retaining D thus breaks
the state interpretation. Changing the imposed spatial covariance S compounds
the problem when recovering G=C^T S^-1 and R=T-G S G^T.

This algebra identifies a structural inconsistency in this tested projection;
it does not quantify every contribution to each failing cell. The full true
transported spatial moments were not saved in this reduced-state diagnostic.

## Decision

Reject this uniform-density reset. Do not reject affine recovery as a class or
alter production checkpoints. A next reconstruction must retain consistent
spatial first/second mass moments together with D and Pi, and realize them
with positive density on the receiver polygon. The existing 17-moment state
already has these fields; adding another conserved scalar is not the remedy.
Recheck the first-remap failure using actual transported spatial moments
before expanding any closure or running a new physical trajectory.

Evidence: `canonical_remap_affine_path_2026-10-18.json`.
Script: `audit_remap_affine_path.py`.
Remote: `/gpfs/kjhan/LagForce/Codex_remap_affine_path_20261018`.
No backend parity, MPI, coupled time evolution or long-time stability claim.
