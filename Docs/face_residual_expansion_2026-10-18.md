# Exact face residual expansion identifies which information is missing

On retained donor segments from the saved t2 conditional remap, write relative
donor velocity as u+d, where u is the recovered affine trace in the frozen
receiver-mean frame. Expand transport exactly in residual d. This uses the
known donor field; it is not a predictive closure from cell moments.

With un=u.n and dn=d.n, corrections to affine flux per unit density are:

    mass: dn
    momentum linear: u*dn+d*un
    momentum quadratic: d*dn
    kinetic linear: |u|^2*dn/2+(u.d)*un
    kinetic quadratic: (u.d)*dn+|d|^2*un/2
    kinetic cubic: |d|^2*dn/2

Integrate each term with the same positive donor density. All corrections
together recover donor flux within3.26e-19 absolute, about2.6e-16 relative L2.
Mass closes with linear residual; momentum closes with quadratic residual.

| Fraction / region | K relative L2 difference after linear residual | After quadratic residual |
|---|---:|---:|
| .25 inner | .921427 | .752122 |
| 1 inner | .796260 | .623272 |
| .25 edge | .999883 | .989269 |
| 1 edge | .997989 | .955988 |

Thus the cubic residual term is not small in this frozen-frame donor trace.
Terms can correlate and cancel: these norms are not additive causal shares.
Also d is a spatial residual relative to the affine trace, not a random
fluctuation with zero mean at each face. Its cubic term must not be equated
automatically with an independent central third-moment field.

## Design consequence

Correcting only face mean/mass transport will not reproduce momentum or
kinetic transport in these cases. Supplying only cell-integrated covariance
also cannot determine the spatially varying d on each face. A useful next
comparison is a compact local face representation of donor velocity segments
against storing additional volume moments. Avoid assuming that ever more
global moments are the cheapest or best way to preserve face information.

This is an exact algebra/diagnostic decomposition of the same conditional
field, not evidence that exact donor trace matching is required for every
convergent numerical scheme. It does not supply a common Riemann flux, full
closed-cell balance, pressure work, moving physical face speed, original-model
reduction or a Kepler cure. Original production runs remain untouched.

Evidence: `canonical_face_residual_expansion_2026-10-18.json`.
Script: `audit_face_residual_expansion.py`.
Remote: checkpoint remap directory, `face_residual_expansion.json`.
