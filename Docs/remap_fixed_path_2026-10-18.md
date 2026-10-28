# Fixed weight path: projection mixing and represented angular momentum

Use the same t2 uniform spatial lift and piecewise-constant overlap projection
as the cycle diagnostic, now on one common path 0 -> .25 -> 0, with one
reversal. Subdivide each leg into 1,2,4,8 increments. No physical time advances.

Alongside m,p,Pi,U, transport an auxiliary angular-momentum scalar ledger
initialized from the uniform spatial field. At the final original geometry,
independently calculate the angular momentum actually represented by the
reconstructed field: sum(c_x p_y-c_y p_x), where c is polygon centroid.
The ledger is deliberately not claimed to be a spatially consistent moment.

For 1,2,4 subdivisions per leg, internal K / original U is respectively
.0425559, .0426403, .0426824. Represented relative angular-momentum change is
respectively -1.65773e-5, -1.65848e-5, -1.65886e-5. The scalar ledger changes
only by 1.30e-14, 2.60e-14, 5.19e-14 absolute. Thus preserving a scalar sum
does not prove that reconstructed momentum carries the same angular momentum.
At eight subdivisions, internal K / original U is .0427035 and represented
relative L change is -1.65905e-5; scalar ledger drift is 1.034e-13 absolute.
All four cases completed successfully. Coverage error stays below 6.60e-14.
Complete results: `canonical_remap_fixed_path_2026-10-18.json`.

## Implication for reconstruction

For mass centroid c, mean velocity v and central cross moment
C=integral rho (x-c)(u-v)^T, actual angular momentum is

    Lz = m(c_x v_y-c_y v_x) + C_xy-C_yx.

A piecewise-constant velocity reconstruction sets C=0. It cannot generally
represent an independently transported angular-momentum ledger. The existing
experimental affine recovery G=C^T S^-1 represents this cross moment if the
same S and density realization are used. This is a reason to test moment-aware
reconstruction, not to add an arbitrary angular-momentum correction or steering.
It still leaves residual velocity covariance, spatial realizability and
consistent face transport to establish.

## Limits and next gate

Fixed spatial resolution and a single weight path: subdivision tests numerical
composition of remapping, not convergence of the physical RK4 evolution. The
uniform lift does not preserve the original point Lz/Phi by construction;
the baseline for this audit is the lifted field's own angular momentum.
Pressure, gravity, shared physical fluxes, MPI and long Kepler stability are
not tested. Internal kinetic covariance is retained, not thermalized.

Next compare a reconstruction constrained by C as well as p on this same
path, rejecting any state that cannot realize S,C,T with positive density
and thermal energy. A conserved angular-momentum ledger alone is insufficient.
Do not change original solver trajectories or adopt this projection in them.

Script: `audit_remap_fixed_path.py`.
Remote: `/gpfs/kjhan/LagForce/Codex_remap_fixed_path_20261018`.
