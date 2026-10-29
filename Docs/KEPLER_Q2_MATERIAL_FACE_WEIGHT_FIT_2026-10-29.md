# Can scalar weight rates make the Q2 faces material?

At the saved W0 Q2 t30 and t47.5 positions, give every generator and every
face point the analytic softened Kepler velocity field. For each face, the
difference between material normal velocity and W0 geometric face velocity
is evaluated at two Gauss points. Power-weight rates change the normal face
speed by the **constant along that straight face**

`delta_w_n = (Wdot_i - Wdot_j)/(2*distance_ij)`.

A length-weighted least-squares problem fits one unconstrained signed Wdot
per generator, including all neighbors of the 3,931 inner cells. It includes
12,508 faces and 4,646 independent generator variables. No rate bounds,
positive-weight constraint, pressure work, or later weight evolution are
imposed. Even allowing an independent constant per face provides a lower
bound on the pointwise residual. At each frozen W=0 snapshot, a nonzero
Wdot would immediately leave W=0; this is not an evolution that keeps the
approved fixed-mass W0 limit.

| Diagnostic | t=30 | t≈47.5 |
| --- | ---: | ---: |
| Pointwise speed-mismatch RMS remaining after weight fit | 99.989% | 99.970% |
| Pointwise lower bound, independent constant per face | 99.969% | 99.884% |
| Face-mean mismatch RMS remaining | 80.81% | 86.06% |
| Cell-net area-rate mismatch RMS remaining | 32.70% | 21.94% |

The denominators are the corresponding uncorrected norms. Face means
represent only 0.063% and 0.232% of the total length-weighted squared
pointwise mismatch in the two samples: most of that norm comes from the
variation *along* a face, which a scalar Wdot cannot remove. However, the
opposite signed variations largely cancel when integrated over faces and
cells. Therefore the 99.9% pointwise residual must **not** be interpreted
as a 99.9% residual density or cell-volume error. In these snapshots the
same fit reduces the cell-net mismatch substantially, even though it was
not optimized directly for cell-net rates. The norm there is an unweighted
cell RMS of net area-rate mismatch divided by area; it is not a finite-time
density result.

As controls, rigid-body rotation gives maximum face slip 2.22e-16, adding
a common Wdot gives zero change, and a random sparse-Jacobian probe agrees
to 5.33e-15 absolute. LSMR reports convergence with normal-equation residual
below 2.2e-7, but this is numerical evidence rather than a singular-value
proof. Input hashes, all norms and the exact script hash are in the
[receipt](KEPLER_Q2_MATERIAL_FACE_WEIGHT_FIT_2026-10-29.json); the
[reader](audit_kepler_q2_material_face_weight_fit_20261029.py) uses existing
snapshots only.

This constrains a purely material-face fixed-mass model: scalar weights alone
cannot make the sampled shear faces material at both quadrature points in
these states. It does not disprove weight changes for density control,
nonuniform subcell states, the user's weight-linked mass relaxation, or
spin. It also does not establish that enforcing local materiality is required
for a different conservative discrete solver. A practical candidate still
needs a coupled pressure/energy/momentum/angular-momentum law and native
finite-time and speed validation.
