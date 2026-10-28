# Kepler material-face weight-rate limit at accepted native t=6.5

CPU-only Slurm jobs **437358** (5% site jitter) and **437359** (aligned)
completed with exit `0:0` in 8 seconds each. The
[auditor](../audit_kepler_bounded_material_rate_t65_20261027.py) rehashes the
17-check accepted native t=6.5 endpoint and compares two weight-rate choices
on its actual evolved generator positions. Both use the already solved static
Laguerre weights that restore each disk cell's initial area:

1. **Volume rate:** solve `J Wdot = -Vdot_fixed`, so the signed sum of face
   sweeps in every cell approaches zero.
2. **Material-face least squares:** minimize the sum of
   `L (u_face,n - v_K(c)·n)^2` over the 145,671 disk-disk faces, using one
   scalar `Wdot` per disk cell. Here `v_K(c)` is the analytic softened Kepler
   velocity evaluated at the face centroid. This is the unconstrained optimum
   for this precise face-speed objective and fixed generator velocities.

| Quantity at t=6.5 | Aligned | 5% site jitter |
|---|---:|---:|
| Fixed-weight face-slip RMS | `1.63399e-4` | `2.02582e-4` |
| Volume-rate face-slip RMS | `1.435995e-4` | `1.811843e-4` |
| Material least-squares face-slip RMS | `1.414590e-4` | `1.791143e-4` |
| Extra RMS reduction, volume rate to least squares | `1.49%` | `1.14%` |
| Least-squares residual relative to fixed-weight RMS | `86.57%` | `88.42%` |
| `||Vdot||₂/||Vdot_fixed||₂`, volume rate | `9.12e-11` | `9.81e-11` |
| `||Vdot||₂/||Vdot_fixed||₂`, material least squares | `0.42886` | `0.29576` |
| Inner `0.5≤r₀<0.75` p99 area change after one Euler step, volume rate | `3.014e-5` | `3.117e-5` |
| Same, material least squares | `7.272e-5` | `8.408e-5` |

The least-squares solve took 1,404 and 1,467 cold Jacobi-preconditioned CG
iterations. Thus a more exact scalar-weight face-matching solve gives only a
small improvement over the volume-rate solution on these snapshots while
giving up substantial cell-volume control. The remaining face-slip RMS is
the least-squares projection residual for this graph and analytic velocity;
it cannot be removed by choosing different scalar `Wdot` values alone while
the generator velocities are fixed. This does **not** rule out a different
fluid reconstruction, a conservative flux across moving faces, or a
time-dependent effect at another Kepler state.

The [aligned](aligned_material.json) and [jitter](jitter_material.json)
receipts include source checkpoint and input hashes, the three slip measures,
cell-volume rates, and the one-step power-diagram geometry. The saved
[aligned](aligned_material.npz) and [jitter](jitter_material.npz) arrays contain
the solved rates. SHA-256 values, in that order, are
`94db552d1932c8f60eedffd1b817d8c1b7649057c7572faba1c4a5aaf350bc47`,
`8e770ad521c994800ccec742bb9855724440cc3aff7422ddbed0fe0d2c6b40a7`,
`61117871cccfaed234c0c26bddb9370507cd6abf090bb2bcaa066e1efa254c16`,
and `a270bd39eb834b5a06455408cb9569e82427d5365d033afdca17985604cec2d5`.

**Implementation decision:** retain the inexpensive truncated volume-rate
solve as the geometry candidate. Before enabling evolving weights in native
hydro, derive and audit a conservative treatment of the residual face
transport and of pressure work caused by changing weights. Preserve the
approved exact `W=0` fixed-mass Voronoi branch and unsteered generators.
This experiment is geometry-only: it has no native weighted hydro step,
CUDA timing, cell-spin state, or long-time stability result.

## Pressure work and transport scale at the same frozen state

The [pressure-work auditor](../audit_kepler_bounded_weight_work_t65_20261027.py)
evaluates the exact face/cell algebraic identity
`Σ_i P_i Vdot_i|W = Σ_(i,j) [L_ij/(2d_ij)](P_i-P_j)(Wdot_i-Wdot_j)`.
Its face and cell sums agree to floating-point precision, and all-cell
weight-induced volume rates sum to zero. With the **native** t=6.5 pressure,
the signed term is `-1.240e-13` aligned and `+6.310e-12` jittered per unit
physical time. The sums of absolute face terms are `2.396e-11` and
`3.714e-10`. Recomputing active-disk pressure adiabatically for the solved
static weighted areas gives signed proxy terms `-6.368e-13` and
`-6.218e-13`; inactive-cell pressure is left at its native value. These
numbers are small for this cold, near-uniform-pressure state, but a weighted
hydro operator still needs the weight-work term for general conservation.
The [aligned](aligned_weight_work.json) and [jitter](jitter_weight_work.json)
receipts include the input hashes and normalizations.

The [transport proxy](../audit_kepler_bounded_material_flux_t65_20261027.py)
uses arithmetic face density and the analytic Kepler velocity. It reports
the **gross absolute** disk-disk mass crossing that would arise from
`rho_face L (v_K,n-u_face,n)` if a moving-face flux were enabled. The fully
converged volume-rate solution gives `0.010997` aligned and `0.015203`
jittered disk masses per unit time, or `2.133e-5` and `2.949e-5` of the
disk mass over one current step. Pairwise fluxes sum to zero over these
interior faces. The p99 **net** cell mass change per step is much smaller:
`6.920e-7` aligned and `6.127e-6` jittered. Gross crossing is not net
mass loss; transport of tangential momentum and energy across shear faces
must still be tested. The [aligned](aligned_material_flux.json) and
[jitter](jitter_material_flux.json) receipts omit disk/core and
disk/exterior faces and use analytic rather than native reconstructed face
states, so they are not a native ALE result or a long-time flux estimate.

The extended [aligned](aligned_material_flux_torque.json) and
[jitter](jitter_material_flux_torque.json) receipts also carry the **advective
momentum** `mdot v_K(c)` through those same faces. For the volume-rate
candidate, the sum of absolute off-center face torques is `6.032e-4`
aligned and `1.642e-3` jittered angular-momentum units per physical time;
their signed global sums are only `4.046e-9` and `5.626e-8` because of
cancellation. A face-centroid spin-torque ledger cancels each pair's
orbital-torque change to at most `2.65e-23` and `5.29e-23` in the proxy.
This establishes the **accounting identity** a future cell-spin branch
must satisfy if it transports off-center momentum. It neither shows that
spin cures Kepler shear error nor measures a native hydro torque.
