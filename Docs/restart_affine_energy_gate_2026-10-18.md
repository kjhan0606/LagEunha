# Existing point-state to affine-cell initialization gate fails for naive replacement

Using actual saved t0 masses, velocities, thermal energies and closed polygon
geometry, compute each uniform-in-cell mass covariance S. Infer neighboring
velocity gradients at geometric mass centroids. If the old mean momentum is
retained and the within-cell velocity is replaced by an affine field, the
additional kinetic energy is K_aff=tr(G S G^T)/2.

This is a hypothetical representation change, not a modification of any state.

| Original material region | sum(K_aff)/sum(U_old) | Cells K_aff>U_old |
|---|---:|---:|
| Inner background | .34959 | 63/804 |
| Disk inner edge | 6.92340 | 390/390 |
| Disk core | .80692 | 2347/10304 |

At fixed cell momentum and total energy, U_new=U_old-K_aff. Therefore a direct
replacement produces negative thermal energy in all tested inner-edge cells
and22.78% of core mass. Keeping U_old instead adds total energy. Neither is an
acceptable unannounced restart conversion or proof of the approved reduction.

This is not a new physical heating measurement in the original run. The old
point convention has no stored internal velocity variance. A piecewise-constant
lift has central cross/velocity moments C=T=0; adding a new resolved gradient
changes its kinetic content.

As already noted in the 2026-10-09 formulation, a neighbor gradient can instead
be merely a basis for the same field, but then cross and residual terms are
essential:

    C_res=C-SG^T,
    R_actual=T-GC-C^T G^T+GSG^T,
    T=GSG^T+G C_res+C_res^T G^T+R_actual.

For the old constant lift, these terms cancel the apparent added kinetic
energy. That bookkeeping identity does not justify using the new affine
velocity as the actual face trace while discarding the residual field. Also,
negative T-GSG^T alone does not prove all reconstructions impossible.

## Decision

Do not connect the new affine face/moment components to existing checkpoints
using naive affine replacement. Component backend tests remain valid, but a
positive, conservative state mapping and the approved original Voronoi
reduction are prerequisites for a coupled trajectory. Do not fix this by
raising pressure, adding energy, discarding moments at zero weights, or
restarting the benchmark with silently changed initial conditions.

Alternative mappings or cross/residual-aware closures remain research options;
this audit does not prove that all moment models are impossible. Evidence:
`canonical_restart_affine_energy_2026-10-18.json`.
