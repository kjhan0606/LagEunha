# Weight cycles expose cumulative piecewise-constant projection mixing

Actual t2 geometry, masses, velocities and thermal energies, with a uniform
spatial lift per original cell. This differs from the previous fitted-density
lift. Retain extensive m,px,py,Pi_xx,Pi_xy,Pi_yy,U throughout remapping.
No physical time evolution, forces, new production implementation or MPI test.

For each amplitude alpha, construct the positive overlap matrix A between
original and weighted cells. Forward mapping is A diag(1/V_old); reverse
mapping is A^T diag(1/V_new). Reconstruct piecewise-constant extensive density
after each leg. Apply N cycles, with alpha=.25/N. Total absolute weight-fraction
path length is .5 in every case; the final geometry always equals the original.

| Cycles | Amplitude | Final internal K / initial U | Velocity change mass RMS |
|---|---:|---:|---:|
| 1 | .25 | .0425559 | .000436828 |
| 4 | .0625 | .0425859 | .000436753 |
| 16 | .015625 | .0425936 | .000436748 |
| 64 | .00390625 | .0425955 | .000436748 |

Coverage errors are below 5.7e-14. After 64 cycles, total mass drift is
1.054e-12, total Pi component drift at most 8.70e-13, and total U drift
4.35e-17 (absolute code units). All thermal energies remain positive.
The observed mixing therefore does not originate from converting kinetic
energy to thermal energy in the diagnostic. Mean kinetic energy is redistributed
into retained internal velocity covariance. If one instead defined U solely
as total energy minus mean kinetic energy, it would count that covariance as
heat and raise pressure spuriously relative to this retained-moment model.

Retaining the original donor subcell field without projection would restore
that same field when the geometry returns. The overlap-average projection
loses this reversibility. Smaller individual remaps at fixed total variation
do not remove the measured accumulated effect for this sequence.

## Scope and consequence

This is evidence against relying on small weight increments alone to control
projection mixing. It is not a proof that affine or higher-order schemes fail.
The number of reversals changes across cases: this is NOT temporal convergence
along a single fixed smooth weight path. Internal K here is not entropy
production, and 4.26% refers to the initial thermal energy of the entire r<.7
diagnostic patch, not a percentage of total orbital energy.

The uniform spatial lift preserves the original m,p,Pi,U aggregate data but
does not prove preservation of original point angular momentum or potential
energy. Consequently these states must not be used as replacement production
checkpoints. There is no steering and no physical generator displacement in
this diagnostic. Equal final weights do not require resetting masses under
the approved persistently-equal-weight reduction condition.

Next: compare a common fixed weight path with and without affine/residual-aware
reconstruction, retaining internal kinetic bookkeeping and testing angular
momentum compatibility. A projection must be judged by convergence and its
pressure/energy interpretation, not exact pointwise donor trace matching alone.
Long-time physical Kepler stability remains open.

Evidence: `canonical_remap_cycles_2026-10-18.json`.
Script: `audit_remap_cycles.py`.
Remote: `/gpfs/kjhan/LagForce/Codex_remap_cycles_20261018`.
