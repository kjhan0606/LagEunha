# Spin torque implied by one reciprocal variational-pressure decomposition

The active bounded CUDA pressure operator accumulates a **directed** face
contribution. It does not itself store a cell spin or specify a unique
reciprocal traction. This analysis constructs one reciprocal decomposition
that leaves every closed cell's total native pressure force unchanged, then
measures its moments at the saved face centroids. The
[auditor](../audit_kepler_variational_spin_pair_torque_20261027.py) uses
independently rebuilt Voronoi/power faces and hash-verified native snapshots.

For a face from cell `i` to `j`, let `A=L n` point from `i` to `j`, `d` be
the generator distance and `c` the actual face centroid. The native directed
force on `i` is `F_i^dir=(P_i-P_j) L(c-x_i)/d`. Define the pair force

`F_ij = F_i^dir - P_i A`, `F_ji = -F_ij`.

The removed `-P_i A` sums to zero over a closed cell. The reciprocal pair
therefore preserves the original cell force. Its moment about generator `i`
is

`(c-x_i) × F_ij = -P_i (c-x_i) × A`.

For a closed polygon, `sum_f (c_f-x_i) × A_f = 0`. With one scalar `P_i`
per cell, this particular decomposition gives **zero net pressure torque
about every cell generator**, even when neighboring cell pressures differ
and even for the weighted power faces. Individual face torques need not
vanish. Adding their moments as a spin source under this decomposition
therefore cannot by itself change the computed Kepler evolution.

| Native checkpoint | Sum of absolute net cell torques / sum of absolute face torques |
|---|---:|
| [W0 aligned t=6.5](w0_aligned_t65.json) | `6.84e-13` |
| [W0 jitter t=6.5](w0_jitter_t65.json) | `4.57e-13` |
| [W0 jitter t=10](w0_jitter_t10.json) | `4.35e-13` |
| [Weighted aligned t=10](weighted_aligned_t10.json) | `6.19e-13` |
| [Weighted jitter t=10](weighted_jitter_t10.json) | `4.35e-13` |

Each calculation covered 49,038 initial-disk cells and all finite faces
touching them. Independent power areas agree with the native areas; the
closed-cell moment residual is at most `2.04e-17`. Reciprocal-pair and
orbital-plus-spin torque identities pass. On the W0 `t=6.5` states, the
reconstructed net pressure acceleration agrees with the independently
captured native GPU values for 64 selected cells to maximum absolute vector
errors `5.48e-12` aligned and `1.01e-11` jitter. The weighted reconstruction
uses the source formula but has not yet been compared with a weighted native
GPU force capture. The two `t=10` jitter snapshots failed their respective
native density-related gates; they are diagnostic inputs, not accepted
physical endpoints.

The pair decomposition is **not unique**: a reciprocal, divergence-free
force circulation on mesh cycles can leave all generator forces unchanged
while changing the moments assigned to individual cells. Thus this result
does not prove physical pressure torque is zero, or that cell spin cannot
help. It shows that the existing generator force alone does not determine
the new spin source. A physical spin branch needs a specified face traction
or conservative momentum flux, a cell centroid and inertia, and an energy
ledger coupling translation, spin and thermal energy. At `Re=∞`, any
tangential momentum exchange must follow the inviscid transport/traction
law rather than an added friction coefficient. Those ingredients and the
long-time native test remain open.
