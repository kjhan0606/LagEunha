# Static Laguerre disk-volume feasibility at native t=6.5

Independent CPU-only Slurm jobs **437352** (jitter) and **437353**
(aligned) completed `0:0` in 14 seconds each. The [auditor](../audit_kepler_bounded_weight_feasibility_t65_20261027.py)
starts from each 17-check accepted bounded native state, rehashes its `t=0`
and exact `t=6.5` snapshots, and uses the **actual** evolved generator
positions. It holds the core and exterior weights at zero and solves for
power weights of the same 49,038 fixed-mass disk cells so each cell's
power area equals its own `t=0` native Voronoi area. No generator is steered
and no mass is changed. The unweighted regular-triangulation areas agree
with native `t=6.5` disk areas within `5.95e-8` relative at worst.

| Actual t=6.5 IC | Unweighted median area error | Unweighted maximum | Weighted maximum | Maximum `|W_i-W_j|/d²` | Faces above native `0.8` guard |
|---|---:|---:|---:|---:|---:|
| Aligned | `1.01822e-4` | `2.06836e-3` | `1.47808e-11` | `0.00322` | 0 |
| 5% site jitter | `1.38710e-3` | `9.76907e-3` | `3.25877e-11` | `0.01184` | 0 |

The [aligned](aligned_weight.json) and [jitter](jitter_weight.json) receipts
record every Newton step, cell-area error, signed weight range, and the ten
most stressed power faces. Two accepted Newton corrections took 1,101
and 1,194 total Jacobi-preconditioned CG iterations, respectively, with
about 13.3 seconds of CPU wall time per case. The maximum absolute disk
weight is `2.11e-6` aligned and `2.83e-6` jittered. A common positive
gauge shift applied to **all** cells could make signed weights nonnegative
without changing any power cell; the present solve itself keeps exterior
weights at zero. The compressed [weight arrays](jitter_weight.npz) retain
the actual numerical candidate, and [file hashes](copied_inputs.sha256.json)
verify the copied receipts and arrays.

This establishes the **existence of a static power diagram** matching all
initial disk volumes for these two actual positions. It improves on the
older C010 geometry in one limited respect: at this `t=6.5` q2 state, no
disk-touching face approaches the current native `0.8` weight-difference
guard. It does not establish a continuous weight trajectory, weighted-face
velocity/work consistency, correct weighted pressure force, mass or angular
momentum conservation, GPU candidate completeness, or Kepler evolution
past `t=6.5`. The full CPU Newton solve is far too expensive to repeat at
every RK4 stage. The next useful gate is a warm-started **weight-rate**
solve using the current face graph, with truncated GPU matvec timing and
an exact native W=0 reduction check before a long calculation.
