# Best disk-face normal-speed match from scalar weight rates at t=6.5

The previous volume-rate solver made the signed volume change per cell
small while leaving substantial face slip. This independent CPU analysis
asks a different, optimistic question: with the native unsteered generator
positions and velocities frozen at the accepted `t=6.5` W0 checkpoint,
how much can **any scalar weight rate per disk cell** reduce the
length-weighted RMS difference between geometric face speed and the
analytic softened Kepler velocity at the face centroid?

For each disk-disk face, `s_f=(v_K(c)-w_f^W0)·n` and weight rates alter
the speed by `(Wdot_i-Wdot_j)/(2d_ij)`. The
[auditor](../audit_kepler_optimal_face_material_rate_20261027.py) minimizes
the sum of `L_f [s_f-(Wdot_i-Wdot_j)/(2d_ij)]²` on the actual 49,038-cell,
145,671-face disk graph. It fixes one common gauge, solves the sparse
positive graph system to relative gradient below `1e-8`, and imposes no
pressure-work, finite-step, boundary or GPU constraints. Thus its RMS is a
lower bound for this single-state disk-disk least-squares problem, not a
native hydrodynamic result.

| Quantity | Aligned W0 | Aligned optimum | Jitter W0 | Jitter optimum |
|---|---:|---:|---:|---:|
| Length-weighted RMS normal slip | `1.62913e-4` | `1.41220e-4` | `2.00694e-4` | `1.77424e-4` |
| Relative RMS remaining | `1` | `0.86685` | `1` | `0.88405` |
| Gross face-crossing mass proxy / disk mass per `dt=0.00193983` | `2.50751e-5` | `2.14351e-5` | `3.23225e-5` | `2.86814e-5` |
| p99 net cell mass-crossing proxy / cell mass per step | `5.90020e-5` | `2.30306e-5` | `1.03084e-4` | `2.88688e-5` |

Even this unconstrained best rate leaves **86.7% aligned and 88.4% jitter**
of the initial face-slip RMS. Net crossings fall much more than gross
crossings because opposite-face fluxes cancel. A scalar rate per cell
therefore cannot make all disk faces material at this state. A physically
conservative method must account for the residual transport or demonstrate
its convergence; freezing cell volumes alone does not remove it.

CPU-only Slurm **437445** completed `0:0` in seven seconds. Conjugate
gradient required 1,224 and 1,353 iterations respectively; that count is
not a CUDA/MPI performance benchmark and is not a proposed production
solver. The maximum one-step `|Delta W_i-Delta W_j|/d²` was `2.98e-4`
aligned and `3.29e-4` jitter, but this does not establish finite-time
feasibility. [Aligned receipt](aligned.json) · [Jitter receipt](jitter.json)
· [staged input hashes](stage_manifest.json).

The reference face fluid velocity is the analytic Kepler benchmark, not a
native Riemann state. Exterior-interface faces, local momentum/energy/spin
fluxes, pressure-work constraints, finite weight evolution, and later times
were excluded. More constrained rates cannot improve this same objective,
but a different flow reconstruction or objective must be tested separately.
