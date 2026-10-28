# One-stage Laguerre weight-rate geometry at native t=6.5

Independent CPU-only Slurm jobs **437354** (jitter) and **437355**
(aligned) completed `0:0` in 26 and 25 seconds. The [auditor](../audit_kepler_bounded_weight_rate_t65_20261027.py)
starts from the accepted actual `t=6.5` positions, velocities, and
[static target-volume weights](../KEPLER_Q2_BOUNDED_LAGUERRE_FEASIBILITY_T65_2026-10-27/README.md).
It computes the fixed-weight cell-volume rate from power-face normal
speeds, then solves the current face graph Laplacian
`J Wdot = -Vdot_fixed`. The exterior and core weight rates remain zero;
generator velocities are unsteered. A centered finite difference of rebuilt
power areas agrees with the face rate to relative L2 `1.84e-5` aligned and
`1.00e-5` jittered. Fully converged CG leaves weighted face-volume-rate
relative L2 below `1e-10` in both cases.

The table is a **single forward-Euler geometry prediction** with generator
displacement `dt*v` and `dt=0.00193983`, about one current global RK4 step.
Values are 99th-percentile absolute cell-area changes in the inner initial
`0.5<=r<0.75` band. They are not evolved density errors.

| Weight-rate approximation | Aligned | 5% site jitter |
|---|---:|---:|
| Hold weights fixed | `1.29929e-4` | `2.01319e-4` |
| 8 CG iterations | `2.91382e-5` | `2.99494e-5` |
| 16 CG iterations | `3.09826e-5` | `3.10339e-5` |
| 32 CG iterations | `3.01476e-5` | `3.11518e-5` |
| Fully converged | `3.01379e-5` | `3.11705e-5` |

Sixteen iterations are already close to the fully converged result for
**this one prediction**: median whole-disk area change falls from
`2.69e-6→6.20e-7` aligned and `1.01e-5→6.78e-7` jittered. The full
linear solves still require 574 and 627 **cold** CG iterations; their CPU
wall times, including repeated power-diagram builds, are about 25 seconds.
The remaining inner-band `~3e-5` error is nearly identical at 32 and full
convergence, consistent with the first-order Euler geometry step setting
the floor here. An RK4 stage-coupled update would require its own check.

The [aligned](aligned_rate.json) and [jitter](jitter_rate.json) receipts give
all 8–1024 iteration snapshots, linear residuals, area percentiles and
source snapshot hashes. The [saved rates](jitter_rate.npz) allow later
replay; [copy hashes](copied_inputs.sha256.json) verify the artifacts.
This does **not** validate native weighted CUDA tessellation, a weighted
pressure/thermal-work operator, conservation, individual timesteps,
production GPU speed, or long-term Kepler stability. The next native gate
must keep `W=0` as the fixed-mass Voronoi reduction and measure a
stage-coupled weight update at a sustained cost compatible with the user's
CUDA+MPI+OpenMP solver.
