# Bounded Kepler: native KDK and individual-step readiness

This is a source audit of the user's CUDA+MPI+OpenMP solver, not a new
simulation result. The active bounded Kepler jobs use global RK4. The
following observations refer to the preserved accepted near32 source and
must be rechecked if that source changes.

## What the existing KDK path actually does

`Exam/KH/kh.c` selects `exam2d_vph_kdk_int_blend` only when
`GAS_EVOLMETHOD` requests the KDK path. The implementation at
`Exam/exam.c:11549` has two global `updateDenW2Pressure2DBlend` calls,
two global `getAccVoro2DBlend_GPU` calls, and an additional end-of-step
`exam2dUpdateVol`. The half-kick, drift, second half-kick, pressure
refresh, and diagnostic loops all traverse `VORO_NP` with no active-cell
time-bin mask. Every step also makes global `MPI_Allreduce` calls for the
timestep and physical diagnostics. This is a **global-step KDK**, not
individual time stepping. The source comment's "2 force evals versus 4"
does not measure whole-step speed because each force/update call can
rebuild geometry, and the final volume refresh remains.

The current KDK path computes one global `Dtime` from signal speed,
native face candidate `dt_dummy`, and `kepler_dt_limit`. Therefore the
same innermost core generators that set global RK4 `dt≈0.00193984` would
still constrain a KDK run until a separate conservative block scheduler
exists. Merely changing the integrator selector cannot realize the
frozen-state [26.24-fold face-work arithmetic opportunity](KEPLER_Q2_BOUNDED_FACE_BUDGET_2026-10-27.json).

## Minimum conservative block-step contract for the present W=0 mode

1. At a face evaluation, drift both active and inactive generators to the
   **same physical time** and build the Voronoi face from those positions.
   The face speed and bounded-wall work must use that same time.
2. Evaluate each active-neighborhood face once. Apply its opposite
   pressure impulses to the two cells, with a pending impulse for any
   inactive neighbor. In the fixed-mass W=0 mode there is no mass flux;
   energy/work increments and central-gravity work still require their
   own time-centered bookkeeping. A wall face contributes a wall impulse
   and torque, with zero wall speed until an impact model is validated.
3. Before an inactive cell becomes active or participates in a new face
   prediction, its pending momentum and thermal changes must be reflected
   in that prediction. Otherwise the pair can conserve a ledger while
   producing a force from stale velocity, pressure or volume.
4. Restrict CUDA tessellation and force work to active cells plus a measured
   neighbor halo. A CPU scheduler wrapped around full-domain tessellation
   would keep nearly all of the current wall cost.
5. Keep a common-time snapshot gate for area tiling, fixed mass, W=0,
   positive EOS, total energy including wall work, angular momentum plus
   wall torque, and the disk's exact-orbit density/entropy/velocity residuals.
   Compare global KDK with the accepted RK4 branch before accepting a
   block-step physical result.

The [native local-step budget](KEPLER_Q2_BOUNDED_LOCAL_DT_BUDGET_2026-10-27.json)
has only eight of 55,351 cells within twice the current minimum step at
both `t=.001` and `t=1.651374`. Its 26.71-fold active-cell and 26.24-fold
face arithmetic ratios are **not** runtime predictions. They exclude
active-neighborhood mesh construction, asynchronous registers, pressure
prediction, GPU launch overhead, MPI synchronization, and topology change.
The ongoing bounded `t=6.5` and `t=10` physics gates should establish
whether this branch is worth accelerating before the local-step operator
is promoted into a long physical run.

## Measured long-horizon cost after the first t>6 gate

The unchanged-binary A10 job 437317 advanced the accepted bounded state
from `t=.825` to exact `t=6.5` in 1h20m48s and passed its native 17-check
gate. This is about 854 seconds per physical time unit over that interval.
Holding that measured throughput fixed would put `t=6.5→600` at roughly
5.9 days and `t=6.5→2222.36` at roughly 22 days of A10 compute, before
queue delays, restart overhead, a changing timestep, or any failure. These
are resource scale estimates, not predictions of physical stability.

The queued fast-CPU and selective-wall variants are measured short-window
improvements of order 5–20%, not the order-of-magnitude reduction needed
for routine longest-duration comparisons. A conservative active-face
block-step implementation remains the main available path to that speed
scale, contingent on the ongoing jitter and later-time physics gates.
