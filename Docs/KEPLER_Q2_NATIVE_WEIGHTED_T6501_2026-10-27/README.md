# Native pressure-work-neutral weighted Kepler step, t=6.5→6.501

The accepted aligned and independent 5% site-jitter native checkpoints were
continued with the user's CUDA+MPI+OpenMP solver. The new mode keeps site
steering off, physical viscosity zero (`Re=∞`), and every cell mass fixed.
It solves 16 Jacobi weight-rate iterations on the resident GPU FaceCSR at each
RK4 stage and projects the rate so global `Pᵀ J Wdot≈0`.

The first native run, Slurm **437379**, failed at stage 0 because the existing
`GFS_VARPRESS_W0` launch guard rejected every weight rate. No evolved endpoint
was accepted. The isolated guard repair allows only the explicit
`GFS_GPU_VOLW2` mode with weighted face kinematics; the W0 branch remains
gated as before. The repaired build **437381** completed with binary SHA
`e90a287d25b0d82705f70ad645b5eadade7a18cd3aebff1622606960a333e7a1`.
The resumed one-step job **437383** completed `0:0`; both cases passed all
16 native gates. Original build, failed run, and logs remain preserved.
[The exact native source patch](../KEPLER_GPU_ENERGY_NEUTRAL_NATIVE_V2.patch)
reproduced all seven patched/new file hashes when applied to an isolated
copy of the recorded parent source files.

| Equal-time metric | Aligned W0 | Aligned weighted | Jitter W0 | Jitter weighted |
|---|---:|---:|---:|---:|
| Disk mean absolute `ΔV/V` | `3.749e-6` | `4.292e-7` | `8.606e-6` | `4.103e-7` |
| Inner `0.5≤r<0.75` mean absolute `ΔV/V` | `1.863e-5` | `1.211e-6` | `2.643e-5` | `1.162e-6` |
| Whole-disk RMS `Δlog P` | `1.277e-5` | `1.053e-6` | `2.342e-5` | `1.006e-6` |
| Absolute total-energy change / `|E₀|` over this step | `4.995e-14` | `4.580e-14` | `4.978e-14` | `4.945e-14` |

In jitter, weighted disk volume drift is 4.77% of the matched W0 result,
and inner volume drift is 4.40%. Existing pressure scatter remains: jitter
`std(P)/mean(P)` is about `0.003615` in both endpoints. Thus the one-step
result establishes suppression of *new geometric density/pressure change*,
not erasure of earlier pressure noise or a long-time cure.
The [matched jitter density field](density_jitter_t6501.png) shows absolute
surface density above and one-step relative change below. Each weighted
pixel is assigned to its certified power-cell generator, and both columns
use identical color scales.

The independently reconstructed weighted cell areas agree with native areas
to maximum relative error `3.03e-8` aligned and `5.95e-8` jitter. Maximum
absolute new weights are `4.74e-8` and `5.66e-8`. Four recorded RK4 stages
per case have projected global pressure-weight work at or below `6.56e-28`
in this short run. The GPU projection reports roughly 3–5 ms per stage;
full-step performance must be assessed after startup and output overhead.
For this single restart/output step, the original W0 wall times were 4.214 s
aligned and 4.112 s jitter, versus 4.318 s and 4.244 s weighted. These
differences include initialization and snapshot I/O, so they are not a
steady-state AREPO-speed comparison.

[Matched comparison with the original W0 CUDA binary](../KEPLER_Q2_WEIGHTED_VS_ORIGINAL_W0_437383.json),
[comparison method](../audit_kepler_q2_weighted_vs_original_w0_437383.py),
[native step gate](../audit_kepler_q2_weighted_short_replay_20261027.py), and
the archived receipts in this directory bind the conclusion to the input
and endpoint hashes. The two independent `t=6.501→10` continuations were
submitted as **437386** aligned and **437387** jitter.
The separately instrumented W0 step and the original W0 binary have
bitwise-identical `x,y,vx,vy,mass,volume,ie,den,pressure,w2` arrays at
`t=6.501` for both cases; their complete snapshot hashes differ only outside
these compared physical fields.

The global projection is not a local first-law derivation for a moving
power face with fixed cell mass. It does not prove pairwise energy/momentum
transport or GIZMO Fig. 7 stability to `t=600`. Keep this mode experimental
until the local conservation law and longer comparisons pass.

A separate midpoint snapshot estimate measures the **local** weight-driven
pressure work that the global sum hides. Its signed one-step domain work is
about `7e-15` of domain internal energy in both cases, but the p99 absolute
cell work relative to that cell's internal energy is `2.08e-5` aligned and
`3.60e-5` jitter; the jitter maximum is `1.08e-4`. The midpoint face/cell
identity closes below `1.5e-25` in absolute power. The midpoint construction
differs from the four native RK4 stage geometries, so these are estimates,
not a stage-exact local conservation verdict. [Aligned local work](aligned_local_work.json),
[jitter local work](jitter_local_work.json), and
[the method](../audit_kepler_q2_weighted_local_work_t6501_20261027.py)
record that limitation and provide a target for a later stage-level ledger.

The [native applied-rate face-slip estimate](../audit_kepler_q2_weighted_native_face_slip_t6501_20261027.py)
compares geometric face-normal speed with the analytic Kepler fluid speed on
145,671 disk-disk faces. The applied endpoint-average `Wdot` lowers the
length-weighted RMS slip to 88.6% of the W0 value aligned and 89.8% jitter;
holding the same small weights fixed barely changes it. Thus nearly all of
the 95% cell-volume stabilization is a change in *cell volume assignment*;
it does not make the moving faces material. The estimate uses midpoint
states rather than the four exact RK4 stages. [Aligned slip](aligned_face_slip.json),
[jitter slip](jitter_face_slip.json).

An [analytic face-transport proxy](../audit_kepler_q2_weighted_native_flux_proxy_t6501_20261027.py)
adds the density and face length to this slip. In jitter, the 99th-percentile
**net** mass crossing per cell per `0.001` step falls from `5.31e-5` of cell
mass at W0 to `3.40e-6` with the applied weights. Yet the gross sum of
absolute face crossings falls only to 91.3% of W0 (`1.67e-5→1.52e-5` of
total disk mass per step). This reconciles the apparent tension: weights
strongly reduce each cell's net volume/mass imbalance while most opposite
face crossings remain. Transported tangential momentum and off-center
torque are therefore still a conservation question. The pair orbital plus
spin torque identity closes to `5.3e-23` jitter in the proxy; **cell spin is
not implemented** by this calculation. It excludes disk boundary faces and
native reconstructed fluid states. [Aligned flux](aligned_flux_proxy.json),
[jitter flux](jitter_flux_proxy.json).
