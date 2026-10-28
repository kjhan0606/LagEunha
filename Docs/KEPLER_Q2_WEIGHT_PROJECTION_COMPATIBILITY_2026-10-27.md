# q2 moving-weight volume and energy compatibility at t=6.5

The accepted fixed-mass CUDA solver chooses a weight rate by 64 Jacobi steps
for the active disk cells. In the linearized power geometry, let
`v` be the volume rate from moving sites with weights held fixed, and let
`J` be the volume derivative with respect to the active weights. The local
volume target is `v + J wdot = 0`. Its current global correction then
projects the computed `wdot` onto `h·wdot = 0`, where `h` is the native
pressure-work gradient. These are generally independent constraints.
When `J` is invertible with inactive weights fixed, the local volume target
already determines `wdot`; exact global pressure-work neutrality can also
hold only when `h·(-J⁻¹v) = 0` at that state. A projection along `h`
restores neutrality but leaves local residual `-alpha J h` in addition to
the Jacobi residual.

At the accepted q2 `t=6.5` snapshot, an independent 148,557-face power
geometry reproduces the native stage-0 projection coefficient within
`7.3e-5` relative. For the first ring (`0.5≤r0<0.52`, 294 sites), its
predicted one-step volume-change RMS is `1.249985e-7` with the raw Jacobi
rate and `1.450015e-6` after the current pressure-work projection. The
actual native projected endpoint is `1.427963e-6` RMS; the projected
prediction has weighted cosine `0.999577` with the observed cellwise
change. A separate native `dt=1e-4` endpoint shows the same roughly linear
step-size scaling. This supports the projection as the immediate cause of
the first-ring volume residual at this state. It does **not** prove that
removing the projection yields lower long-term orbital error.

The raw rate's estimated pressure-weight work over this one step is
`-6.68257e-16`, or `-1.86522e-11` of the domain internal energy. Small
one-step work is insufficient to establish acceptable accumulated energy
drift at `t=600`. A separate correction along `J⁻¹h` keeps global work zero
but leaves first-ring predicted RMS about `1.430e-6`; that simple
constrained correction does not solve the local issue.

An isolated v11b native CUDA build adds
`GFS_GPU_VOLW2_PROJECTION=0` as a diagnostic switch. Unset or `1` retains
the accepted v5 projection arithmetic. Its default-projection Ada replay
from the same `t=6.5` checkpoint passes 21/21 gates and has all ten physical
snapshot fields bitwise equal to the accepted v5 `t=6.501` endpoint.

The raw-rate Ada replay also passes 21/21 *existing* gates. Its first-ring
one-step volume-change RMS is `1.328587e-7`, versus `1.427963e-6` with
projection: a **10.75-fold reduction**. Across all disk sites the RMS is
`1.446154e-7` versus `2.094808e-7`. The first-ring mean absolute orbital
angular-momentum increment is essentially unchanged (`6.80639e-11` raw,
`6.80451e-11` projected), as expected from the pressure field inherited at
the same starting instant. The native RK4 integral of raw pressure-weight
work is `-6.55717e-16`. The measured raw-minus-projected change in total
energy over the step is `-4.44e-15`, too close to output precision to
resolve this work cleanly. [Paired native receipt](KEPLER_Q2_RAW_RATE_NATIVE_PAIR_T6501_2026-10-27.json).

The existing `projected_weight_work` audit field only checks an absolute
`1e-10` tolerance. It passes the small but nonzero raw work rate and must
not be read as an exact-neutrality certificate in this experiment. A
separate one-rank H200 raw-rate continuation from the same accepted
`t=6.5` checkpoint is submitted as job `437608` through `t=10` to test
accumulation. The accepted projected reference reached `t=10` as job
`437530` with all 21 gates passing; its follow-on to `t=20` is job
`437531`. The current projected `t=10` disk density error is
`6.82125e-5` and orbit error is `1.466997e-7`.

The same-binary Ada comparison through `t=6.6` covers 52 native RK4
steps in each branch and passes 21/21 existing checks in both. First-ring
volume-change RMS from the shared `t=6.5` state is `5.51067e-5` with
projection versus `4.00942e-6` raw (**13.75-fold smaller**). But the raw
endpoint's first-ring pressure error relative to the initial equilibrium
is `0.00256987`, versus `0.00250273` projected; its mean absolute orbital
angular-momentum increment is `6.78106e-9`, versus `6.72768e-9`
projected. Raw volume control did not remove the inherited pressure
imbalance or its short-term torque. Native integrated pressure-weight work
over this interval is `-2.486e-14` raw; the raw-minus-projected total
energy difference is `-3.908e-13`. Neither is a t=600 energy guarantee.
[Paired t6.6 audit](KEPLER_Q2_PROJECTION_NATIVE_PAIR_T66_2026-10-27.json).

The accepted projected trajectory's first-ring pressure error grows from
`1.46257e-4` at `t=.05`, to `2.77678e-4` at `.825`, `2.56998e-3` at
`6.5`, and `3.86835e-3` at `10`. Its first-ring density error at `t=10`
is `2.32089e-3` although the diskwide density error is only
`6.82125e-5`; the current diskwide 21-check gate therefore hides this
localized growth. [Accepted radial timeline](KEPLER_Q2_WEIGHTED_INNER_RING_TIMELINE_T005_TO10_2026-10-27.json).
The pressure perturbation tracks adiabatic compression extremely closely:
at `t=10` the first-ring mean absolute entropy-log error is only
`9.62e-6`, about 400 times below its mean absolute pressure error, and
the least-squares slope of `log(P/P0)` against `gamma log(rho/rho0)` is
`1.00005`. The directly computed adiabatic pressure residual is also
`9.62e-6`. Thus this first-ring signal is mainly geometry-driven
compression, not measured shock or viscous entropy generation.
[Radial thermal timeline](KEPLER_Q2_WEIGHTED_INNER_RING_THERMAL_TIMELINE_T005_TO10_2026-10-27.json).
At `t=10` the weighted first-ring density error `2.32089e-3` is only
`1.12×` smaller than the same-IC W=0 value `2.60099e-3`, whereas the
diskwide weighted density error is `34.8×` smaller. First-ring orbital
angular-momentum error is `4.38448e-7` weighted versus `1.19716e-6`
W=0 (`2.73×` lower). The weighted method's large diskwide gain therefore
does not yet resolve its inner-boundary bottleneck.
[Same-IC radial t10 comparison](KEPLER_Q2_WEIGHTED_VS_W0_INNER_RING_T10_2026-10-27.json).
Its cyclic first-ring angular pressure-jump RMS at `t=10` is actually
`1.02015e-2`, versus `8.78160e-3` for W=0. This simple roughness proxy
does not order the observed orbital errors by itself; face geometry,
force direction and pressure-jump phase also matter.
[Same-IC angular t10 comparison](KEPLER_Q2_WEIGHTED_VS_W0_FIRST_RING_ANGULAR_T10_2026-10-27.json).
The actual native [full t10 surface-density field](KEPLER_Q2_FRESH_WEIGHTED_VS_W0_DENSITY_T10_2026-10-27.png)
and [inner-ring close-up](KEPLER_Q2_FRESH_WEIGHTED_VS_W0_DENSITY_INNER_T10_2026-10-27.png)
use the same initial snapshot and common density/residual color scales.
They show substantial suppression of diskwide W=0 checkerboard noise,
while a red-blue first-ring residual persists in the fresh weighted run.
These are q2 bounded-polygon fields, not a matched GIZMO Fig. 7 image.
In the captured `t=6.5` native first-ring face sample, all 128 directed
same-ring faces join sites at identical initial radius (maximum gap
`5.55e-16`). Their **initial pressure jumps are exactly zero**. On the
identical `t=6.5` face geometry, the evolved mean absolute pressure jump
is `0.953%` of the local initial pressure (99th percentile `2.45%`), and
its reconstructed same-ring orbital torque RMS is `7.34035e-8`.
Substituting the actual initial pressures into these faces makes that
category's direct torque exactly zero. The torque-driving same-ring
pressure differences were generated during the evolution, not present as
an initial same-ring pressure gradient. The sample includes the 32 most
pressure-perturbed sites among its 64 selected owners, so its jump
statistics are not a full-ring average.
[Native face counterfactual](KEPLER_Q2_SAME_RING_PRESSURE_ORIGIN_T6501_2026-10-27.json).
The proposed second H200 job `437609` was canceled while pending, before
using GPU time. An Ada raw continuation from its accepted `t=1` parent
followed the same branch to `t=6.5` and passed 21/21 native gates.
The separate raw-from-projected `t=6.5` H200 job `437608` was canceled
while still pending at zero runtime: it would continue an already noisy
parent, whereas the accepted early-raw trajectory is a better long-run
candidate.

A shorter same-binary Ada pair from that accepted `.825` checkpoint to
`t=1` passed 21/21 checks in each branch. Raw weights reduce first-ring
volume-change RMS from `2.15125e-4` to `5.66046e-6` (**38.0-fold**),
and keep its mean pressure error nearly at the parent value:
`2.77678e-4` at `.825` to `2.78754e-4` raw at `1`, versus
`3.17856e-4` projected. This is direct early evidence that the projection
contributes to later pressure-error growth. The first-ring mean absolute
orbital angular-momentum increment nevertheless rises slightly from
`1.57827e-9` projected to `1.72060e-9` raw; the already inherited
pressure pattern and other force errors still matter. Raw cumulative
pressure-weight work over `.825→1` is `-1.180e-13`, and the native
raw-minus-projected energy difference is `-5.954e-12`.
[Early paired audit](KEPLER_Q2_PROJECTION_NATIVE_EARLY_PAIR_T1_2026-10-27.json).

### Fresh t=0 causal pair

To separate newly generated pressure noise from the pattern already present
at `.825`, jobs `437629` (projected) and `437630` (raw) evolved the *same*
accepted jitter IC from `t=0` to `.001` on an A10 with the *same* v11b
binary. Both native audits pass 18/18. All 16 parsed initial state fields,
including IDs, stress and `w2`, agree exactly between branches; native
initial-file hashes differ in bytes outside these parsed fields. The first
ring's relative volume-change RMS is `1.80885e-6` projected versus
`8.48484e-8` raw (**21.3×**); cyclic angular pressure-jump RMS is
`4.83516e-6` projected versus `1.22258e-7` raw (**39.5×**).
The first-step orbital angular-momentum changes are only `~10^-14`, too
small for a robust long-time torque conclusion. The native energy changes
are similarly `~3e-14` relative to total energy, and the raw stage work is
`~10^-15` in absolute rate units. This establishes that the projection
creates a substantial *initial geometric/pressure perturbation* without
an inherited `.825` pressure pattern; it does not establish which branch
will have the better orbit at `t=1` or beyond.
[Fresh same-binary t=.001 pair](KEPLER_Q2_FRESH_PROJECTION_PAIR_T001_2026-10-27.json).

The stored IC pressure is exactly the same float value in all 55,351 cells,
`9.999999974752427e-7`. In exact arithmetic the face-difference gradient
`h_i=Σ_f g_{if}(P_i-P_j)` therefore vanishes and imposes no weight-rate
constraint. The isolated host helper instead projects for *every* nonzero
`h·h`, however small. At the first native stage the projected coefficient
is `alpha=-60996.1738`; the matching raw work rate is `-4.983e-17`,
implying `h·h≈8.17e-22` if their common initial GPU stage has the same
gradient and unprojected rate. That would be a per-active-cell RMS `h` of
`1.29e-13`, about one float32 epsilon times the uniform pressure. This is
a near-null direction with a finite rate correction, rather than evidence
that exact uniform pressure needs
work neutralization. Reconstructing pressure from stored `ie/volume`
produces three float32 levels spanning `2.27e-13`. The archived native
`Exam/exam.c` source explicitly recomputes `pressure=ie/volume*(Gamma-1)`
before GPU SoA upload, providing a concrete path for this small contrast
to enter the face-difference gradient. The actual device-stage pressure
has not yet been captured, so its exact values remain unverified. The
completed `t=.001→1` same-lineage pair below measures how the
initial perturbation affects orbital angular momentum.
[First-stage near-null audit](KEPLER_Q2_FRESH_PROJECTION_NEARNULL_T001_2026-10-27.json).

Both fresh same-binary branches also passed 21/21 native gates at exact
`t=1`. Relative to projected, the raw first-ring mean pressure error is
`1.01402e-5` versus `3.18317e-4` (**31.4× lower**); cyclic neighbor
pressure-jump RMS is `9.23922e-6` versus `2.03199e-4` (**22.0× lower**);
and mean absolute orbital specific angular-momentum error is
`2.62002e-9` versus `6.72061e-9` (**2.57× lower**). Raw accumulated
pressure-weight work over `t=0→1` is `-3.48126e-13`, while the two
relative total-energy drifts differ by about `3.47e-12`. Native A10
evolution wall times are comparable (`833.180 s` raw, `838.896 s`
projected for the same continuation). The clean fresh-raw `t=1` state
is therefore the next long-run parent; the
accepted early-raw branch below remains a separate mechanism control.
Against the raw branch that inherited projected evolution until `.825`,
the fresh-raw first-ring pressure error at `t=1` is **27.5×** lower and
orbital angular-momentum error **2.62×** lower. This is why the new
continuation starts from the fresh t=1 checkpoint rather than extending
the older raw t=6.5 checkpoint.
[Fresh same-binary t=1 pair](KEPLER_Q2_FRESH_PROJECTION_PAIR_T1_2026-10-27.json).

The accepted early-raw `.825→1→6.5` Ada trajectory and the accepted
projected `.825→6.5` A10 trajectory share exact `t=0` and `.825`
snapshots. At exact `t=6.5`, first-ring cyclic neighbor pressure-jump
RMS is `2.40753e-4` raw versus `6.73680e-3` projected (**28.0× lower**);
mean absolute pressure change is `3.10896e-4` versus `2.57125e-3`
(**8.27× lower**), density-change RMS is `2.34081e-4` versus
`2.08002e-3` (**8.89× lower**), and orbital specific angular-momentum
error is `5.81785e-8` versus `2.50667e-7` (**4.31× lower**).
The raw branch's integrated pressure-weight work over `.825→6.5` is
`-3.67719e-14`; its relative total-energy change is `-5.25672e-9`
versus `-5.24667e-9` projected. This supports the raw rate as the
better *current q2* long-run candidate despite its tiny nonzero weight
work. The projected baseline uses v5 and raw uses isolated v11b (the
default v11b one-step projected result was bitwise verified against v5),
so this is not a whole-trajectory same-binary pair. It also does not
show survival to `t=10` or the GIZMO `t=600` endpoint.
[Accepted t=6.5 early-raw comparison](KEPLER_Q2_RAW_FROM0825_VS_PROJECTED_T65_2026-10-27.json).

The first-ring orbital force responds to **angular pressure jumps**, not
only the absolute pressure error. At `t=1` the projected branch has a
larger mean pressure error (`3.176e-4` versus `2.788e-4` raw) but a
smaller cyclic neighbor-jump RMS (`1.99666e-4` versus `2.34003e-4`
raw). The raw branch largely retains the `.825` neighbor roughness
(`2.33654e-4`), while the projected branch initially smooths it.
This explains why the short raw orbit increment can be slightly worse
despite its better mean pressure and volume control. Much later, the
accepted projected branch's neighbor-jump RMS grows to `6.7368e-3`
at `t=6.5` and `1.02015e-2` at `t=10`; the approximate Fourier power
in modes `m≥30` exceeds `99%` at both times. The q2 first ring has
nearly uniform initial angular spacing (maximum/minimum gap `1.214`),
but the Fourier band remains only a descriptive proxy.
[Angular pressure timeline](KEPLER_Q2_FIRST_RING_ANGULAR_PRESSURE_TIMELINE_2026-10-27.json).

At the later `t=6.5→6.501` native step, the projection changes the
first-ring angular neighbor pressure jumps by RMS `4.689e-6`, whereas
the raw rate changes them by only `1.197e-7`. The projected branch's
*existing* neighbor-jump RMS decreases slightly from `6.736797e-3` to
`6.732142e-3`; the raw branch nearly preserves it (`6.736790e-3`).
Thus the projection's local volume perturbation also acts as a small
instantaneous smoother of the already large high-frequency pressure
pattern. Its removal by itself cannot be claimed to suppress orbital
torque, even though it restores much tighter volume control.
[Same-parent angular increment audit](KEPLER_Q2_FIRST_RING_ANGULAR_PRESSURE_INCREMENT_T6501_2026-10-27.json).

The projection coefficient is a global signed scalar. Locally, with
`h=JᵀP` and fixed-mass constant-entropy response, its rate correction
`δwdot=-alpha h` yields `δVdot=-alpha J h` and a pressure-rate term
proportional to `alpha J JᵀP`. For a positive pressure mode, negative
`alpha` has a smoothing sign; positive `alpha` has the opposite sign.
The native `.825→1` pair has negative RK4-effective `alpha` in all 91
steps, consistent with its observed smoothing. In `6.5→6.6`, 13 of 52
steps have positive effective `alpha` (25%). This sign change gives a
plausible *route* to intermittent pressure-noise amplification, but the
face geometry, cell thermodynamics and raw flow also evolve; sign counts
alone do not prove modal growth.
[Hash-backed native sign audit](KEPLER_Q2_PROJECTION_ALPHA_SIGN_EARLY_LATE_2026-10-27.json).

If raw rates improve local volume but develop unacceptable cumulative
energy error, the model must account for the work of changing weights or
allow consistent mass and energy flux through moving faces. The present
one-step result alone does not select between those physical closures.

Evidence: [independent Jacobi/projection audit](KEPLER_Q2_WEIGHT_PROJECTION_LOCAL_RESIDUAL_JACOBI_T6501_2026-10-27.json),
[inverse-gradient audit](KEPLER_Q2_WEIGHT_PROJECTION_INVERSE_GRADIENT_T6501_2026-10-27.json),
[accepted native step](KEPLER_FRESH_WEIGHT_LOCAL_T6501_2026-10-27/audit_t6501.json),
[first-ring native torque](KEPLER_Q2_FRESH_WEIGHTED_FIRST_RING_NATIVE_TORQUE_REVIEW_2026-10-27.md).
