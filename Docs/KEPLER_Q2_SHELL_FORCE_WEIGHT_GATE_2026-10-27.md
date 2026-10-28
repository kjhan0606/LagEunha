# Bounded q2 Kepler: shell forces and the weighted-face gate

## What the native W0 force probe says

The accepted original-binary `t=6.5 -> 6.501` probe reconstructs the
pressure acceleration of the same 64 inner-disk cells to within 1.8%
(aligned) and 0.8% (5% independent site jitter) of their measured
gravity-subtracted velocity change. The [shell decomposition](KEPLER_Q2_BOUNDED_SHELL_FACEFORCE_T65_2026-10-27.json) classifies each
captured face by its two cells' **initial** polar radii and sums vector
forces before projecting them on the measured acceleration. It reproduces
the earlier face-force vectors to below `1.4e-22` absolute per component.

| 64-cell selected inner sample | Aligned | 5% site jitter |
|---|---:|---:|
| Same-shell face rows | 128 | 128 |
| Cross-shell face rows | 249 | 246 |
| Same-shell force projection on measured acceleration | 23.3% | 86.8% |
| Cross-shell force projection | 74.5% | 13.1% |
| Core-interface force projection | 2.1% | 0.1% |

The minimum initial cross-shell radial separation in this face sample is
`0.0220971`; the same-shell chord near `r=1` is `0.0108846`. Thus the
adjacent-shell centers are not closer than same-shell neighbors. The face
projection reverses between aligned and jittered ICs even though their
radial shell spacing is the same. The simple close-cross-shell explanation
is ruled out for this sample. Cross-shell differential rotation could still
create **along-shell pressure differences**, whose force is then delivered
through same-shell faces; the table is a force-path audit, not a causal
decomposition of pressure formation or irreversible heating.

## What a volume-rate weight response can and cannot establish

The separate [weight-rate calculation](KEPLER_Q2_BOUNDED_WEIGHT_RATE_T65_2026-10-27/README.md) solves
`J Wdot = -Vdot_fixed` at the actual accepted `t=6.5` state. For one
forward-Euler geometry step of `0.00193983`, 8 CG iterations bring the
inner-band p99 area-change estimate from `1.299e-4` to `2.914e-5` aligned,
and from `2.013e-4` to `2.995e-5` jittered. Sixteen iterations give
`3.098e-5` and `3.103e-5`, already close to the fully converged results.
The remaining `~3e-5` floor is consistent with the first-order geometry
prediction. This is a promising bounded-cost candidate, but the full
**cold** solves take 574/627 iterations and cannot be used per RK4 stage.

The [face-slip calculation](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/weighted_face_slip_t65_cpu_v1/jitter_slip.json)
shows a necessary limit: the 16-iteration volume-rate solution reduces the
length-weighted mean absolute difference between geometric and analytic
Kepler-material normal face speeds only from `9.669e-5` to `8.254e-5`
aligned and from `1.246e-4` to `1.138e-4` jittered. In the jittered
cross-shell subset it changes `1.459e-4` to `1.390e-4`. The rate equation
constrains each cell's **signed sum** of face-volume fluxes; it cannot force
every face speed to equal material speed. Therefore volume preservation
alone is not evidence that facewise shear/slip has disappeared. It may
nevertheless suppress the fixed-mass pressure source if the complete
weighted dynamics are conservative.

## Local face-velocity fit is a weaker candidate for this state

An independent [linear face-rate auditor](audit_kepler_q2_bounded_material_rate_t65_20261027.py)
uses **only the current native generator velocities** to set each face's
target normal speed to the projected pair mean. It fits the weight-rate
differences by length-weighted least squares on all 148,557 finite faces
touching the disk, with core/exterior rates fixed at zero. The softened
analytic Kepler velocity is evaluated afterward only as an independent
comparison. Results for the same `t=6.5` geometry are in the
[aligned](KEPLER_Q2_BOUNDED_ALIGNED_MATERIAL_RATE_T65_2026-10-27.json)
and [jitter](KEPLER_Q2_BOUNDED_JITTER_MATERIAL_RATE_T65_2026-10-27.json)
receipts.

| One linearized step, 16 CG iterations | Aligned volume constraint | Aligned face fit | Jitter volume constraint | Jitter face fit |
|---|---:|---:|---:|---:|
| Disk-disk face RMS slip to local pair-mean speed | `9.266e-5` | `9.092e-5` | `1.031e-4` | `9.656e-5` |
| Disk-disk face RMS slip to analytic Kepler speed | `1.436e-4` | `1.466e-4` | `1.812e-4` | `1.844e-4` |
| Disk p99 `|dt Vdot/A|`, `dt=0.00193983` | `1.538e-6` | `1.952e-5` | `1.512e-6` | `5.108e-5` |

The direct face fit buys only a **1.9% aligned / 6.4% jittered** reduction
in slip to its own local target relative to the volume-constrained rate.
It worsens the independent analytic face-slip metric and increases the
linearized disk p99 area change by **12.7x / 33.8x**. The 8/16/32-step
face-fit slip values are nearly saturated, although its 32-step linear
residual remains `0.83% / 0.46%`; no claim of exact least-squares
convergence is made. All candidate weights stay below the existing 0.8
face guard in this one-step prediction.

These `dt Vdot/A` values are first-order rates; the previous full
forward-Euler power-diagram rebuild left `~3e-5` inner-band p99 area
change even with a converged rate. The local target is a first-order
pair-mean velocity, not a
reconstructed Riemann contact speed. These limitations prevent a hydro
accuracy claim, but they make a **volume-constrained 8/16-iteration solve**
the better next native experiment than a standalone face-speed fit.

## Next native gate

Attach a bounded 8/16-iteration device-side weight-rate solve to an
isolated native RK4 restart from the accepted `t=6.5` states. Use the
same `Wdot` in power-face normal speed, cell-volume rate and thermal
work; derive pressure force from that same weighted geometry. The
[weight-power audit](KEPLER_Q2_T10_AND_WEIGHT_POWER_GATE_2026-10-27.md)
now shows that this direct fixed-mass update alone leaves an unbalanced
`-P^T J Wdot` energy term. Account for that term with a derived conservative
transport or weight-energy law before enabling weighted hydro. Audit
one-step power area and work identities, pairwise
momentum/energy and mass accounting, W=0 reduction, positivity, and
per-stage GPU wall time against the original binary. Continue only the
accepted weighted checkpoint chain toward `t=10`; do not infer long-time
stability from the geometry-only numbers here. A cell-spin mechanism may
then be tested against residual facewise shear without using physical
viscosity in this Re-infinity Kepler case.

The shell-force evidence covers one 64-cell sample; the weight-rate and
face-slip evidence use actual whole-disk geometry but **no weighted gas
evolution**. The closed polygon and polar layout remain controlled ICs,
not a matched replay of the published long-time GIZMO disk.
