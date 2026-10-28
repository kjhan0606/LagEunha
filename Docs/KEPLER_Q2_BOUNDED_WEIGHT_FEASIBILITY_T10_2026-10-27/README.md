# Static Laguerre volume feasibility at native t=10

Independent CPU Slurm jobs **437376** (aligned) and **437377** (5% site
jitter) both completed `0:0` in 14 seconds. They apply the same static
power-cell target-volume Newton solve used at accepted `t=6.5` to the
actual native `t=10` generator positions. Each of the 49,038 disk cells
targets its own original `t=0` area. Core and exterior weights stay zero;
no generator is steered and no cell mass is changed.

The aligned t=10 native parent passed all 17 physics gates. The jittered
t=10 parent **failed** `disk_ballistic_residual`; the static geometry
solve does not erase or override that failure. Both input checkpoint
SHA-256 values are rechecked before solving.

| Static t=10 geometry | Aligned | 5% site jitter, failed native parent |
|---|---:|---:|
| Maximum relative target-area error | `1.47833e-11` | `1.78252e-10` |
| Maximum absolute disk weight | `2.57786e-6` | `3.46814e-6` |
| Maximum `|W_i-W_j|/d²` over disk-touching faces | `0.002869` | `0.017056` |
| Faces over current native `0.8` weight guard | `0` | `0` |
| Newton/CG iterations over two corrections | `535+586` | `610+586` |

The [aligned](aligned_weight.json) and [jitter](jitter_weight.json)
receipts contain the full nonlinear history, source checkpoint hashes,
face-guard extrema and native failure labels. Their saved
[aligned](aligned_weight.npz) and [jitter](jitter_weight.npz) arrays contain
the static weights. The result shows that the jittered t=10 failure is
**not** due to loss of a static target-volume power diagram: a positive,
well-conditioned geometric solution still exists at those actual
positions. It does not show that such weights can be evolved cheaply or
conservatively, nor that a weighted native hydro run would avoid the
failure. The cold CPU solves are far too slow to repeat per RK stage.

## Frozen-state pressure-force counterfactual

The [auditor](../audit_kepler_bounded_static_weight_pressure_t10_20261027.py)
uses the native GPU-validated variational face-force formula on the same
`t=10` site positions. It compares original W=0 geometry and pressure,
static weighted geometry with the **old** pressure, and static weighted
geometry with pressure recomputed adiabatically as
`P_new=P_native (V_native/V_power)^gamma` at every disk cell and its 1,443
directly adjacent inactive cells. Mass and entropy labels are held fixed.

| Disk pressure-acceleration RMS | Aligned | 5% jitter |
|---|---:|---:|
| W=0, native pressure | `1.68799e-8` | `3.88092e-7` |
| Static weights, native pressure | `1.68798e-8` | `3.88006e-7` |
| Static weights, adiabatically reweighted pressure | `2.20240e-9` | `1.05352e-8` |

For interior `0.75≤r₀<1.75` disk cells the jitter force falls from
`4.01628e-7` to `1.69093e-9` with the reweighted pressure, about 238×.
The disk pressure standard deviation/mean falls from `0.0050875` to
`0.00002679`. Changing face geometry while leaving the inconsistent old
pressure unchanged does almost nothing. This strongly localizes the
**frozen-state** force error to the EOS response to area noise. The
[aligned](aligned_pressure_counterfactual.json) and
[jitter](jitter_pressure_counterfactual.json) receipts preserve native
failure labels and input hashes. They do not implement weighted RK4,
pressure-work conservation, mass flux, or a long-time cure.
