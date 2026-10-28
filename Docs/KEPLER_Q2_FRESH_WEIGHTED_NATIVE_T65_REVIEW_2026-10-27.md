# Fresh weighted native Kepler reaches t=6.5

The original 5%-independent-jitter bounded-q2 IC was evolved from its
accepted `t=0.825` weighted checkpoint to **exact `t=6.5`** with the user's
CUDA+MPI+OpenMP solver (one MPI rank, one A10). Slurm **437499** completed
`0:0` in `1:20:22`. The [native weighted endpoint audit](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_capacity64_v5_recover_to65_a10_v1/audit_t65.json)
passes **21/21** checks: fixed cell masses, positive state, independent
power-area reconstruction (`5.95e-8` maximum relative error), valid
bounded polygon, no recorded energy injection or wall crossing, active
weight evolution, and relative total-energy error `5.25e-9`. Its endpoint
SHA-256 is `bb94041e862efd8a16dd86c2ae8bc0d93475aa03f1a1ac9f8e333b0483f2af23`;
the maximum observed weight difference divided by squared site separation
is `0.01019`, below the solver guard of `0.8`.

The [same-time comparison](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/fresh_weighted_compare_t65_v1/comparison.json)
verifies the identical initial snapshot SHA-256, exact same physical time,
both accepted native endpoint audits, and entropy lineage. Job **437546**
completed the comparison `0:0`. All measures below use the initial
`0.5≤r<2` disk sites and initial cell masses as weights.

| Native t=6.5 disk measure | W=0 Voronoi | Fresh weighted Laguerre | Improvement |
|---|---:|---:|---:|
| Mean absolute cell-density change | `1.69828e-3` | `5.01080e-5` | `33.9×` smaller |
| Mean absolute circular-orbit position error | `5.21458e-6` | `3.60814e-8` | `144.5×` smaller |
| Mean absolute entropy-log change | `1.49389e-5` | `6.91683e-6` | `2.16×` smaller |
| Relative total-energy error | `5.24739e-9` | `5.24667e-9` | approximately equal |

The [common-scale native density image](KEPLER_Q2_FRESH_WEIGHTED_VS_W0_DENSITY_T65_2026-10-27.png)
colors pixels by the actual W=0 Voronoi or weighted power-cell winner.
Its top panels show surface density; lower panels show each cell's
`Σ/Σ_initial−1` on an identical color scale. Both runs pass their native
gates. This is the first **freshly evolved** weighted result at the
historical `t>6` boundary, not a reweighting of an already damaged W=0
state. The conditional, hash-guarded Slurm **437530** was released for
`6.5→10`, followed by **437531** for `10→20` only if each predecessor
passes its own physics audit.

The comparison supports a real reduction of early orbital as well as cell
density error. It does not yet establish survival through `t=10`, the
[GIZMO Figure 7 `t=600` benchmark](KEPLER_GIZMO_FIG7_8_BENCHMARK_CONTRACT_2026-10-27.md),
or general conservative moving-face transport. The one-step
[local-work/PdV audits](KEPLER_Q2_WEIGHTED_FRESH_LOCALWORK_T4952_REVIEW_2026-10-27.md)
show tight adiabatic cell-energy closure, while the
[material-face proxy](KEPLER_FRESH_WEIGHTED_MATERIAL_FACE_T4952_REVIEW_2026-10-27.md)
still leaves about 90% of the W=0 face-speed slip. The present polar-ring
IC and fixed polygon wall also differ from the published GIZMO setup.

Applying the paper's [Figure 8 current-position annulus formulas](KEPLER_GIZMO_FIG8_METRIC_PAIR_T65_2026-10-27.json)
to both accepted states gives `L1(v_phi)=7.78e-8` W=0 versus `1.00e-9`
weighted in `0.8<R<1.2` (10,453 sites each). The literal paper
`L1(Σ)=mean|Σ−1|` is about `0.01` for both because the initial plateau is
`Σ=1.01`; the additional own-initial-density relative mean is
`1.85e-3` versus `2.75e-5`. These annular early-time numbers help keep
future comparison formulas consistent, but the annulus can miss inner-ring
damage and this q2 IC is not a matched reproduction of the paper's setup.

## Local orbital angular momentum

The [same-IC, hash-checked cellwise audit](KEPLER_Q2_FRESH_WEIGHTED_ANGULAR_MOMENTUM_T65_2026-10-27.json)
uses each initial disk site's orbital specific angular momentum
`ell=(x-3)*vy-(y-3)*vx` as its own circular-orbit reference. It checks that
the initial and endpoint hashes match both accepted receipts, that both
endpoints have the same physical time, and that cell masses are unchanged.
For the initial-ID disk `0.5≤r0<2`, the mass-weighted mean absolute
`|ell/ell0-1|` is `8.74743e-7` for W=0 and `1.13578e-8` for moving weights,
a **77.0×** reduction. The unweighted p99 falls from `3.27869e-6` to
`3.03866e-7`; in the inner `0.5≤r0<0.65` bin, the weighted mean absolute
error is `3.95024e-8`, **16.3×** below W=0. The disk's signed integrated
changes of order `10^-11` in both runs would conceal these local errors if
used alone. This is orbital center-of-mass angular momentum only: the audit
does not measure cell spin, face angular-momentum transport, or long-time
stability. The inner-bin tail and later `t=10/20` behavior remain relevant
when deciding whether a spin degree of freedom improves the native method.

The [inner-ring split](KEPLER_Q2_FRESH_WEIGHTED_ANGULAR_MOMENTUM_INNER_RINGS_T65_2026-10-27.json)
locates most of the weighted inner-bin tail in the **first disk ring**:
`0.5≤r0<0.52` has 294 cells and weighted mean absolute relative orbital
angular-momentum error `2.50667e-7` versus W=0 `6.03420e-7` (only 2.41×
lower). In `0.52≤r0<0.55`, the weighted error is `3.63862e-8`; in
`0.55≤r0<0.65`, `4.07770e-9`. This is a real orbital error in the native
endpoint, unlike the face-flux proxy's omitted core-interface contribution.
The first disk ring's contact with the sparse core and the active-weight
mask is a focused mechanism to audit before attributing this residual to
bulk shear or deciding to enable spin globally. The
[all-face mass-flux proxy](KEPLER_FRESH_WEIGHTED_FLUX_PROXY_T4952_REVIEW_2026-10-27.md)
shows that an apparent first-ring mass leak disappears when disk-core faces
are counted, so orbital error needs a separate native force/torque audit.
