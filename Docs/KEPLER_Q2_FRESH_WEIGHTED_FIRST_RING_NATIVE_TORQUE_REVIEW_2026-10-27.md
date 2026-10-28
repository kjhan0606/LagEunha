# Native weighted first-ring orbital torque at t=6.5

The [accepted face-dump replay](KEPLER_Q2_FRESH_WEIGHTED_FIRST_RING_NATIVE_FACE_TORQUE_T6501_2026-10-27.json)
restarts the fresh weighted CUDA solver from its hash-verified `t=6.5`
checkpoint and evolves one native RK4 step to exact `t=6.501` on the
`lageunha` Ada GPU. The 64 first-disk-ring sites (`0.5≤r0<0.52`) include
the 32 largest local pressure perturbations and 32 angularly spaced sites.
The first GPU force call captures 381 directed real-face rows. The
endpoint passes all **21/21** native gates, and its ten physical snapshot
fields are **bitwise identical** to the already accepted one-step run
without face dumping. Captured positions, velocities and weights match the
restart exactly; pressure differs from serialized restart pressure by at
most `1.14e-13`, and captured `pstar` equals owner input pressure exactly.

The auditor reconstructs each selected owner's fixed-weight geometric
pressure force from its actual native power face,
`(P_i−P_j)[dS/2 + L(x_face−(x_i+x_j)/2)/d_ij]`.
It independently subtracts a high-accuracy pure central-gravity orbit
from the observed native one-step velocity. The mass-weighted vector L2
force discrepancy is **0.7175%**; the corresponding orbital-torque L2
discrepancy is **0.1344%**. Thus the captured first-call geometric
pressure faces account quantitatively for the measured short-step torque
on this sample, despite the moving-weight projection elsewhere in the
solver.

| Directed face category | Rows | Reconstructed orbital torque RMS | Projection onto observed torque |
|---|---:|---:|---:|
| Same initial disk ring | 128 | `7.34035e-8` | `0.998525` |
| Different disk ring | 129 | `1.82841e-9` | `0.000500` |
| Disk/core interface | 124 | `5.08267e-9` | `0.001125` |

The projection is a dot-product fraction, not an independently additive
energy or force budget; category vectors may cancel. Together with the
small reconstruction residual, it shows that **same-ring pressure faces
are the dominant direct source of the first-ring orbital torque** in this
one-step sample. The core faces exert sizable *total force* but little
orbital torque. The core/weight-mask interface could still generate the
pressure contrast that appears on same-ring faces; this experiment does
not identify the upstream source of that contrast.

The separate [all-cell one-step angular-momentum audit](KEPLER_Q2_FRESH_WEIGHTED_ONE_STEP_ORBITAL_TORQUE_GRAVITY_REF_T6501_2026-10-27.json)
shows the weighted first-ring mean absolute relative increment is
`6.80451e-11` per `0.001` step, with `96.6%` mass-weighted sign alignment
to the existing error. The independent pure-gravity RK4 reference is only
`6.91e-17`; the observed increment is not plausibly a central-gravity
time-integration truncation effect. A spatial pressure-force channel is
consistent with the separate finding that halving the W=0 timestep from
`t=6.5` to `t=10` did not improve its density or orbital error.

The physical decision is to audit and suppress the first-ring **along-ring
pressure imbalance** while keeping the native mass, momentum and energy
ledgers closed. A spin degree of freedom can account for off-center pair
torque, but conservation of orbital-plus-spin angular momentum by itself
does not remove a spurious pressure imbalance or prove long Kepler
stability. This result covers 64 selected sites and one step; it does not
establish the full-ring integrated torque, its temporal origin, or `t=600`.
