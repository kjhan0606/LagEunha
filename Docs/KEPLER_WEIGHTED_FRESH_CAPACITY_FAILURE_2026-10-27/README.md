# Fresh weighted run: 48-corner abort and accepted restart point

CPU comparison job **437565** independently audited the completed same-IC
W0 and weighted checkpoints at `t=5.77547349` and `5.77547374`, respectively.
Their initial hashes match and all 19 W0 checks, all 21 weighted checks, and
the comparison lineage checks passed. Disk mass-weighted mean absolute
density change was `1.63149e-3` W0 and `6.03901e-5` weighted; mean absolute
entropy log change was `1.37945e-5` and `6.48025e-6`, respectively. Total
energy relative error was about `4.66e-9` for both, with zero injection.
[Comparison receipt](comparison.json) · [W0 receipt](w0_receipt.json) ·
[weighted receipt](weighted_receipt.json). These are fixed-mass trajectories;
weighted geometry changes cell density, and the result does not prove
conservative moving-face transport or the `t>6`/GIZMO stability target.

The fresh `t=0.001→6.5` weighted jitter job **437456** stopped at RK4 step
24: the CUDA tessellator reported `[GPU_CLIP_CAPACITY] nc=48 capacity=48`
while clipping a bounded cell. This is an intermediate polygon capacity
failure, not evidence of a disk instability or wall crossing. The failed
RK4 stage is not a restart state.

The code saved the step-start state at exact `t=0.045616533734626756` as
`khdiag.step000024.start_geometry.dat`, SHA-256
`35c9e11ff60471e7a68a991d16acd5b9571926cceb8e17c8d310141dc9039c7e`.
The [independent weighted checkpoint audit](accepted_step_start_audit.json)
checks the completed log only through the start-of-step marker, excluding
the failed later stage. All 21 checks pass: fixed mass, positive state and
EOS, independent power areas (`5.95e-8` maximum relative difference),
weight edge guard (`0.00533 < 0.8`), zero injection, density and orbit
gates, and total-energy relative change `3.62e-11`. Ninety-two recorded
weighted stages belong to the 23 completed steps. The checkpoint and failed
run remain preserved.

The weighted source archive's `exam_gpu.cu` already has the previously
validated near32 preclip order, but its `exam_gpu.h` retained
`GPU_MAX_CORNERS=48`. This mismatch was missed when the weighted binary was
assembled. The prior bounded W0 near32 path used a 64-corner header and
successfully replayed the analogous early failure point; that earlier W0
result does not by itself prove the new weighted state will pass.

The [isolated capacity patch](../apply_kepler_weighted_capacity64_20261027.py)
changes only the header's corner capacity from 48 to 64 while preserving
the weighted interfaces and the uniform-pressure projection fix. The
[v5 build runner](../build_kepler_weighted_capacity64_v5_20261027.slurm)
explicitly rebuilds every object that directly includes the changed header
to avoid a host/device polygon-layout mismatch. CPU Slurm **437495** builds
the candidate. It must be checked on a short native restart from the saved
step-start state before a longer continuation is released. No fresh restart
from `t=0` is needed for the failed interval.

The capacity64 v5 build **437495** completed with exit `0:0`; its binary
SHA-256 is `feae0a48f737a7eab215f2e31be1d7938f45f894376dcdfda177f7e09f57d26b`.
The projection unit tests pass and the build log records the rebuilt
host/device geometry objects. The [recovery staging script](../stage_kepler_capacity64_recovery_20261027.py)
verifies the 21-check parent receipt, exact checkpoint hash and time, initial
state, mask and v5 binary. Recovery job **437497** runs the accepted
`t=0.045616533734626756` step-start to `t=0.05`, then performs the weighted
physics audit. It must pass before the conditional continuations run.

The [chain staging script](../stage_kepler_capacity64_weighted_chain_20261027.py)
prepared **437498** (`0.05→0.825`, `afterok:437497`) and **437499**
(`0.825→6.5`, `afterok:437498`). Each continuation independently verifies
the parent audit and exact endpoint hash and time before copying the restart,
and audits its own endpoint. At submission the A10 user GPU quota kept 437497
pending; the chain has not yet demonstrated recovery or long-term stability.

Job **437497** subsequently completed `0:0` in 15 seconds, crossing the
failed stage and reaching `t=0.05` from the saved step start. Its
[recovery audit](recovery_t005_audit.json) passes all 21 checks; disk
mass-weighted mean absolute equilibrium-density change is `3.580e-6` and
total energy relative change is `3.933e-11`. The accepted endpoint SHA-256
is `7c9b3cac53e7a815c6671bc95eba60186cf9edfeb7cdf04b3f7484cff07b26c9`.
The successful dependency released **437498**, which has begun evolving
to `t=0.825`; the final `t=6.5` continuation remains conditional on the
intermediate endpoint passing its own audit. This is evidence for the short
recovery only, not a proof against later geometry failures.

The continuation **437498** reached exact `t=0.825` from that accepted
`t=0.05` state. Its [weighted audit](recovery_t0825_audit.json) passes all
21 checks, with disk mean absolute equilibrium-density change `1.20338e-5`,
mean absolute entropy-log change `9.16317e-7`, relative total-energy change
`6.65248e-10`, and zero recorded injections. Its endpoint SHA-256 is
`b929fc51043675f0bdb09eb1858dcf374cd0dbc8f9bac42e7dffc6bd12d7a6dc`.
Job **437499** remains conditional on the Slurm success of this audited
predecessor. The result does not establish the `t=6.5` or `t=600` endpoint,
physical moving-face transport closure, or cell rotation.

After **437498** completed `0:0`, **437499** verified the `t=0.825`
checkpoint and began the `0.825→6.5` segment. Two further guarded jobs were
submitted using the [long-chain stager](../stage_kepler_capacity64_weighted_long_20261027.py):
**437530** runs `6.5→10` only after `afterok:437499`, and **437531** runs
`10→20` only after `afterok:437530`. Each independently checks its parent
receipt, checkpoint hash and time, initial state and binary, and audits its
own exact endpoint before allowing the next. They retain the v5 binary and
write checkpoints every 0.825 simulation-time units, about 32 MB each.
Failing the density, geometry or conservation gates stops the chain rather
than treating an inaccurate endpoint as success.

While **437499** was running, its complete intermediate step 427 was
independently audited at `t=1.6513743395369227`. The
[live-checkpoint receipt](live_step427_audit.json) passes all 21 checks:
disk mean absolute equilibrium-density change `1.66920e-5`, mean absolute
entropy-log change `1.30298e-6`, total-energy relative change `1.33305e-9`,
and zero recorded injections. Its endpoint SHA-256 is
`796d9b4c561d1f67e461be4c5c8a3512e9c4e1275f1c4b35bb4769df6a00edea`.
The audit examines the native log only through the saved step marker; this
does not assert the state of any later RK4 stage or final endpoint.

For the same 5%-jitter initial snapshot, the older W0 job **437327** has
an independently audited step-427 checkpoint at `t=1.651374332206546`;
the time differs by `7.33e-9` from the fresh weighted checkpoint. The
[comparison receipt](w0_vs_weighted_step427.json) verifies matching initial
SHA-256 and both physics audits. Mean absolute cell-density change is
`1.46650e-3` for W0 versus `1.66920e-5` for weighted, and mean absolute
log-entropy change is `6.22665e-6` versus `1.30298e-6`. Whole-disk signed
log-entropy change is `-2.60503e-8` versus `-8.65582e-8`; neither signed
state difference is a measurement of irreversible heating.

This comparison uses different binaries and restart histories. More
fundamentally, moving weights change the power-cell geometry itself; a small
`m/V` difference at fixed cell mass does not by itself prove accurate
physical density or consistent mass/energy transport across moving faces.
It is a promising short-time numerical observation, not an established
solution of the `t>6` or GIZMO-long-time problem.

The same comparison at the next complete checkpoint (`t≈2.475808`), again
using the identical initial snapshot and time difference `2.75e-8`, is in
[step 852](w0_vs_weighted_step852.json). Both endpoint audits pass. Mean
absolute cell-density changes are `1.47095e-3` W0 versus `1.66935e-5`
weighted; absolute log-entropy changes are `7.95049e-6` versus `2.41668e-6`.
The signed disk means are `-1.00479e-7` and `-5.12498e-7`, respectively.
This corroborates the short-time trend in the recorded cell variables but
does not resolve the physical-density or local transport caveats above.

At the next audited checkpoint `t≈3.300241`, both the fresh weighted and
same-IC W0 controls pass their respective gates; the time difference is
`6.11e-8`. The [comparison](w0_vs_weighted_step1277.json) gives mean absolute
cell-density changes `1.48560e-3` W0 versus `2.80499e-5` weighted, and
absolute entropy-log changes `9.52652e-6` versus `3.87731e-6`. The weighted
absolute changes have increased from the preceding checkpoint; its
whole-disk signed entropy change is `-7.29390e-7`, compared with
`-1.76664e-7` W0. These are checkpoint diagnostics only: the calculation
continues toward `t=6.5`, and the weighted transport caveat remains open.

At `t≈4.126612`, the next [same-IC checkpoint comparison](w0_vs_weighted_step1703.json)
passes both native gates with time difference `1.09e-7`. Mean absolute
cell-density changes are `1.53155e-3` W0 and `4.66031e-5` weighted;
absolute entropy-log changes are `1.09945e-5` and `5.04696e-6`. The
weighted density difference is still smaller but has grown since the
`t≈3.30` checkpoint. These results do not establish what happens at the
historical `t=6–10` interval or validate local moving-face transport.

At `t≈4.951044`, the next [same-IC comparison](w0_vs_weighted_step2128.json)
passes both native gates. Mean absolute cell-density changes are
`1.56104e-3` W0 and `5.64549e-5` weighted; absolute entropy-log changes are
`1.24236e-5` and `5.86270e-6`. The weighted mean absolute density change has
increased from `2.80499e-5` at `t≈3.30`, but remains below the configured
`1e-3` gate. This finite-time observation is not a proof of physical
transport consistency or a prediction of the `t=6.5` endpoint.

CPU job **437546**, staged by
[this runner](../stage_kepler_fresh_weighted_compare_t65_20261027.py), uses
`afterany:437499` to compare the fresh weighted and same-IC W0 `t=6.5`
states if a complete native weighted receipt exists. It accepts either all
21 checks passing or **only** the equilibrium-density check failing for
diagnostic purposes. In the latter case it records the failed check and
sets `weighted_native_gate_passed=false`; it does not release the
`afterok` long continuation. Missing receipts, changed hashes or any other
failed checks stop the analysis for manual diagnosis. Synthetic receipt
tests confirmed failed density is rejected by default and remains explicitly
false when deliberately allowed for comparison.
