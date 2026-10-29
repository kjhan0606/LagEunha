# Provisional same-restart timestep comparison at t=10.1

The coarse and half-cap native trajectories both began from the identical accepted t10 checkpoint, reached exact t10.1 and passed all original 28 checkpoint checks. The half label refers to half the recent minimum parent step, not half the coarse median: observed median h was 0.000397564 for coarse and 0.0000476003 for half, with 258 versus 2,101 accepted steps. The 1/4-cap trajectory is still running. The same frozen H200-qualified candidate was used for both.

| Original-ID region | Coarse incremental mean absolute log-entropy | Half-cap incremental mean absolute log-entropy |
| --- | ---: | ---: |
| Core r0<.1 | 2.95535e-4 | 4.23066e-6 |
| Adjacent .1–.5 | 1.22798e-6 | 1.64851e-5 |
| Disk .5–2 | 3.01252e-6 | 6.02580e-5 |

On this common, already damaged restart and this fixed 0.1 physical interval, the smaller cap reduces the core's incremental absolute entropy change by about 70 times, while the disk's corresponding discrepancy grows about 20 times. Disk signed entropy change is -7.67e-7 versus -6.24e-6; these are not evidence for positive net heating. Initial-ID core density remains seriously inaccurate in both states. This does not repair inherited errors or establish a formal convergence order. The quarter-cap outcome is needed before assessing the trend.

The accepted legacy `ie` and `volume` fields are float32 (`VoroQ` in native `eunha.h` and the verified 576-byte offset map). Under normal round-to-nearest representation, a single final readout of the two endpoints gives an entropy-log rounding bound of approximately 1.99e-7, well below the measured disk mean 6.03e-5. This excludes **only a single endpoint readout rounding explanation** under that assumption. It does not identify whether the accumulation arises from repeated float operations, geometric evolution, physical work, or their interaction.

The existing `GFS_THERMAL64` stores authoritative U and V in double payload slots, but the current bounded wall-event controller explicitly rejects it and its required periodic transaction. `GFS_RK4_DIRECT_U` also conflicts with the event controller. These protections must not be bypassed; a compatible implementation would require full-payload rollback, energy/work and collision tests. Historical original-q1 direct-U testing reduced a short incremental entropy discrepancy by about 2x while keeping the energy equation, whereas fixed-K suppressed it further but is not a general shock-capable closure. That older experiment used different boundary/model and cannot by itself validate this bounded branch.

`GFS_NATIVE_CONST_CONFIG=1` in the H200 speed test only caches environment configuration; it does not enable double-precision thermal state. The 0.1-time comparison is an incremental diagnostic and does not qualify t120/t600 or GIZMO-long evolution.
