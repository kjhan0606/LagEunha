# Native centerless q2c Kepler control at t=0.825

## Run and audit outcome

H200:1 job **437179** recovered the accepted q2c step-start state at
`t=0.69101477440242343` after the earlier halo-8 clip. With halo 12 it
evolved to a complete, exact `t=0.825` checkpoint, native `EXIT:0`, in
`1971.066 s` of executable wall time. The final checkpoint is
`/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2c/halo12_after_failure_to0825/khout.006821.dat`,
SHA-256 `7381e879015ec1269c4049a4e2cdb6c50f0bf2a7bd2a0775fdb45e9975cc3124`.
The final RK4 ledger has positive internal energy, zero floor/dual/entropy
injection and relative total-energy change `-1.576e-11` in its runtime
ledger.

Slurm reports **FAILED**, exit `1:0`, because the last staged radial auditor
imported `audit_kepler_seam_checkpoint_20261026`, which was omitted from
that run's `scripts.sha256` bundle. This is a packaging error after the
native simulation and its first two auditors succeeded; it is not a
failed RK4 step or a failed physical gate. No GPU rerun was made.
The current independent [15-check strict audit](KEPLER_Q2C_T0825_STRICT_2026-10-27.json)
and [5-check radial audit](KEPLER_Q2C_T0825_RADIAL_2026-10-27.json) were run
post hoc on the preserved checkpoint and pass. They establish complete
147,456-cell IDs, exact fixed masses and W=0, positive state, box-area
sum `36.000000001775334`, `rho=m/V` within `5.91e-8`, zero injections,
disk/orbit limits and total-energy drift `-2.93803e-10` relative to the
independent initial-state energy. An independent EOS read gives maximum
relative `P=(gamma-1)U/V` error `9.77e-8` at the endpoint. The runtime
ledger and independently recomputed energy use different summation and
reference conventions; both are small and are reported separately.

## Same-time IC comparison

The original full q2 [strict receipt](KEPLER_Q2_HALO12_T0825_STRICT_2026-10-27.json)
and q2c receipt use the same physical time, solver and 384² cell count, but
q2c removes the exact central site and changes the ring allocation. It has
49,007 initial disk cells versus q2's 49,038. These are IC comparisons,
not same-cell convergence tests. Q2c reached this endpoint in 6,821 global
steps versus q2's 10,856, so its entropy difference also includes a
substantial timestep-count difference.

| Initial-disk mass-weighted measure at t=.825 | q2 | q2c |
|---|---:|---:|
| Absolute relative density change | `1.16521e-4` | `1.57610e-4` |
| Absolute log-entropy change | `1.29736e-4` | `1.03994e-4` |
| Signed log-entropy change | `-1.89861e-5` | `-7.01135e-6` |
| Generator orbit-position error | `2.50253e-9` | `3.22909e-9` |

The q2c density error is **1.35 times** q2's, while its absolute entropy
residual is **0.80 times** q2's. At `r<0.1`, q2c's absolute entropy residual
is `6.43685e-5` versus q2's `1.01538e-5`, so removing the center site did
not remove the early core error. Both disk signed means are negative; neither
is a measured positive net heating rate. The q2c exact-orbit
[geometry control](KEPLER_Q2C_BALLISTIC_2026-10-26.json) predicts disk
density drift `1.57660e-4` at `.825`, within `5.0e-8` of the native value.
The [common-scale native density fields](KEPLER_Q2_Q2C_NATIVE_DENSITY_T0825_2026-10-27.png)
show the changed ring/spoke density pattern. Their upper panels are the
absolute density field; lower panels are signed relative density change.

This early control favors q2 for lower geometric density error, but it
cannot decide whether either IC survives the historical `t>6` interval.
The queued full-q2 and corrected sparse-q2 continuations remain the
long-time decision tests.
