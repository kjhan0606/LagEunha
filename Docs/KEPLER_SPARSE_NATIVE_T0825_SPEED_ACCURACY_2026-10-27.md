# Corrected sparse q2 versus full q2 at native t=0.825

## Completed native gates

The periodic-area-corrected narrow (54,312 sites, Slurm **437160**) and wide
(55,931 sites, **437161**) H200 runs both reached exact `t=.825` in 426
global steps. Both ended `COMPLETED`, exit `0:0`, and passed all 16 early
native checks, including global area 36, fixed mass, W=0, positive state,
zero floor/dual/entropy injection, first-ring orbital/density/entropy gates,
and zero fixed-mass generator crossings of the inner disk edge. The separate
[narrow strict](KEPLER_SPARSE_NARROW_T0825_STRICT_2026-10-27.json),
[wide strict](KEPLER_SPARSE_WIDE_T0825_STRICT_2026-10-27.json), and state
consistency receipts independently confirm `rho=m/V` and the ideal-gas EOS.
This is an early gate, not a `t>6` stability result.

## Measured H200 cost on the same physical interval

The full q2 path from accepted `t=.001` to `.825` used 10,842 global steps:
original job **436890** to its accepted `t=.691008` step start, then the
validated halo-12 recovery **436995** to `.825`. Their native executable wall
times were `12,347.875 + 3,060.390 = 15,408.265 s`. The original job's
halo-8 geometry abort and restart are included; no accepted interval was
replayed. Narrow and wide used halo 12 from `.001` to `.825`, with native
wall times `197.845 s` and `197.260 s` respectively, about **77.9–78.1
times less elapsed H200 time for these complete early paths**. The sparse
runs took 425 continuation steps (step 1 to step 426). Fewer central
timestep constraints and fewer tessellation sites both contribute; this is
an observed path-cost comparison across *different ICs and halo histories*,
not a controlled per-kernel speedup or a forecast for `t=600`. The sparse
`t=10` jobs **437174/437175** are now independently eligible after their
passed `.825` parents and must measure evolved timestep, halo-32 cost, and
physical stability.

Decomposing the observed wall ratio, the full path took **25.51 times more
accepted steps**. Its wall divided by accepted steps was `1.421 s/step`,
versus `0.466 s/step` narrow and `0.464 s/step` wide, a further factor of
about **3.05**. These are end-to-end path averages that include the full
run's halo failure/recovery overhead; they are not matched-kernel timings.

## Matched 49,038-cell original disk

The full, narrow and wide initial disk generator coordinates and circular
velocities agree exactly; per-cell initial disk masses agree within
`6.06e-8` relative. The [same-ID auditor](audit_kepler_full_sparse_matched_disk_t0825_20261027.py)
compares those cells at the same native time. Its
[narrow receipt](KEPLER_FULL_SPARSE_MATCHED_DISK_T0825_2026-10-27.json)
and [wide receipt](KEPLER_FULL_WIDE_MATCHED_DISK_T0825_2026-10-27.json)
record exact checkpoint hashes and all radial bins.

The [native density-field comparison](KEPLER_FULL_SPARSE_NATIVE_DENSITY_T0825_2026-10-27.png)
shows all three complete `.825` checkpoints on identical spatial and color
scales. The upper panels show absolute density; the lower panels show signed
relative density change from each case's initial state. The figure is a
density comparison, not an entropy or late-time stability diagnostic.

| Initial-disk statistic, mass-weighted | Full q2 | Narrow | Wide |
|---|---:|---:|---:|
| Absolute relative density change | `1.16521e-4` | `1.16596e-4` | `1.16521e-4` |
| Absolute log-entropy change | `1.29736e-4` | `2.50389e-6` | `2.51848e-6` |
| Absolute matched density difference from full | — | `1.75508e-7` | `2.62864e-8` |
| Absolute matched generator position difference from full | — | `2.41243e-9` | `6.622996e-10` |

The entropy norm is about 52 times smaller in either sparse run while the
same-disk density field is almost unchanged. This does **not** establish
that sparse sampling removes physical shear heating: full q2 took about
25.5 times as many global steps; rounding of float internal energy, core
and outer gas sampling, and the different halo history are confounded.
The signed mean entropy changes are negative in all three cases. The
already submitted same-checkpoint energy, direct-U and fixed-K RK4 probes
are intended to isolate the thermal arithmetic contribution. The
[same-IC smaller-step control](KEPLER_SPARSE_IC_MASS_CORRECTION_2026-10-27.md#matched-native-disk-and-time-step-control-at-t0825)
**437240**, with `EUNHA_KEPLER_DTETA=0.006`, directly tests the
accumulated-step contribution without changing sparse site sampling.

The strongest same-ID contrast is the initial `0.75<=r<1` band (6,050
cells): full and narrow density drift agree within `1.94e-8`, while their
absolute log-entropy errors are `1.5652e-4` and `2.1666e-6`. In **98.98%**
of those matched cells the sparse absolute entropy residual is smaller.
The full signed mean is `-8.2354e-5`, so this early discrepancy is not a
measured positive net heating rate. Its spatial coherence makes the
controlled arithmetic probes worthwhile, but the different global step
counts still prevent a causal attribution.

Wide's matched density and position differences are smaller than narrow's
at essentially the same early wall cost. If both pass `t=10`, wide is the
better first long-run candidate, with narrow retained as a guard/sampling
control. The sparse first-ring profile-weighted mean gravity differs from
the generator value by about 30.8%; this known resolution limitation and
the outer periodic seam still require native `t>6` diagnosis. Neither
early result verifies GIZMO's `t=600` or longest-orbit claim, and no AREPO
throughput comparison has yet been measured.
