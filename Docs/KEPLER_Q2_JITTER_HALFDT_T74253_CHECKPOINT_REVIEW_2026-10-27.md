# Kepler jitter timestep control at the first equal-time checkpoint

The [read-only checkpoint comparison](KEPLER_Q2_JITTER_HALFDT_CHECKPOINT_T74253_V2_2026-10-27.json)
uses the same native W=0 CUDA executable (SHA-256 `3fc66a5f...f75cd55e`),
the same t=0 initial snapshot, and the identical accepted t=6.5 restart
(SHA-256 `d02e15eb...1a9f78d`). The only intended integration change is
`EUNHA_KEPLER_DTETA=0.06` versus `0.03`. The complete saved checkpoints are
at `t=7.4252993734` and `7.4253000173`, differing by `6.44e-7`.
This is an intermediate checkpoint, not either run's t=10 terminal gate.

| Initial-radius band and mass-weighted measure | Standard dt | Half dt | Half / standard |
|---|---:|---:|---:|
| Disk 0.5–2: mean absolute density change | 1.837827889e-3 | 1.837828015e-3 | 1.000000069 |
| Disk 0.5–2: mean orbit position error | 7.448352700e-6 | 7.448354173e-6 | 1.000000198 |
| Inner 0.5–0.75: mean absolute density change | 2.081794493e-3 | 2.081793471e-3 | 0.999999509 |
| Inner 0.5–0.75: mean orbit position error | 9.996794775e-6 | 9.996796825e-6 | 1.000000205 |

Halving the step **after t=6.5** therefore has no measurable beneficial
effect on the already accumulated density or orbit error through t≈7.4253
at the precision of these means. The disk mean absolute entropy-log change
is `1.64065e-5` versus `1.67383e-5` (half/standard `1.0202`); neither
value by itself measures irreversible heating. The common t=6.5 history
means this control cannot assess timestep error accumulated before t=6.5.
The requested t=10 endpoint audit remains necessary to judge subsequent
growth, and one timestep pair cannot establish convergence order.

The ID-aligned cellwise comparison rules out cancellation in the headline
density and orbit means: the disk mass-weighted mean absolute **difference
between the two density fields** is `7.37e-9` relative to each cell's t=0
density, and the difference between their analytic-orbit residual vectors is
`6.66e-11`. Entropy-log states differ by a larger `5.36e-6` on average,
so the almost identical density and orbit metrics do not imply identical
thermal updates. Raw velocity vectors differ by `4.78e-7`; this comparison
does not subtract the expected velocity change over the `6.44e-7` snapshot
time offset and is not used as a timestep-error measure.
