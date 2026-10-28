# Bounded q2 aligned W=0 Kepler reached t=20

Slurm **437343** completed `0:0` at exact `t=20` after 1:49:59 on one A10. It restarted the accepted `t=10` aligned W=0 state and used the native CUDA binary SHA-256 `d92b6cd3b8491678cdfafa7610608378db89ca9bd16585a52a51aaf26f3b3367`. That fast CPU/old GPU binary had previously matched all ten physical fields bitwise in a short matched-state check. The [t=20 native endpoint audit](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/fast_to20_a10_v1/audit_t20.json) passed **17/17** checks, including fixed masses and zero weights, EOS/area consistency, positive states, no floor or entropy injection, energy, orbit, density, and boundary checks. Endpoint SHA-256 is `cf2907c33718a04d1ac6f5db288115afd1d9dfe818286ffe514a099bcb4211db`.

| Mass-weighted disk diagnostic from the same initial state | t=10 | t=20 |
|---|---:|---:|
| Mean absolute density change | 1.96891e-4 | 2.28372e-4 |
| Mean absolute native-minus-ballistic density residual | 6.03158e-6 | 3.46573e-5 |
| Mean absolute entropy-log change | 1.19396e-5 | 2.08748e-5 |
| Mean position error from exact circular orbits | 3.57456e-7 | 1.35457e-6 |
| Relative total-energy change | 8.10306e-9 | 1.61036e-8 |
| Mean absolute exterior density change | 0.187909 | 0.381548 |

The disk remains within the current aligned-case gates, but its native-minus-ballistic residual has increased by **5.75×** from `t=10`; orbital position error is **3.79×** larger. The outer gas has a much larger change than the disk. The final RK4 ledger records `n_pair=20` for its pressure limiter while all injection counters remain zero. These trends do not establish a stable asymptotic state or generality under disordered particle layouts. This special polar-ring IC was already found to suppress most geometric density noise compared with independent 5% site-angle jitter.

An ID-matched radial read of the same initial and endpoint snapshots localizes the worsening. For the innermost disk annulus `0.5≤r₀<0.65` (2,323 sites), mass-weighted absolute density change rises from `3.93145e-4` at `t=10` to `6.10027e-4` at `t=20`; mean circular-orbit position error rises from `1.45294e-6` to `7.03938e-6`. At `1.75≤r₀<2` (11,934 sites), density change remains near `1.04e-4` and orbit error rises from `2.34707e-7` to `4.48950e-7`. These values use initial-radius bins and the same softened analytic circular orbit as the endpoint auditor. The inner disk therefore warrants continued attention even though the full-disk mean gate passes.

The dependency-protected `t=20→30` continuation, Slurm **437348**, started after the passing t=20 audit. It uses the separately bitwise-checked selective-wall binary SHA-256 `87d60c1cb059e1710dd2749c84af5c98ba17fe077e6fe080b9c8b42d5fd40968`. This follow-on is a further test, not yet a result. The same-time ring-coherent native control and nonaligned convergence remain necessary before extending any claim toward GIZMO's long Kepler benchmark.
