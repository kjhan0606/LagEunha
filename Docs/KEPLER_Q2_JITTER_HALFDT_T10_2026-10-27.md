# Independent-jitter Kepler timestep convergence from t=6.5

The accepted 5% independent site-jitter W=0 `t=6.5` endpoint, SHA-256
`d02e15ebd01193f1da002d87d689f83525b0150541a9f78dedbc451702c5a2ff`,
is the common restart for the existing `EUNHA_KEPLER_DTETA=0.06` path and
the newly submitted half-timestep control. The [isolated A10 runner](kepler_q2_jitter_halfdt_to10_a10_20261027.slurm)
sets `EUNHA_KEPLER_DTETA=0.03`, preserves the original native binary and
other physical controls, requests `--gres=gpu:A10:1`, and targets exact
`t=10`. It does not restart from `t=0`. The stage manifest checks the
parent's 17/17 audit, checkpoint hash and exact time before job submission;
the runtime follow-on helper repeats the accepted-state checks.

Slurm **437512** was first submitted to H200 but remained pending under the
per-user `QOSGrpGRES` limit behind older H200 work. It was canceled with
zero elapsed time before any evolution. The equivalent A10 job **437524**
was then submitted from the same accepted checkpoint. Its endpoint will be compared
with the existing standard-step jitter result using disk density change,
native-minus-exact-orbit density residual, orbit error, signed/absolute
entropy and energy. The common `t=6.5` history lets the comparison test
whether the `t=6.5→10` growth is sensitive to timestep. The existing
standard-step comparison was also run on A10, avoiding a GPU-architecture
change. No endpoint result is claimed until the A10 job and independent
audit finish.

The A10 job started on `syn03`. Its runtime stage record verifies the accepted
parent endpoint SHA-256 and native binary SHA-256. Early RK4 ledger steps use
`dt≈9.699e-4`, half the standard run's `≈1.940e-3`, confirming the intended
time-step perturbation rather than merely an environment-variable change.

## Exact t=10 result

The half-step simulation reached exact `t=10` with native `EXIT:0` and
endpoint SHA-256 `9d106ddae0752ee85472bb30d0591ebb1248f839a6a355caaf802839be6703e4`.
Its [native endpoint audit](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/jitter005_halfdt_to10_a10_v1/audit_t10.json)
passes 16/17 checks. As in the previous standard-step W=0 result, the
single failed check is `disk_ballistic_residual`; this intentional physics
gate makes the Slurm job report `FAILED 1:0` despite a normal simulation
exit. It is **not** a CUDA or positivity crash.

| Mass-weighted disk measure at t=10 | Standard dt | Half dt |
|---|---:|---:|
| Mean absolute relative density change | `2.37445662e-3` | `2.37445752e-3` |
| Native-minus-exact-orbit density residual | `1.19668389e-3` | `1.19668631e-3` |
| Mean circular-orbit position error | `1.62485180e-5` | `1.62485395e-5` |
| Mean absolute entropy-log change | `2.03301852e-5` | `2.31079160e-5` |
| Total-energy relative drift | `8.08e-9` | `5.34e-9` |

The density residual differs by only about **2.0 parts per million** and
the orbit error by about **1.3 parts per million** relative to the
standard-step values. The half-step absolute entropy-log measure increases
by about **13.7%**, so halving the step is not uniformly beneficial even
among these diagnostics. The [hash-checked local orbital angular-momentum
comparison](KEPLER_Q2_JITTER_T10_HALFDT_ANGULAR_MOMENTUM_2026-10-27.json)
gives a full-disk mean absolute relative per-cell error of `1.67060338e-6`
standard versus `1.67059213e-6` half dt, a difference of only `6.74e-6`
relative. Its first-ring values are likewise `1.19715755e-6` versus
`1.19731648e-6`.

This control starts from the *same t=6.5 W=0 checkpoint*, so it isolates
the `6.5→10` time-step sensitivity, not errors accumulated before t=6.5.
Within that scope, RK4 temporal truncation is not the dominant source of
the growing density, orbit, or local angular-momentum error. The next
effective levers are spatial geometry, face/stress coupling, and boundary
behavior rather than simply decreasing `DTETA` again.
