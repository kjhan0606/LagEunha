# Fresh weighted jitter Kepler pilot from t=0

The [native t=10 comparison](KEPLER_Q2_WEIGHTED_T10_EQUILIBRIUM_2026-10-27/README.md)
enabled weights only after a W0 history to `t=6.5`. The 5%-jitter disk's
mass-weighted density change was `1.698279e-3` at that parent state and
`1.716845e-3` at weighted `t=10`, while its W0 `t=10` reference reached
`2.374457e-3`. The continuation reduced further growth but cannot tell
whether evolving weights from the original IC would prevent the inherited
error. A fresh start tests that distinct question; it does not discard or
replace the accepted W0 path.

The isolated [stage script](stage_kepler_q2_weighted_fresh_jitter64_t001_20261027.py)
copied the accepted jitter glass, initial areas, parameter file, disk mask,
and exact `t=0` native snapshot into a new run directory. It uses the
64-iteration opt-in CUDA weighted binary SHA-256
`7358c54996b820a7b45bd53bcf6f788f125754b3158845d2b4265afead570f38`.
The input SHA list, stage manifest and import check passed before submission.
Site steering and physical viscosity remain off, and one A10/one MPI rank
is requested. The [runner](kepler_q2_weighted_fresh_jitter64_t001_20261027.slurm)
stops at exact `t=0.001`. It writes the true initial dump and the one-step
endpoint, then invokes the [independent pilot audit](audit_kepler_q2_weighted_fresh_pilot_20261027.py).
The audit requires matching initial physical fields, exact fixed masses,
positive state and EOS, independent power areas, zero floor injection,
four weighted RK4 stages and small one-step energy/density change.

Slurm **437431** completed `0:0` in 10 seconds. Its
[native t=0.001 audit](KEPLER_Q2_WEIGHTED_FRESH_JITTER64_T001_2026-10-27.json)
passed all 18 checks: the 10 initial physical arrays match the accepted W0
initial state exactly, masses stay fixed, all four weighted RK4 stages run,
the independent power-area maximum relative error is `5.96e-8`, the disk
mean absolute density change is `9.44e-8`, and relative energy change is
`3.50e-14`. Complete initial snapshot hashes differ outside the compared
physical arrays. The work directory is
`/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_weighted_fresh_jitter64_t001_a10_v1`.
Only a passing one-step gate supports a longer fresh weighted continuation.
Neither branch yet supplies a local moving-face flux/spin conservation law.

The [conditional continuation runner](kepler_q2_weighted_fresh_jitter64_to65_20261027.slurm)
was independently staged with 12 hash-checked inputs and submitted as
**437432** with Slurm dependency `afterok:437431`. At runtime it also checks
every pilot audit flag and the exact checkpoint hash before copying the
`t=0.001` state. It then targets exact `t=6.5` with the same 64-iteration
weighted CUDA binary. The [equilibrium auditor](audit_kepler_q2_weighted_equilibrium_20261027.py)
uses the fresh `t=0` initial state, the exact `t=0.001` parent, independent
power geometry, and the `mean |rho/rho(0)-1| < 1e-3` density gate. Its
work directory is
`/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_weighted_fresh_jitter64_to65_a10_v1`.
After the exact-zero pressure failure was repaired and independently tested
on CUDA, the unstarted **437432** was held and canceled at zero elapsed
time. Its staged files are preserved. The replacement
[v4 runner](kepler_q2_zero_pressure_v4_fresh_jitter64_to65_20261027.slurm)
was submitted as **437456** from the passing v4 `t=0.001` checkpoint; its
runtime gate additionally verifies the v3/v4 physical-field comparison.
The [v4 prestager](stage_kepler_q2_zero_pressure_v4_fresh_jitter64_to65_20261027.py)
verified the GPU zero-gradient test receipt, inputs and exact v4 binary, and
an isolated restart-guard smoke test passed. At 12:06 KST **437456** was
pending the per-user A10 GPU limit. No simulated interval was repeated.

Job **437456** subsequently failed at `t≈0.0456165` because its weighted
build retained a 48-corner intermediate polygon buffer. The
[failure and restart record](KEPLER_WEIGHTED_FRESH_CAPACITY_FAILURE_2026-10-27/README.md)
preserves the independently audited accepted step-start checkpoint and
the 64-corner repair in progress.
