# Weighted Jacobi iteration comparison through t=10

The native 64-iteration aligned and 5%-jitter jobs **437393** and **437394**
branched from the saved, audited 16-iteration `t=7.426299367` states and
evolved to exact `t=10`. The 16-iteration reference continued from that same
trajectory without the branch restart. The 64-iteration branch used a
different build, so the comparison is informative but is not a perfectly
controlled same-binary, same-restart experiment.

| Mass-weighted disk metric at t=10 | Aligned 16 | Aligned 64 | Jitter 16 | Jitter 64 |
|---|---:|---:|---:|---:|
| Mean absolute equilibrium-density change | `1.45954e-4` | `1.44529e-4` | `1.71685e-3` | `1.71013e-3` |
| Mean absolute entropy-log change | `9.06615e-6` | `9.62911e-6` | `1.53116e-5` | `1.53706e-5` |

The independently corrected audits are [aligned](aligned_corrected_audit.json)
and [jitter](jitter_corrected_audit.json). All 21 aligned checks pass. The
jitter run has 20/21 passing: equilibrium density exceeds the `1e-3` target.
Its total-energy relative change is `8.075e-9`; floor and other recorded
injections are zero, and the independent power-cell area check passes.
Slurm **437394** reports `FAILED 1:0` because its original post-run physics
audit returned failure, although the native MPI process completed and
produced a complete checkpoint. The corrected independent audit also fails
the equilibrium-density gate; it is a physical-gate failure, not a successful
long-time stability result.

Increasing iterations gives a small density improvement in these endpoints,
while entropy is slightly worse. This comparison supplies no evidence that
Jacobi iteration count alone fixes the jittered disk or reaches GIZMO
`t=600`. These fixed-mass time-varying-weight runs also do not establish the
missing moving-face mass/energy closure or implement the proposed cell spin.
