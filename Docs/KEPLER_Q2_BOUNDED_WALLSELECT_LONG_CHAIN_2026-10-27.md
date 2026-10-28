# Conditional bounded W0 Kepler continuation to t=100

The aligned bounded q2 native run passed its 17-check `t=20` audit and
completed `20→30` as **437348**; the exact `t=30` receipt passed. The planned
A10 follow-on jobs were:

| Job | Interval | Slurm dependency | Required runtime gate |
|---|---|---|---|
| **437441** | `30→60` | `afterok:437348` | Passing exact `t=30` audit and matching checkpoint, initial state and binary hashes |
| **437442** | `60→100` | `afterok:437441` | Canceled before allocation after 437441 failed; no `t=60` parent exists |

Both retain fixed mass, `W=0`, the existing RK4 variational pressure
operator, zero physical viscosity, and no generator steering. They use the
same selective-wall CUDA/fast CPU binary SHA-256
`87d60c1cb059e1710dd2749c84af5c98ba17fe077e6fe080b9c8b42d5fd40968` as
the live `t=30` branch. Each requests one A10, two CPU threads, 32 GB and a
15-hour cap. The runtime stages a complete accepted parent only after its
audit passes, verifies all copied inputs, and applies the same native
17-check audit to the exact new endpoint.

Snapshots are scheduled every 2.5 **simulation time units**. An existing
55,351-cell snapshot is 31,882,244 bytes; approximately 12 and 16 interval
snapshots therefore require about 0.38 GB and 0.51 GB respectively, plus
maps and final/restart files. `/gpfs` reported 222 TB available when the
jobs were staged. The [runner](kepler_q2_bounded_wallselect_chain_20261027.slurm),
[runtime preparer](prepare_kepler_bounded_wallselect_chain_20261027.py),
and [prestager](stage_kepler_q2_bounded_wallselect_long_chain_20261027.py)
passed shell/Python syntax and input-hash checks. A temporary staging smoke
test used a copied accepted `t=20` snapshot and the same candidate binary
to exercise checkpoint time/hash, audit, input copy, and hash verification;
it ran no simulation and was not a `t=30` result. Slurm confirmed both
`afterok` dependencies after submission.

This is a favorable aligned, finite-boundary 55,351-site control. Success at
`t=100` would not establish the jittered weighted method, physical cell
spin, full 384² periodic behavior, or the GIZMO paper's `t=600` density
benchmark. Those remain part of the objective.

## Observed outcome

Job **437441** reached `Time=48.3753` (step 11580) and then exited with code
98. The exact log entry is
`[FIXED_POLYGON_CROSSING] id=55283 xy=(0.3026259641,2.4645811651) side=543 penetration=3.15984e-6`.
The preserved branch's `fixed_polygon_postStage_blend` checks sites against
the static 1,024-gon with apothem 2.75 and calls `MPI_Abort` as soon as a
site crosses it; it does not apply a reflecting impulse or account for wall
work/torque. The triggering site is at the outer domain boundary, not in the
active inner disk. The last completed RK4 ledger has
`dEtot/|Etot0|=-1.108e-8`, `n_ie_le0=0`, and zero floor injection. Therefore
the job demonstrates finite-boundary exhaustion at `t=48.3753`; it does not
provide an exact audited endpoint at that time and does not establish a
physical instability of the Kepler disk.

The dependent `t=60→100` job **437442** was canceled with zero runtime after
Slurm reported `DependencyNeverSatisfied` on failed job 437441. No output
was deleted or replayed. The static-wall branch cannot be continued toward
`t=100` without a separately justified and conservative boundary-impact
model. The full-domain periodic q2 control 437076 remains the active route
through the historical `t>6` interval.
