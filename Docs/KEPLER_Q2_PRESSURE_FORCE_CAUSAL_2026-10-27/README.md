# Causal pressure-force control from the failed jittered Kepler branch

The original CUDA+MPI+OpenMP bounded q2 jitter run passed its exact
`t=6.5` native audit but failed the `disk_ballistic_residual < 1e-3`
gate at exact `t=10` with `1.196684e-3`. The independent, GPU-cross-checked
[pressure-force reconstruction](../KEPLER_Q2_BOUNDED_PRESSURE_FORCE_T10_2026-10-27/README.md)
finds a monotonic increase in its tangential pressure acceleration and
pressure variation during `6.5→10`.

For a causal diagnostic, the [source patch](../KEPLER_DIAG_PRESSURE_FORCE_OFF_2026-10-27.patch)
sets only the CUDA `GFS_VARPRESS_W0` pressure-force vector to zero. It leaves
the geometric face-volume work, internal-energy update, EOS, gravity,
generator velocity, timestep and mesh rules in place. The new binary is
identified by SHA-256
`0a65c373613a5b31f8449b456387ebcf80d72136a81e46eec10c9692a8bed5e7`;
the [build receipt](BUILD_RECEIPT.json) records the original and candidate
source/object hashes. This intentionally breaks the physical pressure
force/work adjoint. Its output can isolate the force's effect but cannot
be accepted as a hydrodynamic solution.

The one-step A10 job **437370** restarted the byte-identical accepted
5%-jitter `t=6.5` checkpoint SHA-256
`d02e15ebd01193f1da002d87d689f83525b0150541a9f78dedbc451702c5a2ff`
and evolved to exact `t=6.5005`. The original executable's previously
saved one-step run used the same physical inputs. Job 437370 completed
`0:0` in 7 seconds; both one-step safety audits pass fixed mass, `W=0`,
positive state, EOS, polygon tiling, no floor injection and no wall
crossing. The [paired receipt](one_step_pair.json) gives:

| One-step measurement | Value |
|---|---:|
| Disk RMS velocity difference, original vs force-off | `1.2534508e-10` |
| Difference divided by `dt=0.0005` | `2.5069016e-7` |
| Independent original `t=6.5` pressure acceleration RMS | `2.5068547e-7` |
| Acceleration ratio | `1.0000187` |
| Disk RMS position difference | `3.1340e-14` |
| Disk mean absolute density difference | `2.4582e-12` |

Thus the implemented switch removes precisely the independently measured
pressure acceleration at one step to about `1.9e-5` relative, while the
geometry and density have not yet changed appreciably. The short test does
not establish that removing pressure force repairs the later density error.

The same accepted `t=6.5` state was staged again for the A10
`t=6.5→10` diagnostic **437372**, rather than restarting from `t=0` or
from the clipped short endpoint. It completed `0:0` in 48m33s. Its
[runner](../kepler_pressure_force_off_jitter_to10_a10_20261027.slurm),
[native audit](audit_t10_native.json), and
[paired comparison](causal_comparison_t10.json) preserve the original failed
result. All 10 diagnostic safety checks pass, but both the original and
force-off runs fail the same `disk_ballistic_residual < 1e-3` physics gate.

| Exact `t=10` measurement | Original pressure force | Pressure force off |
|---|---:|---:|
| Disk mass-weighted absolute native-minus-ballistic density residual | `1.196684e-3` | `1.135361e-3` |
| Disk mass-weighted absolute density drift | `2.374457e-3` | `2.321232e-3` |
| Disk RMS radial speed / circular speed | `1.297121e-6` | `6.440487e-7` |
| Disk mass-weighted absolute entropy log change | `2.033019e-5` | `2.033633e-5` |

Removing the pressure force lowers the disk-wide residual-density metric
by only **5.12%**, despite roughly halving the radial-speed metric. The
[radial audit](radial_tails_pressure_off_t10.json), compared with the
[original radial audit](../KEPLER_Q2_BOUNDED_JITTER005_RADIAL_TAILS_T10_2026-10-27.json),
shows a 19.38% residual reduction at `r=0.50–0.75`, falling toward zero
near `r=1.50`; at `r=1.75–2.00` the residual instead rises 3.40%.
Pressure force therefore influences the inner orbit error but does not
account for most of the disk-wide density residual in this intervention.
The [independent point-force continuation from the actual t=6.5 state](../KEPLER_Q2_JITTER_INHERITED_ORBIT_T10_REVIEW_2026-10-27.md)
now shows that the pressure-force-off `t=10` orbit error is almost entirely
inherited: its native positions differ from point-force integration by only
`3.14e-12` on average, and independently reconstructed Voronoi volumes
predict its `1.135361e-3` density residual. Therefore this late force-off
intervention does not identify the force that created the error before
`t=6.5`. This comparison does not uniquely identify the remaining cause. Since
the pressure force/work adjoint is deliberately broken, this branch
cannot be promoted as a physical hydrodynamic solution. The physical
Kepler solver still needs a conservative mechanism, potentially evolving
Laguerre weights and a cell-spin ledger, before any long-time claim.

The separate [same-t0 intervention to t=0.825](../KEPLER_Q2_EARLY_PRESSURE_CAUSAL_TRIAD_T0825_2026-10-27.md)
answers a different, earlier causal question. Removing pressure force from
the exact accepted initial state leaves the predicted geometric density
change near `1.62e-3`, while mean orbital position error falls from
`3.43e-8` to `4.09e-13`. Thus the late `6.5→10` force-off result must not be
read as evidence that pressure force was unimportant during the earlier
development of the jittered orbit error.
