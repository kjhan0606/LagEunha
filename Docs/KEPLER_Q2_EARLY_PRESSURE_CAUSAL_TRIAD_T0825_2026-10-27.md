# Early Kepler error separated by geometry, pressure force, and weights

The [short pressure-force-off native receipt](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/pressure_force_off_jitter_from0_to0825_a10_v1/audit_t0825.json)
starts from the **same byte-identical accepted t=0 snapshot** SHA-256
`8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b`
as the standard W=0 and fresh weighted 5% independent-jitter runs. Its
only physics source change is the previously one-step-validated diagnostic
removal of the W=0 pressure-force vector; geometric PdV, EOS, gravity,
fixed mass, mesh and timestep rules remain. The diagnostic binary SHA-256
is `0a65c373613a5b31f8449b456387ebcf80d72136a81e46eec10c9692a8bed5e7`.
Slurm **437540** was canceled at zero elapsed time because the user A10
GPU limit kept it pending. The same staged inputs ran manually on the
previously idle lageunha RTX 5000 Ada in 449.43 seconds. The
[runner derivation receipt](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/pressure_force_off_jitter_from0_to0825_a10_v1/ada_runner_receipt.json)
records the single Slurm-bootstrap removal; all 12 staged input hashes
passed. Its native run reached exact `t=0.825`, exited zero, and passed
all 9 diagnostic safety checks, including fixed mass, W=0, positive state,
EOS, area, no energy-floor injection and no geometry abort.

| Disk measure at t=0.825 | Standard W=0 | Pressure force off | Fresh weighted |
|---|---:|---:|---:|
| Mean absolute density change | `1.619961e-3` | `1.621040e-3` | `1.203383e-5` |
| Mean exact-orbit position error | `3.427501e-8` | `4.094354e-13` | `2.209657e-10` |
| Mean exact-orbit velocity error | `9.817118e-8` | `1.262710e-12` | not measured here |
| Native-minus-exact-orbit density residual | `2.157255e-6` | `1.855482e-8` | W=0 geometry inapplicable |

The standard W=0 values come from the accepted
[native t=0.825 audit](KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T0825_2026-10-27.json),
the weighted values from its 21/21 passing
[fresh t=0.825 audit](KEPLER_WEIGHTED_FRESH_CAPACITY_FAILURE_2026-10-27/recovery_t0825_audit.json).
The pressure-force-off endpoint SHA-256 is
`5edff9019af2792e1c786ddba2dddb48bf099aab434ca70aa80ad446e3440d41`.

**Causal reading.** The [exact circular-orbit Voronoi construction](KEPLER_Q2_BOUNDED_JITTER005_BALLISTIC_T0825_T65_2026-10-27.json)
predicts `1.621040e-3` mean density change without any native forces.
Pressure force off leaves that density pattern almost unchanged while
reducing orbit error by **8.37×10⁴** and native-minus-ballistic density
residual by **116×**. Together with exactly zero pressure force at t=0
and the measured 15.5× jitter/aligned pressure force by t≈0.825, this
strongly supports the sequence **differential orbital shear → irregular
Voronoi cell areas → pressure contrast → spurious tangential acceleration
and orbital deviation**. Fresh evolving weights suppress the first
geometric/pressure stage while retaining the pressure-force term, giving
135× smaller density drift and 155× smaller orbit error than standard W=0.
The [same-scale native density image](KEPLER_Q2_FRESH_WEIGHTED_VS_W0_DENSITY_T0825_2026-10-27.png)
shows the W=0 geometric pattern and its suppression by weights.
For the present point-and-scalar-weight tessellation, adding an internal
cell spin without changing the tessellation cannot remove this leading
cell-area disturbance; spin may still improve face-state shear transport
or torque accounting at later times.

The pressure-force-off trajectory is **not an accepted hydrodynamic
solution**: it breaks the pressure force/work adjoint. The standard endpoint
used A10 and the force-off endpoint used Ada, so the numerical magnitudes
are a cross-architecture comparison; the exact initial snapshot and source
intervention are controlled, and the five-order orbit contrast is much
larger than the remaining point-orbit integration error. The experiment
does not establish a locally conservative weighted operator or stability
beyond `t=0.825`. Assess the running fresh weighted `t=6.5` endpoint and
half-timestep `t=10` control before longer production runs.
