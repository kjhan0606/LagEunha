# Pressure-force growth across the native t=6.5→10 Kepler interval

The [force auditor](../audit_kepler_bounded_pressure_force_t10_20261027.py)
rebuilds the saved `W=0` Voronoi faces independently and evaluates the
native variational pair pressure-force formula on every disk site. On the
accepted `t=6.5` states, it matches the separate native GPU face-force
capture for the same 64 selected cells to maximum absolute acceleration
vector differences `5.48e-12` aligned and `1.01e-11` jittered. The
independent Voronoi areas also agree with native saved disk areas within
`1e-5` relative, as required by the auditor. These checks support the
reconstruction as a frozen-state readout of the existing CUDA force.

| Whole-disk quantity | Aligned t=6.5 | Aligned t=10 | 5% jitter t=6.5 | 5% jitter t=10 |
|---|---:|---:|---:|---:|
| Mass-weighted pressure-acceleration RMS | `1.746e-8` | `1.688e-8` | `2.507e-7` | `3.881e-7` |
| Tangential component RMS | `8.227e-9` | `8.264e-9` | `2.380e-7` | `3.696e-7` |
| Pressure standard deviation / mean | `4.679e-4` | `4.668e-4` | `3.617e-3` | `5.088e-3` |

The t=10 jitter endpoint is explicitly **failed** by the separate native
`disk_ballistic_residual` gate. Its force readout is diagnostic, not a
successful physics validation. The jittered pressure force is dominated by
the tangential component. In contrast, aligned pressure fluctuations and
force stay nearly constant over this interval.

CPU-only Slurm job **437368** completed `0:0` in 34 seconds. It read six
existing complete snapshots from each t=10 run, without restarting either
simulation. The [series receipt](pressure_force_series_t10.json) records
the exact checkpoint times and SHA-256 values:

| Native time | Aligned force RMS | Jitter force RMS | Jitter / aligned |
|---:|---:|---:|---:|
| `6.600871` | `1.696e-8` | `2.537e-7` | `15.0` |
| `7.425299` | `1.638e-8` | `2.800e-7` | `17.1` |
| `8.251666` | `1.653e-8` | `3.057e-7` | `18.5` |
| `9.076092` | `1.840e-8` | `3.407e-7` | `18.5` |
| `9.900517` | `1.802e-8` | `3.836e-7` | `21.3` |
| `10.000000` | `1.688e-8` | `3.881e-7` | `23.0` |

The jittered force and pressure variation rise together across all saved
states. The [t=10 density diagnosis](../KEPLER_Q2_BOUNDED_T10_DIAG_2026-10-27/README.md)
also finds growing density residual and site-position error, with near-zero
**signed** disk entropy change. These observations support investigation
of a geometry-induced pressure-force feedback instead of treating this
failure as demonstrated net shear heating. They do not prove a causal
force-history integral, and the six intermediate snapshots have no
separate full 17-check native physics gate. A controlled weighted or
pressure-force intervention at the same checkpoint is needed before
claiming a repair.

The four endpoint JSON files and twelve intermediate JSON files in this
directory contain the full radial summaries and source snapshot hashes.
The series receipt SHA-256 is
`a0d20cf4c993fcadce3e02f20c23a213bd8a8af086cf4709c3b3b2c6e60b0641`.
