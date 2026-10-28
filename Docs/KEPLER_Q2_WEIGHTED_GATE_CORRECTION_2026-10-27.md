# Weighted Kepler endpoint: equilibrium density is the physical gate

The bounded q2 Kepler run has a fixed-mass Voronoi (`W=0`) reference and
experimental time-dependent Laguerre-weight branches. The W0 native audit
compares its final cell volume with the Voronoi volume of sites advanced on
exact circular trajectories. That comparison isolates trajectory error from
unweighted tessellation error **only when both geometries are Voronoi**.

The staged auditor used by running jobs 437386/437387 instead computes
`ballistic_Voronoi_area / final_power_area - 1` and requires its disk mean
absolute value below `1e-3`. The numerator and denominator are different
cell partitions. A weighted solver could preserve every initial cell area
and the intended equilibrium density exactly, yet fail this gate because
the exact-orbit Voronoi cells shear. At `t=10` the existing W0 reference's
ballistic density change already has disk mean absolute magnitude
`1.471922e-3`. The old check is therefore a geometry counterfactual, not a
weighted-physics acceptance condition. It remains useful as a diagnostic.

[The corrected independent auditor](audit_kepler_q2_weighted_equilibrium_20261027.py)
keeps the strict native safety, EOS, independent power-area, edge-weight,
orbital, angular-momentum, energy, and projected-work checks. Its physical
density gate is the mass-weighted mean absolute change from the **initial
equilibrium density**, with threshold `1e-3`. It reports the W0 geometric
counterfactual separately without using it for pass/fail. It accepts the
two observed pressure-force startup messages only when the weighted GPU
stage sequence and nonzero final weights are also present.

The running jobs carry immutable staged copies of their original auditors;
this correction does not alter their inputs or interrupt their evolution.
After their exact `t=10` outputs exist, apply the corrected auditor to the
saved endpoints and compare matched aligned/jitter branches. A corrected
endpoint pass would establish only this finite-time equilibrium gate. Local
moving-face mass, momentum, energy and spin exchange, MPI scaling, the
paper-matched initial condition and the GIZMO long-time density image remain
separate requirements.

The corrected auditor was exercised on the saved native `t=7.426299367`
checkpoints, without changing either evolving job:

| IC | Equilibrium density change | Difference from W0 ballistic area | Corrected gate |
|---|---:|---:|---|
| [Aligned](KEPLER_Q2_WEIGHTED_EQUILIBRIUM_T7426_ALIGNED_2026-10-27.json) | `5.7030e-5` | `1.5699e-4` | All checks pass |
| [5% jitter](KEPLER_Q2_WEIGHTED_EQUILIBRIUM_T7426_JITTER_2026-10-27.json) | `1.094632e-3` | `1.380552e-3` | Only the strict density check fails |

The independent power-area maximum relative errors are `3.03e-8` and
`5.95e-8`. These are intermediate checks, not the requested `t=10` or
long-time result. The jitter result shows why the old geometric difference
and the actual equilibrium density error must be reported separately.
Both receipts verify the exact initial, restart, and executable hashes.
The stage-work ledger is cut at the checkpoint's own `step=478` log marker
(1,908 stages), so later writes to the still-running job cannot change this
intermediate verdict.
