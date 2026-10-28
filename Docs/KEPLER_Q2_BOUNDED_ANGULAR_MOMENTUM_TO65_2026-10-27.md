# Bounded q2 angular momentum through t=6.5

The [read-only auditor](audit_kepler_bounded_angular_momentum_series_20261027.py)
rechecked every native snapshot hash in the [wall-margin series](KEPLER_Q2_BOUNDED_WALL_MARGIN_SERIES_TO65_2026-10-27.json),
then summed `m[(x-3)vy-(y-3)vx]` over fixed groups selected by their **initial**
radius. The [JSON receipt](KEPLER_Q2_BOUNDED_ANGULAR_MOMENTUM_TO65_2026-10-27.json)
contains all five sampled times and exact group counts.

At accepted `t=6.5`, whole-domain angular momentum differs from `t=0` by
`-9.1753e-10`, or `-6.0186e-11` relative. The disk `0.5<=r(0)<2` changes
by about `+1.18e-12`. The 56 initial `r>=2.6` sites gain `1.1159e-5`
(`+1.8445e-4` relative within that group), while the adjacent 201 sites
with `2.35<=r(0)<2.6` lose `1.1140e-5` (`-1.0801e-4` relative).
These large local terms nearly cancel in the global sum. The observed
outermost-site motion therefore does not yet imply a comparably large **net**
angular-momentum defect.

This endpoint difference is not an integrated wall-torque or pressure-flux
ledger. It cannot separate physical intercell exchange from discrete force
and wall contributions, and it does not predict whether the shrinking
generator-to-wall clearance will remain positive at later times. The
ongoing `t=10`, `t=20`, and `t=30` restarts retain the boundary risk.
