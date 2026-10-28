# From-t0 timestep cap: t10 physical accuracy still fails

Existing saved states, not restarted or newly generated simulations:

| Region and metric | Adaptive HG t=10.00207 | From-t0 maxdt=.0025 t=10.00027 |
|---|---:|---:|
| Inner background relative pressure RMS | 6.73221 | 6.30597 |
| Inner background phase RMS, degrees | 55.1856 | 54.8005 |
| Inner background circular velocity error RMS | .0382934 | .0377774 |
| Disk inner-edge relative pressure RMS | 5.40990 | 5.00117 |
| Disk inner-edge phase RMS, degrees | 1.96809 | 1.72367 |
| Disk core relative pressure RMS | .00558641 | .00512015 |

Regions follow the original t0 material IDs. Pressure RMS is relative to each
cell's initial pressure (6.3 means 630%, not 6.3%). Phase is measured against
the original exact circular orbit at each snapshot's own time and wrapped to
[-pi,pi]; it cannot count missed complete turns. Circular velocity error uses
the equilibrium velocity at the current position. Mass arrays are unchanged
and weights are zero in both snapshots.

This extends the earlier t4/t8 comparison through the historical t6–10 failure
interval. Reducing timestep from the initial state does not rescue accuracy.
The small difference between snapshot times and audited binary changes forbid
claiming a formal time-convergence order or attributing every difference to dt.
The result supports continuing the spatial/material-volume investigation.

Snapshot paths and hashes: `canonical_from_t0_t10_comparison_2026-10-18.json`.
Analyzer: `analyze_canonical_from_t0_t10.py` (run on lageunha).
No production parameters, checkpoints, or active processes were altered.
