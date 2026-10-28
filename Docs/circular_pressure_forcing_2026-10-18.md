# Exact-circular geometry does not explain the evolved t4 pressure force alone

Reconstructed independent finite Voronoi faces for actual t2/t4 snapshots
and exact initial circular trajectories evaluated at the same times.
Reference pressure is P0*(V0/V_circular)^gamma. Compared by persistent ID
in the lab frame. Actual volume error is below8.6e-15; actual t2 pressure
force matches saved hydro-minus-Q within4.52e-8 relative L2 on r<1.9.

| Initial inner material bin | t2 | t4 |
|---|---:|---:|
| Actual pressure acceleration RMS | .00188102 | .02760875 |
| Exact-orbit isentropic acceleration RMS | .00097368 | .00104181 |
| Actual pressure work on actual departure | 1.40147e-9 | 9.23168e-7 |
| Exact-orbit force dot actual departure | 4.92780e-10 | 3.15540e-9 |
| Actual isentropic pressure work | 1.44419e-9 | 8.72372e-7 |

At t4 actual acceleration RMS is26.5 times the circular-reference value,
and94.5% of actual pressure work is the isentropic algebraic component.
Thus substantial evolved geometry/density response precedes late entropy
dominance. This does not identify a unique growing mode or assign causal
fractions: lab-frame subtraction compares forces at different positions,
and integrated response to earlier forcing may be large even if the
instantaneous circular-reference forcing is modest.

The first calculation used sqrt(x*x+y*y) for initial masks; roundoff at
r=.5 changed some membership relative to earlier np.hypot masks. Preserved
that output at Codex_circular_pressure_forcing_20261018 as superseded.
Authoritative v2 uses the same np.hypot material masks as earlier reports
and reproduces their t2 work values. No simulation was rerun for this fix.

Next experiment is continued coupled response through t4, rather than
another isolated .1 window. Started existing reference/+/- branches from
checkpoint224 at t2.1029059 toward t4.0029059, dt cap.0025, OMP4 each.
Periodic .2 checkpoints permit continuation; approximate snapshot storage
is253MB per case, excluding logs. Runtime may be tens of minutes.

Remote continuation: /gpfs/kjhan/LagForce/Codex_coupled_onset_to_t4_20261018.
Runner session62605. Preserve the live handle; do not relaunch on timeout.
Evidence: canonical_circular_pressure_forcing_2026-10-18.json.
Scripts: audit_circular_pressure_forcing.py, resume_coupled_onset_to_t4.py.
Long-term stability, the historical t6-10 failure, and backend scope remain
unresolved; t4 response is a mechanism experiment, not the final gate.
