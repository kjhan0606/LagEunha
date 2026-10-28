# First continued checkpoint at exact matching control time

All three live branches have complete checkpoint263 at t2.2004059177633986,
duration .1975 from the original t2 seed. This is not the .2 duration in the
regular gravity-series table: output scheduling uses absolute simulation
time. Added compare_coupled_checkpoints_orbits.py to integrate controls at
actual checkpoint durations, with no time interpolation, and cache results
for completed checkpoints.

| Squared response ratio from initial | Hydro | Gravity only | Hydro/gravity |
|---|---:|---:|---:|
| Raw velocity | 1.630134 | 1.661207 | .981295 |
| Circular-departure velocity | 1.131144 | 1.043303 | 1.084194 |
| Omega-weighted displacement | 1.712766 | 1.792023 | .955772 |

Maximum ratio change between gravity dt=.0005 and .00025 is8.85e-11.
Hydro velocity even/odd response ratio is.001112. At duration .1, the
circular-departure hydro/gravity ratio was1.004812; it is now1.084194.
This motivates following the fluid contribution beyond the shortest probe,
but does not establish an unstable mode or additive causal energy fraction.
Raw velocity/displacement metrics remain below the gravity-only values.

The original t4 continuation is still running; do not interrupt it.
Next independent checks should target this completed .1975 duration:
half-amplitude +/- seeds and half-dt reference/+/- seeds, retaining the
same initial geometry/entropy reconstruction and comparing full response
fields, not just norm ratios. Existing reference/checkpoints can be reused
where initial state, binary, time and parameters match. Hydro convergence
is not yet established by the accurate gravity-only control.

Evidence: canonical_coupled_t2p2_comparison_2026-10-18.json.
Checkpoint outputs and cached exact-time comparisons are under
Codex_coupled_onset_to_t4_20261018/response_frames and orbit_comparisons.
Long-time physical accuracy, t6-10 recovery and full backend scope remain
unresolved.
