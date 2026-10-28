# Gravity-only reference through t4 for live coupled branches

Integrated the independent gravity-only central response through two time
units from t2.0029059, saving metrics every .1. Used dt=.0005 and .00025.
Only1194 particles have nonzero initial position/velocity seeds. All other
particles have identical branch initial states and exactly zero response
under independent central gravity, so omitting them is exact for these
response norms. This reduction does not apply to hydrodynamics.

At t4.0029059, final/initial squared response metrics are:

| Metric | Gravity-only ratio |
|---|---:|
| Raw velocity | 17.490789 |
| Circular-departure velocity | .904380 |
| Omega-weighted displacement | 18.830106 |

Maximum absolute ratio differences between dt resolutions over all20 sample
times are9.31e-9,2.07e-9,9.63e-9 respectively. Orbital phase/shear can
therefore produce large raw state separation while the circular-departure
velocity metric does not grow at the final sample. Never identify large
raw norm growth alone as a hydro instability.

The live hydro continuation remains session62605, root
Codex_coupled_onset_to_t4_20261018. Last observed reference accepted step65
after restart, elapsed .1625 (absolute t2.2654059), finite/no retries.
Original longer baseline/fixedHG/from-t0-cap runs were also live when
checked at20:30 UTC, at t15.92179/t15.65524/t12.55564 respectively.
These are last observations, not guaranteed current states or physics passes.

Added analyze_coupled_onset_checkpoints.py: intersects case filenames,
checks full snapshot size before reading, matches times/IDs/masses, and
saves per-checkpoint response JSON with source hashes. It reads only the
three specific run directories and skips incomplete files. The source224
checkpoint reproduces previous coupled-response ratios. New continuation
checkpoints must be analyzed as they become available.

Evidence: canonical_coupled_orbit_series_2026-10-18.json.
Control script: audit_coupled_orbit_series.py.
Long-time hydro stability and full backend requirements remain unresolved.
