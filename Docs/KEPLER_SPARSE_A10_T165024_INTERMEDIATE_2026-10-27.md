# Sparse q2 A10 runs: first complete intermediate native checkpoint

The independent wide (55,931 sites, job **437253**) and narrow (54,312
sites, job **437254**) CUDA+MPI+OpenMP runs saved their first post-`.825`
snapshots at actual times `1.6502423913857924` and
`1.6502423887805227`, respectively. These are output checkpoints of the
ongoing `t=10` runs, not restarts or new evolution jobs. Their absolute
time difference is `2.61e-9`; the nominal dump clock is `1.65`.

| Initial `0.5<=r<2` disk, mass weighted unless stated | Wide | Narrow |
|---|---:|---:|
| Native absolute relative density change | `1.24537046e-4` | `1.24962637e-4` |
| Exact-orbit Voronoi density change | `1.24587354e-4` | `1.24587353e-4` |
| Pointwise native-versus-orbit density residual | `3.17051431e-7` | `9.36094480e-7` |
| Native-versus-orbit position error | `1.05130419e-8` | `1.76853385e-8` |
| Signed logarithmic entropy change | `-1.93383e-6` | `-1.69342e-6` |
| Absolute logarithmic entropy change | `4.06547e-6` | `3.90231e-6` |
| Disk angular-momentum relative change | `-4.95e-13` | `-2.47e-14` |
| RMS radial speed / circular speed | `1.20425e-8` | `5.74618e-8` |
| Whole-box energy relative change | `2.32878e-8` | `2.32735e-8` |

The [wide](KEPLER_SPARSE_WIDE_STATE_T165024_A10_2026-10-27.json) and
[narrow](KEPLER_SPARSE_NARROW_STATE_T165024_A10_2026-10-27.json) independent
state receipts pass all eight checks: complete finite positive native
state, exact fixed mass, W=0, periodic area 36, `rho=m/V`, and ideal-gas
EOS. The [wide](KEPLER_SPARSE_WIDE_NATIVE_BALLISTIC_T165024_A10_2026-10-27.json)
and [narrow](KEPLER_SPARSE_NARROW_NATIVE_BALLISTIC_T165024_A10_2026-10-27.json)
geometry receipts verify complete checkpoint hashes and safe bounded disk
polygons. Entropy, angular momentum, radial speed, and energy in the table
were computed directly from the two immutable native snapshots using the
same formulas as `audit_kepler_sparse_area_t10_20261027.py`; they do not
depend on the evolving run's latest log line.

The wide pointwise density residual increased from `8.47e-8` at `.825` to
`3.17e-7`, and narrow from `2.31e-7` to `9.36e-7`. These remain much
smaller than the geometric density change, but the residual growth is
real and warrants checking at `t>6` and `t=10`. The mean signed entropy
changes are negative; this checkpoint does not show positive net heating.
Independent [wide](KEPLER_SPARSE_WIDE_A10_T1650_RADIAL_2026-10-27.json)
and [narrow](KEPLER_SPARSE_NARROW_A10_T1650_RADIAL_2026-10-27.json)
radial audits find about `0.2653` absolute density change in cells within
`0.2` of a box face. The outer periodic seam is already disturbed; its
influence on later disk evolution remains an open causal question.
The [cell-ID provenance audit](KEPLER_SPARSE_A10_SEAM_PROVENANCE_2026-10-27.md)
finds that 20 initially seam-near generators moved inside `r=2.75`, but
none entered the `r=2.35` guard or the disk at this checkpoint.
It does not establish late-time stability or justify the longer GIZMO
comparison. Both `t=10` jobs should continue under their existing audit
contract; no running or pending job was modified for this check.

The [native density-field comparison](KEPLER_SPARSE_A10_NATIVE_DENSITY_T165024_2026-10-27.png)
uses the same `[-2,2]` view and color scales for both cases. Its upper
panels are actual native surface density and its lower panels are signed
relative density change in the original disk. The plotter allowed only
the measured `2.61e-9` endpoint-time difference through an explicit
`--time-tolerance 1e-8`; its default remains `1e-10`. This view omits the
box faces, so it cannot visually exclude the separately measured seam
disturbance.
