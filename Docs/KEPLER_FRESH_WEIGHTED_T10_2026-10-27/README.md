# Same-initial-state Kepler comparison at t=10

The W=0 Voronoi and evolving-weight Laguerre runs share initial snapshot
SHA-256
`8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b` and
both audited endpoints are exactly `t=10.0` (`Δt=0`). The W0 branch has one
failed check: `disk_ballistic_residual = 1.1966839e-3`, above its `1e-3`
threshold. Its other native checks pass. The weighted branch passes all
21 of its own native checks, including fixed individual mass, positive
states, independently reconstructed power area, density/EOS, zero energy
injection, total energy, orbit and angular momentum.

| Disk metric at t=10 | W=0 Voronoi | Evolving weights | Weighted / W=0 |
|---|---:|---:|---:|
| Mass-weighted mean absolute density change | `2.37446e-3` | `6.82125e-5` | `0.02873` |
| Mean circular-orbit position error | `1.62485e-5` | `1.46700e-7` | `0.00903` |
| Mean absolute log-entropy change | `2.03302e-5` | `9.09539e-6` | `0.44738` |
| Relative total-energy error | `8.07655e-9` | `8.07152e-9` | — |

At this checkpoint the weighted branch's density change is about 35 times
smaller and its orbital-position error about 111 times smaller than W0. The
absolute entropy-change measure is lower by about 55%. The signed entropy
changes have opposite signs (`-5.08e-7` W0, `+9.98e-7` weighted), so these
state differences do not by themselves identify irreversible heating. The
weighted run's difference from the W0 ballistic **Voronoi-area** diagnostic
is `1.48670e-3`; that is diagnostic, not a weighted-branch gate, because the
power cells are not required to reproduce W0 areas.

The density-field figure uses common colors and the actual Voronoi or power
cell winner at each pixel. The bottom row plots each cell's density relative
to its own initial density. It shows the large cell-scale W0 noise and the
smaller, organized residual in the weighted state.

This weighted run fixes individual masses and evolves Laguerre weights, but
it is **not** a rotation run: its launcher unsets `GFS_LAGUERRE_ROTATION`,
and that option only corrects face-normal velocity on off-midpoint Laguerre
faces. No per-cell spin, torque or rotational-energy state was evolved.
This comparison also does not establish local moving-face mass, momentum and
energy consistency or long-time stability. The exact-time comparison and
receipts preserve the one failed W0 gate rather than hiding it. The weighted
continuation from `t=10` to `t=20`, A10 job **437531**, completed with exit
code 0 and passed the exact `t=20` audit. The receipt reports fixed mass at
every site, initial `w2=0` followed by evolved nonnegative weights,
mass-weighted disk density change `1.29209e-4`, orbit-position error
`1.17038e-6`, relative energy error `1.60557e-8`, and zero recorded energy
injection. The first `t=20→30` continuation, A10 job **437611**, was
cancelled while dependency-pending. Its replacement, H200 job
**437613**, uses the same t=20 restart provenance and audit gate, requests one
H200 GPU and eight OpenMP threads, and depends on 437531's successful t=20
audit. Since both the GPU and thread count differ from the earlier A10
segment, this measures end-to-end H200 throughput and is not an isolated OMP
scaling comparison. The receipt still does not establish local moving-face
conservation, cell spin, or GIZMO-long-time stability.
