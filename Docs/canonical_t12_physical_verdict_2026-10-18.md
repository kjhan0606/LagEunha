# t12 physical comparison and pressure decomposition

The from-t0 cap run has now produced a complete t=12.0006952650 snapshot.
Compare against adaptive frozen-consistent HG at t=12.0007311876. Initial
persistent-ID regions and the correct (3,3) center are used. Snapshot times
are read from trailers; hashes are in the JSON evidence. Both retain initial
cell masses and zero weights. This is not a clean convergence-order estimate:
times differ by 3.59e-5 and binaries include previously audited differences.

| Initial region/metric | Adaptive HG | From-t0 dt cap .0025 |
|---|---:|---:|
| Inner background phase RMS, degrees (wrapped) | 83.76 | 81.36 |
| Inner background pressure relative RMS | 12.575 | 11.773 |
| Inner background circular-velocity error RMS | .04686 | .04466 |
| Inner disk edge phase RMS, degrees | 4.613 | 4.352 |
| Disk core phase RMS, degrees | .01088 | .00952 |

Thus the smaller time cap from initialization does not cure the inner-region
physical failure at t12. The interior of the dense disk remains substantially
more accurate than its inner edge/background. Do not label the whole domain
accurate based on global energy conservation or backend equality.

Using gamma=float32(1.666667), split log pressure exactly into
gamma*ln(rho/rho0) plus ln[(P/rho^gamma)/(P0/rho0^gamma)]. For the adaptive
inner background, mass means are 2.2120 = .1087 + 2.1033 and RMS values are
2.3648, .4638, 2.2761 respectively. For the capped run, corresponding means
are 2.1675 = .0785 + 2.0890. The t12 pressure rise contains a dominant entropy
proxy contribution, not merely reversible density compression. RMS terms
are correlated and do not add; this is not an AV-specific causal budget.
HG, TCF, floors and integration effects require separate attribution.

Live check at 2026-10-17 19:05:12 UTC confirmed baseline t15.7863,
fixed-HG t15.5585 and from-t0 cap t12.0135 processes still present, last states
finite, no retries. No existing simulation was restarted or modified.

Evidence: canonical_from_t0_t12_comparison_2026-10-18.json and
canonical_t12_pressure_entropy_2026-10-18.json. Reproducible scripts:
analyze_canonical_from_t0_t12.py, analyze_t12_pressure_entropy.py.
Next physical gate: attribute the accumulated entropy change to existing
TCF/HG budgets or matched restart branches; do not infer that deleting AV
will stabilize the orbit, since prior reduced-AV branches worsened velocity.
