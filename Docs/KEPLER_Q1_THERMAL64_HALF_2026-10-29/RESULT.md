# Same-parent thermal precision comparison: t=10 to t=10.1

The H200 native thermal64 control completed at t=10.1 from the exact same accepted t=10 checkpoint as the earlier single-precision half control. Both executed 2101 accepted steps with median h=4.760027150041424e-5. Fixed masses and global IDs agree; W=0, positive finite tagged U64/V64, float mirrors, EOS, bounded domain, zero reported floor/dual-energy injection and continuous accepted times pass. The original 28 Q1 checkpoint gates also pass; those are not a core-accuracy or long-time criterion. The independently recomputed total-energy change from the common restart, using authoritative U64, is -8.30e-14 relative to initial magnitude.

| Original radial region | Cells | Single: mass mean absolute incremental log entropy | Thermal64: mass mean absolute incremental log entropy |
|---|---:|---:|---:|
| r0<0.1 | 134 | 4.23066e-6 | 2.65171e-6 |
| 0.1<=r0<0.5 | 3112 | 1.64851e-5 | 1.83206e-9 |
| 0.5<=r0<2 | 48680 | 6.02580e-5 | 1.99641e-11 |

In the disk the absolute incremental diagnostic is approximately 3.02 million times smaller in this control. The signed disk mass means are -6.23504e-6 and -2.24e-14; these numbers are not a measurement of net heating. The same-parent comparison changes the thermal arithmetic package together: authoritative U/V, pressure and work evaluation, and the RK4 anchor. It strongly establishes sensitivity to that package for this short case; it cannot by itself isolate a single float operation or prove all error is roundoff.

At t=10.1, an independent all-site fixed-polygon reconstruction agrees with the candidate's double volumes to 2.31e-10 maximum relative error, and its two guard radii agree within 1.47e-10. The zero-time U64 initial reference was promoted on Ada; its volume was independently checked at t=10 and may have small cross-GPU differences. These limits and the complete measured results are in `COMPARISON.json` and `FINAL_GEOMETRY.json`.

The checkpoint was already damaged in the core at t=10. This does not repair that inherited trajectory or establish stability from the initial state, formal temporal convergence, shocks, multi-rank behavior, t=120/600, or GIZMO reproduction. A fresh comparison must start from the preserved, already accepted early checkpoint, with core and disk assessed separately.

## Angular momentum over the same interval

The independent snapshot-and-log [angular budget](ANGULAR_BUDGET_T10_T10P1.json) covers all 2101 consecutive accepted steps, with fixed masses and W=0. Total fluid angular momentum changes by +1.6698e-13, while the recorded wall pressure angular impulse is -1.6787e-13; their sum is -8.8817e-16. The sum agrees with the printed integrator residual to 2.35e-24. There are no wall impacts in this interval, so this is not an impact-branch test. The grouping used by this audit selects disk IDs by their positions at t=10, rather than by their original t=0 radii; its disk L change (+1.3019e-11) and complementary-region change (-1.2852e-11) do not individually represent conserved quantities or measured inter-region fluxes. The audit checks accounting over this short bounded-domain continuation, not long-time core accuracy or boundary equivalence.
