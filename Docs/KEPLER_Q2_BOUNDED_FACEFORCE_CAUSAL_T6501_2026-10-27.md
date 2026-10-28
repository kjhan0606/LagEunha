# What drives the bounded Kepler orbit defect at t=6.5?

The [matched t=6.5 native comparison](KEPLER_BOUNDED_ALIGNED_JITTER_T65_REVIEW_2026-10-27.md)
showed growing orbit and density residuals under 5% independent site-angle
jitter, while mean entropy changes remained near zero. Exact circular-site
Voronoi areas predicted most of the pressure noise. The next question was
whether that pressure pattern actually accelerates the gas.

## Original-binary, one-step face probe

The [guarded stager](stage_kepler_bounded_faceforce_pair_t6501_20261027.py)
selected 32 jittered inner-disk cells with the largest absolute relative
pressure changes and 32 deterministic random inner-disk cells. It copied
the two **accepted** `t=6.5` checkpoints, identical physical parameter
files, and identical original near32 CUDA binary SHA-256 `3fc66a5f…` into
an isolated pair. A10 Slurm **437350** ran the [one-step probe](kepler_q2_bounded_faceforce_pair_t6501_a10_20261027.slurm)
to exact `t=6.501`; it completed `0:0` in 11 seconds. Both short native
replay gates passed **11 of 11** checks. GPU face capture returned 385
aligned and 382 jittered directed real-face rows for the same 64 owner IDs.
The first-call captured positions and velocities agree **exactly** with
their accepted `t=6.5` restarts; captured face pressure equals the owner's
GPU input pressure exactly. No selected face is a wall or ghost face.

The [face-force auditor](audit_kepler_bounded_faceforce_pair_t6501_20261027.py)
reconstructs the source's W0 pair force
`(P_i−P_j)[dS/2 + L(x_face−(x_i+x_j)/2)/d_ij]` on each owner. It divides
the summed force by the accepted cell mass. Independently, it integrates
each owner's `t=6.5` position and velocity through `Δt=.001` under only
the softened central gravity, then subtracts that reference velocity from
the native `t=6.501` velocity. The [machine-readable receipt](KEPLER_Q2_BOUNDED_FACEFORCE_PAIR_T6501_2026-10-27.json)
hashes both face captures and checkpoints.

| Mass-weighted RMS acceleration on 64 cells | Aligned | Jitter |
|---|---:|---:|
| Reconstructed native W0 pressure face force | 4.32147e-8 | 2.57633e-7 |
| Observed velocity change beyond gravity, divided by Δt | 4.32791e-8 | 2.57734e-7 |
| Relative vector L2 difference | 1.76% | 0.77% |
| Mass-weighted vector cosine | 0.999845 | 0.999971 |

Jitter raises the sampled native pressure acceleration **5.96-fold**.
The independent random half of the sample gives a 5.64-fold ratio, so
the result does not depend solely on deliberately selected pressure
extremes. In this one-step sample, the source's pressure force quantitatively
accounts for the **non-gravitational acceleration**. This is direct
evidence of the force pathway, not proof that it alone determines the
long-time disk lifetime or entropy production.

## How much force comes from exact-orbit cell-area variation?

The [ballistic-face-force auditor](audit_kepler_bounded_ballistic_faceforce_t6501_20261027.py)
replaces every captured cell and neighbor pressure by
`P_i(0)[A_i(0)/A_i^ballistic(6.5)]^γ`, computed from the exact circular
site trajectories. It retains the **same native t=6.5 face geometry and
W0 force coefficients**. Its [receipt](KEPLER_Q2_BOUNDED_BALLISTIC_FACEFORCE_T6501_2026-10-27.json)
therefore isolates an instantaneous pressure-source change, not a
counterfactual evolution.

| Force comparison on the same 64 cells | Aligned | Jitter |
|---|---:|---:|
| RMS force from exact-orbit adiabatic pressures | 4.13999e-8 | 1.85116e-7 |
| RMS force from native pressures | 4.32147e-8 | 2.57633e-7 |
| Predicted/native vector cosine | 0.997120 | 0.953698 |
| Native-minus-predicted vector relative L2 | 8.53% | 38.18% |

Exact-orbit Voronoi volume changes explain the **direction and a large
part** of the jittered pressure force. The 38% vector discrepancy is
substantial and grows with the native orbit/volume departure; geometric
pressure excitation is a primary mechanism, not a complete explanation.
The separate [sign-resolved entropy audit](KEPLER_Q2_BOUNDED_JITTER005_ENTROPY_SIGN_T65_2026-10-27.json)
does not show positive **mean** entropy growth at `t=6.5`; these force
audits do not measure irreversible heating, and pressure/volume work can
be reversible.

## Development implication

An internal cell-spin term that leaves the tessellation volumes unchanged
cannot remove the fixed-mass cell-volume pressure excitation identified
here. A generic Laguerre weight
response that controls cell-volume changes is the next candidate. Before
changing the native solver, its time-dependent face speed, volume rate,
pressure force and energy work must be derived together so the discrete
power identity and fixed-mass limit remain consistent. The current
selected-face result is a one-step diagnostic on one resolution and one
jitter seed; the gated aligned and jittered `t=10` continuations and a
boundary sensitivity check remain necessary before any long-time claim.

The earlier [C010 Laguerre feasibility study](KEPLER_LAGUERRE_VOLUME_FEASIBILITY_2026-10-26.md)
already established the weighted-cell area Jacobian and the correct
weight-rate/face-speed identity, but its global solve needed hundreds of
CG iterations and four disk-boundary faces violated the native guard.
Those results cannot be treated as a cheap native Kepler fix. The bounded
q2 case has a much smaller geometric pressure perturbation and merits a
**separate, local and costed** weight-response test before any integration
into every RK4 stage. The existing [operator contract](KEPLER_NATIVE_WEIGHT_OPERATOR_GATE_2026-10-26.md)
remains the conservation gate.
