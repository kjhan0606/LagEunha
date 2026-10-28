# Positive hybrid density completes the fixed remap path, but mixing remains

Both N4 and N8 cases complete 0 -> .25 -> 0 with quadratic density when
continuously positive and normalized exponential quadratic otherwise. Retain
actual spatial moments, affine mean velocity and nonthermal residual covariance.
No positivity floors or added energy. Maximum exponential cells per state is
13/1590. Exponential fits use independent quadrature checks; overlap integrals
are checked at two orders. Aggregate moment gates are enforced per transition.

| Subdivisions per leg | Internal K / initial U | Velocity change mass RMS | Represented L change |
|---|---:|---:|---:|
| 4 | .0405747 | .000423183 | 5.20e-14 |
| 8 | .0405659 | .000423211 | 1.04e-13 |

For comparison, piecewise-constant projection on the same path gave internal
K/U=.0426824/.0427035, velocity change .000436798/.000436794 and relative
represented L changes about -1.66e-5. The hybrid preserves the represented
angular momentum and positivity, but reduces mixing only by about five percent.
Internal kinetic covariance remains about four percent of initial thermal
energy; it is not heat and is not a fraction of total orbital energy.

Maximum global 17-slot drift is 6.70e-14/1.33e-13 (potential slot unused).
Minimum continuous rho/mean over all states is .000619/.000795. Residual
covariance passes stated floating tolerances, not strict eigenvalue positivity.
Maximum independent fit check residual is 9.82e-13/1.06e-12. Diagnostic wall
times are 48.3/95.5 seconds; this Python clipping reference is not a production
cost or backend benchmark.

## Decision

This removes the identified density/covariance/represented-L inconsistencies
for one finite remap path. It does not remove projection mixing. Do not adopt
the hybrid as a Kepler cure or add it to every RK stage based on these results.
Two subdivision levels are not a formal convergence-order test. Before further
implementation expansion, distinguish irreducible fixed-resolution projection
mixing from errors that improve with spatial resolution or better retained
velocity structure, and relate the magnitude to the original pressure/AV
feedback. Conservation alone is insufficient for the physical goal.

The initial uniform lift is not equivalent to the original point Lz/Phi state.
Residual covariance transport is an experimental representation, not a derived
gas-pressure closure. The approved persistently equal-weight fixed-mass
reduction remains unproved. No physical time evolution, pressure/gravity,
production CUDA/OpenMP/MPI parity, or long Kepler stability is established.

Existing original trajectories were not restarted or modified. At
2026-10-17 18:11:19 UTC live processes were confirmed at baseline t15.646209,
fixed-HG t15.479874, from-t0 time-cap t11.587316.

Evidence: `canonical_hybrid_affine_path_2026-10-18.json`.
Script: `audit_hybrid_affine_path.py`.
Remote: `/gpfs/kjhan/LagForce/Codex_hybrid_affine_path_20261018`.
