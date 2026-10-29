# Net fixed-mass ALE flux and scalar weight rates on Q2 geometry

For the Q2 saved t30 and t47.5 configurations, use exact softened Kepler
velocities and a hypothetical uniform unit density. On the 3,931 original
inner-disk cells and their 12,508 faces, solve a weighted linear least-squares
problem for a scalar Wdot on each of 4,646 inner and touching-neighbor
generators. The objective is **zero net ALE mass flux per inner cell**,
not zero flux at every face point. Four-point face quadrature is used.
Wdot is unconstrained and outer neighbor rates are free; neither boundary
conditions for those rates nor a global all-cell solve is imposed.

| Diagnostic | t=30 | t≈47.5 |
| --- | ---: | ---: |
| Weighted RMS net-mass residual / uncorrected norm | 1.77e-9 | 6.50e-10 |
| LSMR iterations / stop | 775 / 1 | 888 / 1 |
| Pointwise face-slip RMS remaining | 99.992% | 99.984% |
| Face-mean slip RMS remaining | 86.97% | 92.95% |

The independently accumulated face flux and sparse Jacobian agree to relative
RMS below 4.9e-15; native stored cell areas agree with reconstructed areas
to better than 5.9e-8. The signed inward and outward face flows cancel in
each selected cell, while large **pointwise** exchange remains. Thus fixed
net cell mass is a weaker requirement than a material face. Momentum, total
energy and internal moments can still cross those faces even when their net
mass flux is zero; they cannot be dropped based on fixed mass alone.

This is a frozen analytic-velocity, uniform-density counterfactual, not a
simulation or a Riemann solution. It uses W=0 geometry only at the sampled
instant: nonzero Wdot subsequently leaves W0. It does not establish the
user-required fixed-mass behavior when W0 is maintained, rate regularity
near W0, nonuniform density, positivity, weight boundedness, whole-domain
constraint solvability, or RK4-stage/restart/MPI consistency. Its many
iterations and omitted face state/force costs also do not establish speed
competitive with AREPO. This result does **not** undo the earlier finding
that weight-speed-only mass exchange cancels the instantaneous density
effect of weight-induced volume change in a uniform-density state: this
test evaluates the **full** ALE slip and constrains its net sum instead.

The related earlier work was found in the separate historical Git branch
`origin/codex/physics-validation-20261017`, starting at commit
`9068809f9b7232fc9406abefbea42dfa7335ef29`. Its files
`Docs/staggered_fixedmass_weights_2026-10-09.md`,
`Docs/staggered_hllc_fixedmass_2026-10-09.md`, and
`Docs/staggered_moment_evolution_2026-10-09.md` describe static fixed-mass
face transport, a four-cell HLLC constraint, and a small evolved affine
moment experiment. The last reported unresolved realization and topology
errors; those short historical tests do not establish a usable Kepler solver.
These files are preserved in that branch and were read with `git show`; they
were not merged or used to overwrite the user's current native CUDA source.

The new [receipt](KEPLER_Q2_ZERO_NET_MASS_WEIGHT_2026-10-29.json) and
[reader](audit_kepler_q2_zero_net_mass_weight_20261029.py) record inputs,
solver residuals and the explicitly limited conclusion.
