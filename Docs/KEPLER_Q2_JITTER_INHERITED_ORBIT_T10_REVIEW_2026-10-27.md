# What the pressure-force-off continuation actually isolates

The accepted independent-site-jitter W=0 CUDA state at exact `t=6.5`
(SHA-256 `d02e15ebd01193f1da002d87d689f83525b0150541a9f78dedbc451702c5a2ff`)
was used for both the native `t=6.5→10` pressure-force-off diagnostic and an
independent high-accuracy point-force orbit propagation. This comparison
starts from the **actual perturbed positions and velocities at t=6.5**, not
the original circular-orbit IC.

The [orbit receipt](KEPLER_Q2_JITTER_INHERITED_ORBIT_T10_2026-10-27.json)
uses DOP853 with `rtol=1e-12`, `atol=1e-14`, and `max_step=0.05` for all
55,351 sites. It records the source and endpoint hashes and a tighter-step
sample cross-check. Its disk mass-weighted mean position differences at
exact `t=10` are:

| Comparison | Mean position difference |
|---|---:|
| Point-force propagation vs circular orbit from the original IC | `1.5248222771e-5` |
| Native pressure-force-off vs point-force propagation | `3.1402808220e-12` |
| Native pressure-force-off vs circular orbit from the original IC | `1.5248222776e-5` |
| Standard native hydro vs point-force propagation | `1.6911003181e-6` |

The [independent ordinary Voronoi reconstruction](KEPLER_Q2_PRESSURE_OFF_INHERITED_GEOMETRY_T10_2026-10-27.json)
at the pressure-force-off native positions agrees with native cell volumes
to `1.86e-8` mass-weighted absolute relative error. It predicts the
pressure-force-off native-minus-circular-orbit density residual as
`1.1353606721e-3`, against the native measurement `1.1353607093e-3`.

**Interpretation.** Almost all of the pressure-force-off `t=10` orbit and
density residual is inherited from its already perturbed `t=6.5` restart
state. The point-force propagation reproduces the subsequent pressure-off
site trajectories, so no extra hydro force during `6.5→10` is needed to
explain that branch's residual. Standard hydro does produce an additional
`1.69e-6` mean positional difference from that point-force continuation;
the previous paired density metric falls by `5.12%` when the pressure force
is removed in this interval.

This late intervention **does not identify what created the t=6.5 orbital
perturbation**. In particular it cannot exonerate or blame pressure forces
during `0→6.5`; it cannot establish that pressure is irrelevant to future
growth; and the force-off branch breaks the pressure force/work adjoint,
so it is diagnostic only. The next causal decision should use the accepted
same-state timestep test and the independently running fresh weighted
branch before scheduling any expensive repetition from `t=0`.
