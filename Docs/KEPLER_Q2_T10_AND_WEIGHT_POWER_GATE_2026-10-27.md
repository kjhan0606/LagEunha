# Bounded q2 Kepler: equal-time t=10 result and weight-power gate

## Native t=10 comparison

The original CUDA+MPI+OpenMP aligned continuation **437328** reached exact
`t=10`, exited `0:0`, and passed all 17 native physics gates. The independent
5% site-jitter continuation **437344** also reached exact `t=10` with a finite,
positive state, unchanged cell masses, `W=0`, and no energy floor or wall
crossing. Its auditor deliberately returned failure because the disk's
mass-weighted native-minus-exact-orbit density residual was
`1.1966839e-3`, above the predeclared `1e-3` gate. Its other 16 checks
passed; this remains a **failed physics endpoint**. The read-only paired
analysis **437367** completed `0:0` and explicitly preserves that failure.

| Disk metric at equal t=10 | Aligned | Independent 5% jitter |
|---|---:|---:|
| Mean absolute density change | `1.96891e-4` | `2.37446e-3` |
| Exact circular-site geometry prediction | `1.94290e-4` | `1.47192e-3` |
| Native-minus-ballistic density residual | `6.03158e-6` | `1.19668e-3` |
| RMS radial/circular speed | `2.00775e-8` | `1.29712e-6` |
| Whole-domain relative energy error | `8.10306e-9` | `8.07655e-9` |

The jitter residual is **198.4 times** the aligned residual (versus 119.4
times at `t=6.5`). In the initial `0.5<=r<0.75` band it is `1.59140e-3`
versus `2.89102e-5`, a factor 55.0. The inner jitter signed entropy change
is only `+1.27628e-7`, so this endpoint still does not establish positive
mean irreversible heating. The common-scale [density field](KEPLER_Q2_BOUNDED_ALIGNED_JITTER_DENSITY_T10_2026-10-27.png)
shows the persistent inner disk and the much stronger disordered fine
structure. [Paired metrics](KEPLER_Q2_BOUNDED_ALIGNED_JITTER_PAIR_T10_2026-10-27.json)
and [inner-band tails](KEPLER_Q2_BOUNDED_JITTER005_RADIAL_TAILS_T10_2026-10-27.json)
are tied to the original endpoint hashes. This is one bounded polar IC pair,
not a GIZMO Fig. 7 equivalent result.

## A changing Laguerre weight has a pressure-work cost

For power-cell volumes, the stage derivative is
`Vdot = B v + J Wdot`, where `B` is the generator-position derivative and
each interior face contributes `L/(2d)` to the weight Jacobian `J`.
The existing smooth W0 variational pressure force is `F = B^T P`.
If fixed-mass thermal energy is updated adiabatically from the **actual**
weighted volume, `Udot = -P^T (B v + J Wdot)`, then

`d(K+U)/dt = -P^T J Wdot`.

Thus a solver that merely inserts the validated GPU `Wdot` into the face
velocity and RK4 weight state would generally break total-energy
conservation. The [independent 25-cell periodic audit](KEPLER_WEIGHT_RATE_ENERGY_IDENTITY_2026-10-27.json)
checks `J` against finite-differenced power areas to `7.92e-11` maximum
error. It chooses `Wdot` to cancel every cell's geometric area rate to
`5.55e-16`. The unweighted pressure work cancels to `2.22e-16`, whereas
the weighted fixed-mass kinetic-plus-thermal power is `-0.167117`;
`-P^T J Wdot` predicts it to `5.83e-16`. This is a structural identity,
not evidence that the Kepler run has an energy error of this size.

On the actual accepted W0 `t=6.5` states, the [aligned](KEPLER_Q2_BOUNDED_ALIGNED_WEIGHT_WORK_T65_2026-10-27.json)
and [jittered](KEPLER_Q2_BOUNDED_JITTER_WEIGHT_WORK_T65_2026-10-27.json)
pressure-weight power audits evaluate the validated Jacobi16 rates on
148,557 unique power faces. They find `P^T J Wdot = -3.12e-13` and
`+5.29e-12` respectively; face and cell sums agree to below `1e-25`.
At the single `0.00193983` step scale the unmatched **global** work is
only `1.69e-11` and `2.87e-10` of current thermal energy. That near-zero
value is specific to these cold, nearly uniform-pressure states and
cannot replace a conservation law over long evolution. Absolute face-work
sums exceed the signed totals by about 69 and 65 times respectively,
so local exchanges also matter.

## Fixed-mass candidate: a pressure-work-neutral rate

One inexpensive algebraic screen keeps the 16-step Jacobi geometry solve
and projects its rate onto `P^T J Wdot=0`. With
`h = J_full^T P_full` restricted to the disk's active weights, the
projection is `Wdot_new = Wdot - h (h·Wdot)/(h·h)`. It requires one
global scalar reduction per stage in a parallel implementation. On the
two accepted W0 `t=6.5` states, it changes the rate vector by only
`0.72%` aligned and `0.32%` jittered in relative L2. The pressure-weight
power falls below `1.4e-26` in both cases. An independently rebuilt
one-Euler-step power diagram leaves inner-band p99 relative area change
`3.268e-5` aligned and `3.224e-5` jittered, versus `3.203e-5` and
`3.214e-5` from unconstrained Jacobi16 and `1.300e-4` and `1.993e-4`
with frozen weights. No disk cell is hidden; maximum `|delta W|/d^2`
is below `3.43e-4`. [Aligned](KEPLER_Q2_BOUNDED_ALIGNED_ENERGY_NEUTRAL_RATE_T65_REBUILT_2026-10-27.json),
[jitter](KEPLER_Q2_BOUNDED_JITTER_ENERGY_NEUTRAL_RATE_T65_REBUILT_2026-10-27.json).

This constraint restores the **global first-order energy identity** in
the small mathematical model while largely preserving the one-step
volume benefit. It does not show that moving a fixed-mass power face is a
material motion, that each cell has a valid local first law, or that a
weighted RK4 step preserves momentum and energy at finite timestep.
Those questions remain native release gates. The alternative is a
conservative moving-face flux that transports mass, momentum and energy;
its `W=0` fixed-mass reduction must be checked against the user's
approved condition.

## Consequence for the native experiment

The GPU rate arithmetic and native FaceCSR attachment are validated, but
the next code change cannot be a bare RK4 `W += dt*Wdot` switch. A complete
operator must account for weight work through the pressure-work-neutral
constraint **plus a checked local law**, conservative opposite-face
mass/momentum/energy transport, or a derived weight degree of freedom and
its conjugate energy. It must retain the user-approved `W=0` Voronoi and
fixed-mass reduction, no steering, and `Re=infinity` in Kepler. The
candidate must first pass a weighted small-mesh force/work identity, then
a native one-step pairwise and global conservation audit, positivity,
and GPU stage-time comparison before continuing the same accepted
checkpoint chain. A cell-spin term is a later residual-shear experiment;
it does not balance `P^T J Wdot` by itself.
