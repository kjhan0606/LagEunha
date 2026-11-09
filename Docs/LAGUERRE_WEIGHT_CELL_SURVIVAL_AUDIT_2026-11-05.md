# Laguerre weight and cell survival audit

## Answer to the weight-bound proposal

Weights are not being rejected. The current source contains several different
weight mechanisms, and the proposed neighbor-distance rule is already present
in one experimental path. It has not been carried consistently into the newer
native `W` state, and the legacy face-based ceiling is not the same rule.

For the power distance used by the CUDA tessellator,

\[
\pi_i(x)=|x-x_i|^2-W_i,
\qquad
C_i=\{x:\pi_i(x)\leq\pi_j(x)\ \forall j\}.
\]

`W` has units of length squared. The bisector fraction along the segment
`x_i` to `x_j` is

\[
f_{ij}=\frac12+\frac{W_i-W_j}{2|x_i-x_j|^2}.
\]

The symmetric strict condition

\[
|W_i-W_j|<|x_i-x_j|^2
\]

for every pair (including periodic images, when present) keeps each generator
strictly inside its own power cell. A generator then has a neighborhood of
positive volume in that cell. This is a sufficient cell-survival condition,
not a necessary condition: a power cell can have positive volume even when its
generator is outside it. It is relative-weight invariant: adding the same
constant to all `W_i` changes no cell. A bound on absolute `W_i` is therefore
not the geometric condition.

## What the current code does

1. **Legacy adaptive `w2`:** `applyW2Controls` applies a floor, optional
   relaxation, optional rate limiter, then `w2 <= w2ceil`. The tessellator
   sets `w2ceil_i = min_j(|x_i-x_j|^2+w2_j)` over its current real face
   neighbors. This is a one-sided, current-face ceiling, not the global
   symmetric pairwise rule. It may be conservative or miss a competitor when
   the topology changes.
2. **Experimental `GFS_DYN_W2`:** the source already applies the idea in a
   more explicit form on its current real faces. It limits the step using
   `|w2_i-w2_j| < (1-delta)d_ij^2` (`delta=0.2` by default), then activates a
   weighted barrier target after the ratio exceeds `1-2 delta`. This path also
   has a PCG solve and tunable barrier strength. It is enabled only for AV mode
   5, `kappa=0`, one MPI rank, and CPU execution; it is disabled on CUDA. It is
   an existing experimental model, not a qualified production path.
3. **Native persistent `weight.W`:** the newer state has units `L^2` and
   evolves by RK4 from `dW/dt=pi/inertia`; the RK stages do not apply a neighbor
   cap. CUDA computes a `w2ceil` diagnostic from current faces but does not
   apply it to `weight.W`. The geometry kernel clamps a raw bisector fraction
   into `[0.05,0.95]` and counts any clamp as a geometry error. For each
   candidate pair it visits, that is effectively the stricter test
   `|W_i-W_j| <= 0.9 d_ij^2`. The bounded native weighted validation can reject
   and retry a whole trial through its opt-in wall-event transaction; ordinary
   runs without that transaction abort on the geometry error. The clamped
   fraction is not accepted as a replacement power diagram.

The native `weight.W` state also carries momentum and kinetic energy
`pi^2/(2 inertia)`. Its restricted validation path can roll back a failed
geometry trial, but the rejection currently comes from fractional-face
geometry failure rather than a separately measured all-pairs margin. Mutating
`W` after an RK step while leaving `pi` unchanged would alter the dynamics
without accounting for the constraint work. A silent post-step clip is not
acceptable; the first explicit candidate should be a whole-step geometric
admissibility check and time-step reduction, with any later active constraint
force derived and energy-audited.

## Existing evidence and limits

The exact power-geometry audit at
`Research/FVNative_20261105/fv_power_geometry_gate_v2/stdout.json` passed its
volume and gauge checks with raw face fractions spanning `-0.1` to `1.1`.
The frozen source gives a concrete periodic example: four sites per axis,
spacing `h=1/4`, weights `[0, 1.2 h^2, 1.2 h^2, 0]`, and exact cell-edge
coordinates `[0, 0.1, 0.5, 0.9, 1]`. Every cell has positive width although
some active-face ratios are `|W_i-W_j|/d_ij^2=1.2`. The symmetric cap therefore
excludes an exact, positive-volume Laguerre partition. This demonstrates that
the cap is a stronger modeling invariant, not a necessary condition for
positive cell volume. It may still be selected to keep generators inside
their cells, if that is the intended model.

The jittered 2D material-polygon pressure-pulse fixture stopped at
`t=0.01772452294471101` with **all Laguerre weights equal to zero**. That
convexity failure cannot be repaired by a weight cap. N=32 fixture runs were
interrupted for runtime and produced no physical result. The fixture is not
the production MPI/CUDA solver.

## Next action

First run a geometry-only comparison on retained Laguerre configurations:
measure the exact cell volume and the all-pairs margin
`|x_i-x_j|^2-|W_i-W_j|`, including periodic images. Compare that certificate
with the current-face `w2ceil` and `GFS_DYN_W2` checks, especially at topology
changes. Do not change production `W` values in this audit. If the symmetric
certificate is selected, integrate it as a trial-stage validity condition and
test rollback/time-step behavior before considering a physical constraint law.
If the trajectory stalls at the admissible boundary, report that explicitly;
do not hide it with a clamp. Then continue the already-authorized 2D KH/RT
program only with the resulting rule and test status stated clearly.
