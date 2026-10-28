# Generator velocity recovery changes the interpretation, not the remaining mixing

Analyze saved final N4/N8 hybrid paths without further remapping. Recover
v(X)=p/m+G(X-q/m), using G=C^T S^-1 from the same transported moments.
Compare with original t2 velocities; geometry returns to its original state.
No generator is moved and no steering or force is introduced.

| N8 region | Mean velocity change mass RMS | Recovered generator velocity change mass RMS |
|---|---:|---:|
| Whole r<.7 patch | .000423211 | .000253719 |
| Inner r<.5 | .00211872 | .00125494 |
| Edge .5<=r<.6 | .000157453 | .000109593 |

For the whole patch, the point reconstruction reduces the displayed velocity
change by about40% relative to the cell-mean change. These are different
representations of the same saved state, not a40% improvement in a physical
trajectory or a new numerical method. N4 gives nearly identical values.
Original constant-in-cell donor velocity is the unprojected reference; the
remaining change measures this diagnostic projection's non-reversibility.

N8 whole-patch kinetic split relative to original thermal energy:
resolved affine .00459339, residual .0359725, total .0405659. Thus about89%
of the internal kinetic content remains residual. In the inner region those
ratios are .00885618 and .0671078. Relative global thermal change is only
1.53e-13; tiny regional changes are redistribution, not global heating.
The unchanged outer band has only roundoff-scale effects.

## Consequence

Mean velocity change must not be equated to generator trajectory error in an
extended-cell representation. However point recovery does not eliminate the
remaining velocity change or determine how residual covariance couples to
physical pressure and face transport. The preceding conclusion that the
hybrid reduces internal kinetic mixing only about5% remains valid for that
specific energy measure; it is not a universal accuracy score.

Do not add physical evolution or production CUDA kernels for this closure
until the residual transport and approved original-model reduction are stated
consistently. Next compare the representation on a spatially refined controlled
shear remap with fixed smooth velocity field, retaining the original Kepler
geometry/data tests as the relevance check. This distinguishes spatial
projection error from the time-subdivision plateau already measured.

Evidence: `canonical_hybrid_generator_recovery_2026-10-18.json`.
Script: `audit_hybrid_generator_recovery.py`.
Remote: existing hybrid result directory, `generator_recovery.json`.
No long-time stability or full backend/MPI claim.
