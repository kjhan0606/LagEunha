# Irregular 2D cells: one-remap mixing decreases with spatial resolution

Perturbed N-by-N generator grids on [-1,1]^2, N16,32,64, jitter within .15h,
fixed random seed. Cells use the same diagnostic Laguerre polygon construction
and [-3,3]^2 clipping domain. Measured annuli r=.2-.5 and .5-.7 are interior.
Use the original canonical spline-softened circular velocity coefficient,
piecewise-constant donor velocity and uniform density. There is no disk density
jump and these are not replacement Kepler initial conditions.

Apply one small compact-support weight field proportional to h^2; this holds
deformation in cell units comparable while refining, not absolute weights.
Direct overlap integration retains m,q,Q,p,D,Pi,U. No physical time evolution.

| N | Inner internal K / mass | Inner residual K / mass | Inner recovered point velocity change RMS |
|---|---:|---:|---:|
| 16 | .00151949 | .00145807 | .00462641 |
| 32 | .000368088 | .000353480 | .00249424 |
| 64 | .0000921209 | .0000883109 | .00115320 |

From N32 to64, internal/residual energy exponents are1.998/2.001 in the inner
annulus and1.969/1.970 outside. Point-velocity changes have exponents1.113/.877.
Three nonnested grids do not prove a general asymptotic scheme order. This
supports spatial decay in a curved-flow irregular-cell case, beyond the
symmetry-dominated strip check. Coverage error<=6.80e-14; maximum global raw
moment drift<=1.43e-13. These are absolute code units, not solver tolerances.

## Point recovery qualification

At N64 the inner cell-mean velocity change RMS is .00111739, smaller than the
recovered point change .00115320; the outer values are .000190065 and .000205603.
The earlier40% reduction on the saved hybrid round trip is not a universal
benefit of affine point recovery. Point recovery changes representation and
must be checked against the field intended by each test. No physical trajectory
improvement follows from either comparison alone.

## Remaining relevance gap

This test is one remap, uniform density, newly perturbed grids, and fixed
generator geometry. It does not cover the canonical density contrast, repeated
physical pressure/AV feedback, moving-generator geometry, or long-time Kepler
stability. Do not adopt the experimental moment closure based on this result.
Next establish how spatial representation and residual kinetic bookkeeping
affect pressure at the canonical density interface before extending physical
evolution. Preserve original point-state energy and approved reduction.

Evidence: `canonical_irregular_circular_refinement_2026-10-18.json`.
Script: `audit_irregular_circular_remap_refinement.py`.
Remote: `/gpfs/kjhan/LagForce/Codex_irregular_circular_refinement_20261018`.
No production code, checkpoint, CUDA/OpenMP/MPI execution path or job changed.
