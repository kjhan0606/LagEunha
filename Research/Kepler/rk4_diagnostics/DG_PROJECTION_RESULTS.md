# Conservative Euler projection: P2, P4, P6 on the same native Kepler mesh

Subsequent representation implementation/evidence: `ENRICHED_EULER_RESULTS.md`.
This document retains the failed pure polynomial projections and their cause
analysis; no new time-evolution claim supersedes them.

2026-10-20–21. The preceding goal turn made progress: native affine initialization
and restart were implemented, and a proposed non-Euler kinetic reservoir was
rejected. This bundle implements the next selected Euler spatial projection
and obtains evidence that changes the representation decision. No time
evolution of this new candidate has occurred; t6–10 remains unverified.

## Native implementation and scope

Canonical repository remains
`/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src`,
origin `git@github.com:kjhan0606/LagEunha.git`, HEAD `8e25dfe0c2e41a7374c11e2e3102da2a8b028625`
plus uncommitted driver changes. No old code, input, checkpoint or output was
deleted/replaced, and no push/shared-context update occurred.

New `Exam/gfs_dg_euler.h` implements a shared mean-free Euler polynomial basis,
positive Gauss/Duffy integration, one Gram factorization for rho,px,py,E,
standard Euler pressure/flux and canonical-face/fixed-basis weak-form helpers.
Cell-integral modes are preserved structurally. `gfs_dg_gpu.inc` runs projection
and independent point assessments in CUDA on the existing native polygons.
`gfs_dg_native.inc` exposes opt-in `GFS_DG_PROJECT_ONLY=1`. Degree2 is first;
degree4 and6 are bounded spatial-order controls after degree2 fails.

This mode exports projection/geometry packets and exits. It writes no evolved
checkpoint, applies no positivity limiter, energy repair, mass/weight change,
generator steering, or independent Pi evolution. There is **no new DG time
integrator, moving-basis stage update, weight solve or shock-validation pass**.
The existing native C/OpenMP+MPI+CUDA body is retained.

## Checks completed

Frozen build `dg_projection_v2_build`, source hashes and binary SHA retained.
`contract.json` verifies orders2,4,6: Gauss polynomial-moment error5.56e-17,
constant-pressure recovery2.27e-13 for exact representable fields, conserved
integral error3.56e-15, instantaneous fixed-basis affine-mesh Reynolds/weak-form
identity7.82e-14. This does not validate a moving-basis timestep.

Native4096-cell rigid rotation, boosted shear and varying-density/constant-P
controls pass. For P2 the maximum cold boosted pressure relative error is
9.92e-8 and varying-density isobaric error4.49e-10. The latter also confirms
that the rejected constant-specific-u assumption is not being used.

The P2 actual71080-cell projection has bitwise-identical IDs, geometry, means,
all Euler coefficients, conserved integrals, pressure diagnostics and flags
between MPI1 and MPI4. This is projection consistency, **not a DG-state
migration/restart implementation**.

First v1 actual projection stopped on63 fake signed measures from exactly
inactive edges. With CUDA contraction, a determinant of identical vectors can
retain one multiplication's roundoff residual. v2 explicitly gives identical
endpoints zero measure; it does not remove a resolved edge or change a weight.
The successful full v2 projections replace no previous artifact. All71080
geometry/linear solves complete for every degree.

## Same-IC spatial result: no candidate passes

All use the same original continuous rho,v,P, epsilon=.047,71080 generators
and quadrature family. Pressure is nominally1e-6. No limiter is applied to
make these results positive.

| Degree | Cells with nonpositive sampled P | Global volume L1 relative P error | Largest conservation-relative error |
|---|---:|---:|---:|
| P2 | 1552 | .01336196 | 1.39e-15 |
| P4 | 4 | 2.71955e-5 | 2.60e-15 |
| P6 | 4 | 1.18863e-6 | 2.51e-15 |

The four remaining cells are67943,68000,68527,68584, at initial radius
.00828640759203, with the physical origin as a polygon corner. P6 elsewhere
improves substantially: max relative P error .00021533 for .047<r<.1,
4.8845e-6 for .1<r<.2; the old focus ID68526 has about.00482076. Positive
pressure alone still would not certify a long-time physical solution.

Manual lageunha/Ada wall times for71080-cell native projection are P2 2.634s,
P4 2.751s, P6 3.204s; P2 MPI4 3.889s. These include process startup and native
geometry. Higher-degree accuracy checks did not require long simulations.

Authoritative `dg_projection_v2_p{2,4,6}/{manifest,analysis}.json` retain all
radial/point data and the original input SHA. Earlier controls are in
`dg_projection_v1_p2`; corresponding v2 zero-edge handling leaves their
nondegenerate geometry unchanged.

## Independent integration distinguishes order error from quadrature error

`reference_dg_projection_v1.py` uses NumPy/SciPy-free Gauss integration on the
**exported native polygons**, a local Galilean frame and an independent linear
solve. It projects the same Euler field; no bulk kinetic energy is removed.
Selected16/32/64-point results are retained in each `reference.json`.
Away from the central corner cusp, pressure predictions converge to the
native values: the P2 negative pressure at old focus68526 remains
-1.6528445922e-5 with independent integration. It is not roundoff.

Naive Gauss refinement at the four origin-touching cells converges slowly.
`analyze_dg_core_profile_v1.py` instead anchors the integration at the physical
origin and segments its radial coordinate geometrically around the literal
1e-6 density denominator. The three positive quadrature controls converge
in pressure to about1e-9 relative or better. It also integrates the same
velocity with a **separately labeled constant-density analysis control** of
equal cell mass. Neither control rewrites the actual IC or runs a simulation.

## The existing IC contains a separate narrow central density cusp

`Exam/KH/util.c` defines, for this original smooth branch,

    rho=.001 + [2/(r+1e-6)] * .5[1+tanh((r-2)/.25)]
                                * .5[1+tanh((8-r)/.25)].

The nonzero inner tanh tail multiplied by1/(r+1e-6) gives rho(0)=.2260703,
rho(1e-6)=.1135361, rho(1e-4)=.00323020, while the core background is near.001.
This scale is independent of the gravitational softening epsilon=.047 and is
not resolved by the current minimum spacing .01171875. No IC parameter or
formula is changed in this bundle.

For representative corner cell67943, converged independent pressure error:

| Degree | Original density: max relative P error | Equal-mass uniform-density control |
|---|---:|---:|
| P2 | 45.1533 | 15.6598 |
| P4 | 50.1752 | .349079 |
| P6 | 66.7482 | .00508767 |

P6 original minimum P is -8.69807e-6 in the dense independent assessment;
the control minimum is +9.94912e-7. This isolates an obstacle to **this
polynomial projection** at the four central cells. It does **not** replace the
previous native pressure-work/thermal-depletion diagnosis or prove that this
cusp caused the historical t6–10 failure. Adjacent old-focus errors and
curved-velocity errors also exist without the cusp.

## Angular and compatibility contracts

`dg_angular_v1.json` checks the actual canonical HLLC pair and fixed-basis
spatial moments. Full pair angular residual is zero while generator-orbital
torque alone is -.2217943867; symmetric Euler volume residual is1.39e-17.
This is an operator identity, not evolved DG angular conservation.
`dg_legacy_v2/manifest.json` records a native3-step legacy SSP replay with
all exported state and metadata bitwise identical to its frozen reference.

## Next representation selected by bundled review

Do not start a long run or keep increasing global degree on the strength of
global L1 improvement. Preserve original IC/history and keep the goal intact.
The Euler representation must handle cold pressure and nonpolynomial central
structure while preserving actual conservative means, rotation/angular
moments and proper shock heating. A primitive/energy-compatible nonlinear
representation or disclosed equilibrium enrichment may be needed; a benchmark
definition change, energy deletion, independent kinetic reservoir or forced
orbit is not a substitute. Astra receives this completed bundle once for
the next implementation decision. The completed Astra review selected
`ASTRA_ENRICHED_EULER_PLAN.md`: P2 plus the four reference Euler fields as a
common scalar space, with every coefficient evolving by full Euler equations.
No stationary reset or residual cancellation is allowed. Enriched projection
alone cannot authorize evolution; complete moving-cell RHS, conservative
weight/mass coupling, runtime-state, perturbation and shock gates remain.
