# KEP stress filter, original 1D registry, and physical velocity query

## Status at publication

The native C/OpenMP candidate remains experimental. The filtered version3
passes the registered cold-shear pressure gate and improves the preceding
KEP Gresho error. The version4 physical point query corrects an identified
outside-cell extrapolation error, but neither strong-shock trajectory reaches
the original endpoint. No production CUDA/MPI backend, full 2D suite, new
Kepler trajectory, or matched GIZMO speed gate is qualified by the results.
The existing manuscript and accepted Kepler baseline retain their previous
scope. No new experimental result is promoted into the manuscript.

The source and evidence accompany the branch
`codex/native20pi-contact-audit-20261104`. Run directories are under
`/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1`. The lossless archive is
`Research/KEPFilterPoint_20261105`; executable hashes are retained while
compiled executables are excluded. All work and review builds use GPFS.

## Fixed LSQ compression/residual filter

The filter uses the unlimited velocity gradient G from the existing centroid
least-squares objective before slope limiting. The neighbor weight remains
A/|r|^(d+1). Define Eparallel=sum omega |G r|^2 and
Eperpendicular=sum omega |delta_u-G r|^2. The compression fraction is
c=min(1,max(0,-tr G)^2/||G||F^2), with c=0 for a zero gradient. The stress
fraction is chi=(Eperpendicular+c Eparallel)/(Eperpendicular+Eparallel),
with chi=0 for a zero denominator. Each pair multiplies the preceding
compression-only stress by max(chi_i,chi_j). The original linear/quadratic
coefficients remain 1 and 2. No problem-dependent fit or floor is introduced.

The same definition applies in 1D/2D/3D. Fifteen native fixture checks cover
translation, affine shear, rotation, traceless strain, planar compression,
and alternating velocity where applicable. The fixtures pass. Affine
incompressible flow has a roundoff-level fraction, while planar compression
and alternating velocity retain unit fraction. Pair kinetic work and
semidiscrete entropy budgets pass separately from finite-step accuracy.

| Version3 trajectory | Original target | Measured result | Status |
|---|---:|---|---|
| alternating N16 | one acoustic crossing, 0.052822140920532286 | final/initial kinetic energy 0.0005131344 | pass |
| Sod N200 | 0.2 | full endpoint | pass |
| cold sine shear 64^2 | 0.25 | pressure relative L1 0.0001745527682 (0.0174553%) | pass against 1% |
| Gresho 64^2 | 1 | physical node velocity L1 0.02858829788 | improves preceding KEP 0.03691193183 |
| blast N200 | 0.038 | stops at 0.0008291804393 | fail |
| original collision N200 in physical region | 0.035 | stops at 0.008865373711 | fail |

The Gresho result remains less accurate than the separate HLLC result
0.0126674. The native Gresho evolution takes 42.3852 seconds on the recorded
8-thread resources, with 2530 steps and 1653 retries. The measurement is not
a matched GIZMO speed comparison. Cold shear has only two accepted steps at
the registered endpoint. Higher resolution and time-step convergence remain
necessary before adopting the filter.

The first controller failed while serializing a NumPy boolean. Completed
alternating output was assessed without repeating the PDE. The initial
blast/control checkpoints were continued after a numerical solver fix.
The original controller error, corrected assessments, initial failures, and
continuations remain archived rather than replaced by successful-only logs.

## True residual termination and numerical scope

The initial PCG implementation could stop on a recursive relative residual
while the actual scaled flux residual still exceeded 1e-12. The early blast
failure occurred after 21 iterations, before the registered cap. The common
solver now checks the actual residual at a candidate solution and replaces
the recursive residual when necessary. The system, immutable birth mu,
64-iteration cap, and both 1e-12 convergence criteria remain unchanged.
Actual failed iteration counts and reasons replace misleading zero counts.
The correction alone does not complete the strong-shock tests.

## Original 1D initial conditions and moving boundaries

The native registry now connects Sod, blast, collision, Shu--Osher, Noh,
Lax, double rarefaction, stationary contact, and moving contact to the
original setup parameters. Collision uses the original high-pressure states,
not the separate Mach-8 contrast. Shu--Osher uses sin(5*pi*x). Noh uses
gamma=5/3. An independently compiled receipt reads the original setup.
Independent split-Gauss integration checks native mass, momentum, energy,
and birth mu on actual intervals. All nine IC-only checks pass. IC-only
checks perform zero evolution and must not be counted as completed tests.

Open boundaries extend the same infinite-line IC into guard cells. Each
transmissive exterior plane follows 1.5*q_edge-0.5*q_inner, with the same
linear combination of physical velocities. Reflecting walls retain their
original positions. Cells are never inserted, deleted, or remapped. Boundary
fluxes enter the global conservation, kinetic, and entropy ledgers. Errors
are measured on the original physical interval.

Four contact boundary trajectories reach t=0.25 with guard multipliers 1
and 2. Both padding comparisons pass the registered 1e-9 ROI criterion.
The screened weight solve is elliptic. A sound-speed guard estimate is not
a causal independence proof. Shock padding independence and multidimensional
moving open boundaries remain unqualified.

## Physical velocity at an actual point

The saved version3 blast failure has one generator outside its own weighted
cell. The closest generators 21 and 22 both belong to physical cell 22.
Extrapolating each generator's own cell gives gap velocity -1.429554658.
Evaluating both points in cell 22 gives -2.614453615e-5, matching the same-cell
affine identity. The latter result is a read-only causal diagnostic, not a
new evolved trajectory or a guarantee against all generator collisions.

Version4 locates the power cell containing each generator and evaluates
the existing limited physical velocity reconstruction in the containing
cell. A generator inside its own cell retains its preceding velocity value
exactly. Conserved state, gradients, stress, and immutable birth parameters
remain assigned to their original cells. No positional steering is added.
All gradients are available before the separate point-query loop.

The native locator descends actual relative power halfspaces through mesh
adjacency, handles periodic images, and resolves face ties by stable ID.
Cycles, inconsistent geometry, or boundary failure return an explicit error.
There is no fallback to outside-cell extrapolation. In 2072 native 1D/2D/3D
queries the locator matches exhaustive diagnostic ownership; the maximum
walk is 10 hops. Exhaustive search is used only in the audit.

Fresh version4 trajectories use the original full targets and start at birth.
Blast reaches t=0.005558214436 of 0.038 before exceeding eight step halvings.
The exact retry cause remains to be instrumented; the result is not labelled
as a demonstrated owner-transition instability. Collision reaches accepted
t=0.008613241478 of 0.035 before the next Euler trial exceeds PCG's cap.
The collision trial at t=0.008632956056 is saved separately with a diagnostic
version flag. Production restart loading rejects the flagged trial.
Accepted final.fv remains distinct from the rejected failure_trial.fv.

## Conditional block preconditioner and pre-push validation

The optional GFS_NATIVE_BLOCK_JACOBI build selects disjoint pairs with mutual
strongest aggregate edge coefficients. Remaining rows use diagonal Jacobi.
Each pair inverts the principal block [[a+bi,-a],[-a,a+bj]]. Positive bi and
bj are summed directly from mu and unpaired edges. The determinant uses
a*(bi+bj)+bi*bj with scaling, avoiding subtraction of nearly equal squares.
The preconditioner changes numerical solution cost only. The physical
checkpoint version remains 4.

Before any block-backend PDE launch, both preconditioners were evaluated on
the exact rejected collision trial. Matrix coefficients, RHS, scaling, and
cell diagnostic bytes are identical. All mu are positive and edge coefficients
are nonnegative. Diagonal Jacobi stops at 64 iterations with recursive
relative residual 1.018402736e-12 and true scaled residual 1.025337318e-13.
Block Jacobi converges in 49 iterations with recursive relative residual
6.567553298e-13 and true scaled residual 4.203537373e-14. The registered
cap and tolerances are unchanged. The conditional read-only gate passes.
No full block-backend trajectory has yet been performed.

The pre-push review builds default, filtered, global-query, and block native
executables plus five audit executables. All nine builds succeed. The 2072
point queries, 15 stress fixtures, same-cell affine identity, unchanged-state
check, and exact-trial comparison pass. Compiler warnings remain in the
preserved Voro source and one audit indentation. Python source syntax is
checked separately. Full strong-shock endpoints remain the next physics gate.

## Next work

1. Reproduce the failed blast step from its accepted version4 checkpoint with
   passive retry-reason logging. Do not change caps or tolerances.
2. Apply the common block numerical backend to fresh full-target original
   blast/collision trajectories after the passed exact-trial gate. Do not
   select preconditioners by test name or transfer version3 results to version4.
3. After strong endpoints pass, repeat contact/Sod/cold/Gresho for the common
   definition, add strong-shock refinement and guard-doubling checks, and
   complete the original nine tests before expanding to the full 2D/3D suite.
4. Qualify GPU/MPI evolution and matched GIZMO timing separately. Retain the
   existing accepted Kepler history and manuscript baseline throughout.
