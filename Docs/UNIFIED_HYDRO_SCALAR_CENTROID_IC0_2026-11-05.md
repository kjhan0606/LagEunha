# Scalar velocity limiting, physical strain length, and IC(0)

The common physical candidate is checkpoint version8. The candidate adds
componentwise scalar velocity bounds at centroid midpoints and uses the
positive centroid-normal separation in the artificial-viscosity strain.
Reflecting normal rows retain the constrained wall field. The raw constrained
least-squares gradient still defines the stress sensor. Face kinematics and
the weight-rate matrix retain the true generator separation.

The numerical backend optionally uses unshifted, zero-fill IC(0)-LDL^T in
native row order. IC(0) changes the preconditioner alone. The physical
version, A, b, birth mu, 64-iteration cap, and actual scaled flux-residual
tolerance 1e-12 remain unchanged. Duplicate owned-pair edges are aggregated
only for factor construction. Invalid pivots cause explicit failure.
Triangular sweeps are serial and rebuild the factor on every solve.

The lossless evidence archive is `Research/ScalarCentroidIC0_20261105`.
The original files remain under
`/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1`. The archive retains frozen
source, build commands and hashes, logs, matrices, accepted checkpoints,
and non-restartable rejected trials. Compiled executables and Python caches
are excluded. All fluid evolution uses native C/OpenMP on GPFS.

## Physical operator checks

All 52 native fixtures pass, comprising 22 wall and corner fixtures,
15 frame/contact fixtures, and 15 stress fixtures across 1D/2D/3D.
The checks cover affine consistency, boosts, rank-one aliases, wall traces,
and kinetic/entropy balances. Nonlinear reconstruction and frame errors
decrease on shrinking smooth patches. Exact nonlinear rotational covariance
is not established. The local patch checks are not PDE convergence tests.

Frozen 2x2 controls isolate the two physical changes on the same historical
blast state. Scalar limiting changes the nearest-pair generator gap rate
from -0.82847104 to -2.32509745e-9. Centroid strain length changes the pair
artificial pressure from 0.64298890 to 1.29617235e-8. The controls preserve
A exactly. Only the combined physical candidate is evolved.

## Numerical diagnosis and preconditioner checks

The fresh version8 block-PCG runs stop at t=0.01083643632 for blast and
t=0.02059531253 for collision. Neither exact failed system has an earlier
true-residual pass within 64 iterations. Independent 60-digit tridiagonal
LDL solutions, rounded to double and evaluated by the native flux routine,
give actual residuals 8.94e-16 and 5.09e-16. A more accurate root therefore
exists for both saved systems.

The grouped numerical evaluation selects IC(0) without shifts, fill,
reordering, factor reuse, tolerance changes, or physical modifications.
Four frozen systems pass A/b/scale equivalence, native matrix action,
preconditioner symmetry, positive-pivot, and cost checks. A synthetic graph
also checks duplicate-edge aggregation and explicit invalid-pivot rejection.

| Saved system | Block iterations | IC(0) iterations | IC(0) actual residual | IC(0)/block solve cost |
|---|---:|---:|---:|---:|
| blast N200 failed trial | 64, failed | 1 | 3.13e-15 | 0.104 |
| collision N200 failed trial | 64, failed | 1 | 5.09e-16 | 0.118 |
| Gresho 64^2 historical state | 26 | 13 | 8.24e-14 | 0.645 |
| Sod 16^3 historical state | 36 | 15 | 1.02e-14 | 0.580 |

Costs include assembly, factorization, triangular application, and true
residual checks. Geometry and I/O are excluded. The 1D denominators are
failed block64 work windows, not equal-accuracy speed comparisons. The
historical 2D and 3D states supply frozen operator checks, not new version8
trajectories. The measurements do not establish the GIZMO <=2x goal or
production CUDA/MPI performance.

## Accepted-checkpoint continuations

Both continuations load the last accepted version8 checkpoints. Entry A, b,
scale, and physical generator velocities match exactly across backends.
The copied initial checkpoint bytes equal the parent accepted final bytes.
Production restart rejects diagnostic trial files. No calculation restarts
from birth merely to change the numerical backend.

| Test | Original target | Last accepted t | Total accepted steps | Endpoint status |
|---|---:|---:|---:|---|
| collision, physical N200 plus guards | 0.035 | 0.035 | 2149 | reached target |
| blast N200 | 0.038 | 0.03501883433 | 3363 | failed before target |

The collision continuation adds 888 accepted steps. The whole physical
trajectory has maximum conservation error 9.64e-16 and residual-corrected
mass-weight invariant error 1.70e-12. Endpoint completion alone does not
establish agreement with the exact Riemann solution. The physical error norm
has not yet been evaluated for the completed collision state.

The blast continuation adds 2318 accepted steps. The next rejected Euler
trial has dt=2.794999995e-5 and actual residual 1.576629527e-12 at cap64.
The positive factorization does not fail. The last accepted state has a
minimum neighboring generator gap 5.493666633e-10. A causal explanation for
the late residual failure remains unresolved. The tolerance remains 1e-12.
The whole physical trajectory has maximum conservation error 2.49e-14 and
residual-corrected mass-weight invariant error 5.04e-13.

Both physical trajectories combine a block-PCG prefix and an IC(0)-PCG
continuation under the same version8 equations. The combined trajectories
are not all-IC(0) timing benchmarks. The blast failure closes the conditional
strong-test bundle. N400 refinement, the full original-nine suite, cold64,
and subsequent complete 2D/3D runs have not been released by the failed gate.

## Publication and remaining scope

The native C sources match the final frozen IC(0) build. The two updated
Python runners match the copies actually used for the accepted-checkpoint
continuations. The older build manifest keeps the older runner hashes
rather than being rewritten after the runs.

The accepted Kepler manuscript remains unchanged. Overleaf main is
6cb0e1a471aa46aad65496a30ec23c5da5158dfa and matches the manuscript snapshot
retained in GitHub. The present experimental results are not substituted
for the accepted Kepler trajectory or historical test configurations.

The full goal remains open. Required work includes the late blast diagnosis,
exact-solution physical errors and refinement, all original 1D tests, the
same-method full 2D suite, consistent 3D evolution, production GPU/MPI
integration, and matched GIZMO runtime <=2x.
