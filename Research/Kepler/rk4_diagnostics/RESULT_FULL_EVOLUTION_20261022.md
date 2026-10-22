# Full native RK4, restart and GPU mass continuation — 2026-10-22

The historical artificial heating and t6–10 breakdown are **not solved yet**.
The actual full-disk trajectory is retained at t=8.2355880701376864e-5 after
four accepted steps. These early steps validate integration and continuation;
they are not a substitute for the original long Kepler target.

## Authoritative results

The original user source remains HEAD `8e25dfe0c2e41a7374c11e2e3102da2a8b028625`,
remote `git@github.com:kjhan0606/LagEunha.git`. Variants use isolated native C
entry, MPI4/OpenMP, CUDA geometry/EOS/forces, fixed binary80 mass, and RK4.
The 71,080 source cells / 1,138,432 triangles / 10,816,441 scalar DOFs and the
original 2,040,256,512 positive quadrature points per pass are unchanged.

1. `material_full_evolution_v2_h200_mpi4`, Slurm432870, native exit0:
   two accepted steps through t=3.8755708565353815e-5, plus one actual second-step
   restart branch. The restarted and uninterrupted state/ledger bitmatch.
   The first stored k1 was reused. No initial mass or curvature was repeated.
2. Its Slurm wrapper FAILED1:0 **after native success** because the source
   manifest mistakenly contained the link-mutated `Exam/libexam.a`.
   `material_full_GPU_mass_v1_h200_mpi4/predecessor_adoption.json` verifies every
   immutable source and identifies this sole expected mutable-archive change.
   This bookkeeping failure does not justify replaying accepted steps.
3. `material_full_GPU_mass_v1_h200_mpi4`, Slurm432872 COMPLETED0:0:
   one genuine next step to t=5.8133562848030723e-5, with a v2 accepted checkpoint
   including endpoint energies, thermal work and immutable generation binding.
   GPU DD local block products/inverses use the original binary80 M for the
   authoritative outer residual. Host vectors, dots and halos remain binary80.
4. `material_full_GPU_mass_audit_v1/contract.json`: independent saved-data audit
   PASS; 69 immutable/input/binary hashes and all 71,080 cell ledgers in each of
   four actual GPU-mass stages checked, without physical query or replay.

| Quantity | Actual evidence |
|---|---:|
| First GPU mass true residual x/y | 8.33e-19 / 8.34e-19 (<1e-16) |
| DD inner iterations per component | 54; one outer correction |
| binary80-to-hi/lo reconstruction error | exactly zero |
| Same-RHS CPU/GPU work difference/U | 3.18e-19 |
| Same-RHS CPU/GPU torque difference/(2K) | 8.77e-23 |
| Third-step finite E defect/U | 1.67e-13 |
| Third-step thermal increment defect/U | -6.13e-19 |
| Third-step angular defect/(2K) | -4.03e-17 |
| Maximum cell rate E defect/U, four stages | 4.16e-13 |
| Maximum cell rate L defect/(2K), four stages | 4.20e-16 |

Resolved spin/inertia/Erot are diagnostics of the actual subcell velocity field.
No independent spin reservoir, damping, steering, heating fix or viscosity was
introduced. Kepler remains Re=infinity. Generating entropy stays immutable.

## Actual cost and the next continuation

Third step: 297.485228 seconds on ONE H200 shared by four ranks, including a
one-time independent CPU mass solve comparison. Subsequent mass solves inside
that step took 19.46–19.63 seconds per RHS; changed-state geometry/EOS/force
preparation took about31 seconds, including about27.6 seconds of CUDA kernels.
The measured material dt bound is ~4.322978e-4. With the existing0.8 safety
factor, even ~225sec/step would require ~7.5days to t1; a blind24h run is not
warranted. The current small accepted h is allowed to grow1.25x toward the
physical bound; it is not fixed to a small test interval.

The next implementation uses a **bounded stage-local GPU point cache** for
complete source-cell tiles (1,048,576 points/rank,568bytes/point), reusing the
same shapes/current fields/generating j0 across dependent passes. It does not
lower quadrature order or store the full two-billion-point array. A tile of the
first genuine next RHS compares original streaming and cached reduced results.
The next accepted endpoint records cell density, area-mean pressure, U/K/L,
intrinsic rotation, positivity and a real restart. Its measured cost determines
the next performance work. Global map injectivity monitoring remains unfinished;
positive sampled J alone is not a proof of no folding.

`material_full_cached_continuation_v1_h200_mpi4` compiled successfully. Submitted
H200433103 with `--gres=gpu:H200:1`, MPI4/OMP3/64G, but it stayed pending for
QOSGrpGRES and was **CANCELLED before allocation or physics**. It must not be
restarted while the Ada alternative is executing.

`material_full_cached_continuation_v1_ada_mpi4` is the manual authorized Ada
alternative, MPI4/OMP8, sm89, restoring the actual third-step checkpoint.
Check its live runner/native results; a launch is not an accepted step.
Ada timings cannot establish H200 cache throughput.

### Completed Ada cached continuation and adopted next bundle

The Ada native run finished exit0 (runner83131 is CLOSED). It accepted step4
with h=2.4222317853346135e-5, reusing the actual third-step checkpoint. No
successful interval was repeated. Actual step442.167919sec; the cost gate
stopped at this accepted checkpoint instead of parking a14.8day t1 run
(optimistic estimate using current measured maximum safe bound).

- `material_full_cached_continuation_v1_ada_mpi4/run/full_continuous_result.json`.
- Finite step E defect/U8.03e-15, thermal defect/U6.13e-19,
  L defect/(2K)-7.90e-19; source mass unchanged.
- The first actual complete-source cache tile/rank agrees with original
  streaming reductions exactly. Ada CUDA pass times ~41/4/9sec; this is not an
  H200 benchmark. Resident GPU memory was about18.2GiB for all four processes.
- `material_full_cached_audit_v1/contract.json` independently checks73hashes,
  all71080cells in4stages and the accepted endpoint density/P/rotation output.
  At this EARLY time, R<.5 maximum cell U/U0 change4.93e-13 and area-mean
  P/P0 change1.24e-12. These are not long-time heating results.
- Latest full checkpoint prefix:
  `material_full_cached_continuation_v1_ada_mpi4/run/continuous_step4/accepted_r`.

`ASTRA_LONG_EVOLUTION_COST_DECISION_20261022.md` is the one strategic evaluation
after the actual cost result. It rejects blind long continuation at present
cost and chooses positive reference cubature compression plus device-resident
DD mass iteration. M, x/v, source grouping and generating entropy remain fixed;
rule discrepancies must have a separate ledger. The production gate requires
the actual next RK step compared at its same stage states with the full rule,
then a measured cost sufficient for t1/24h or10x improvement. No full-rule/order/
dt/tolerance sweep or old-state mass/curvature replay is authorized by that plan.

Native `material_cubature_constraints_v1.cpp` began the fixed representative
preprocessing on the latest accepted state, not t0. It exported24triangles and
93184positive original triangle points plus whole-source delta contexts in
4.86199119sec. Constraints include mass products, determinant derivative span,
actual anchor heat/gradient and softened gravity/potential. No weights/state/M
were changed. All24payload lengths and the exported total were verified.

The standalone preprocessor then exited1 because its Ownership communicator
destructor ran after MPI_Finalize. This is preserved in
`material_cubature_feasibility_v1_mpi4/export_adoption.json`; completed durable
matrices are adopted **only for saved matrix analysis**. Root source now
destroys these objects before MPI_Finalize. Frozenv1 is unchanged and queries
are not replayed for this shutdown mistake. No positive reduced rule has been
constructed or committed yet. Saved-matrix analysis40090 finished exit0:
`material_cubature_matrix_analysis_v1/contract.json` verifies the24payloads'
binary80 target sums (max relative-to-absolute-sum discrepancy2.57e-18).
Each has324--400constraint rows; double SVD diagnostically resolves118--124
directions, but this is not a certified mathematical rank and no direction was
discarded. Positive elimination/finite-arithmetic contractions and resident DD
CG remain to implement. Original source and binary untouched.

## Remaining physical scope

Continue this retained trajectory, not a new IC. t1,6,8,10, original/paper-matched
boundaries, the longest correctly normalized GIZMO orbit count, density/pressure/
orbital accuracy and long conservation remain unproved. The current research
mode has moving material boundaries under constant Pext=P0, unlike the original
periodic baseline. Map injectivity, conservative repartition with compatible
recipient thermal context, W=0 dynamic Voronoi/fixed-mass reduction, and shock
entropy remain incomplete. Fixed entropy alone does not prove absence of
spurious compression/heating.
