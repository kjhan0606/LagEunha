# Compact CUDA solve passes; the fixed compact space fails initial curvature

2026-10-22. The selected core-protected P3 plus two global initial-residual
candidate has now been implemented and evaluated on the full saved domain.
Its linear solve passes, but its spatial accuracy gate fails. **This fixed
candidate is rejected and no physical time evolution was started.** The
t6–10 artificial-heating and longest-GIZMO stability objective remains active
and unachieved. The prior authoritative accepted step5 remains unchanged at
t=0.00011263377801805954.

## Actual new native CUDA/MPI implementation

`material_compact_gpu.cu` reuses the original native CUDA physics, DD
arithmetic and resident vector kernels, with an added compact CSR action and
globally assembled Jacobi diagonal. `material_compact_solver.h` supplies
owned/ghost DD exchanges, owned-only scalar reductions, polynomial solves and
a two-column Schur solve. Vectors stay on the GPU during inner iteration;
only sparse ghost pairs and scalars are exchanged. An already assembled
Jacobi output is not summed across peers a second time.

`solve_compact_initial_20261022.cpp` checks the final full equation against
**T^T M_original T**, using the actual original binary80 element M, fine
owned/ghost assembly and the correct uniquely-owned transpose. It does not
use the new CSR matrix itself as the final full-system oracle. The saved
weak RHS remains T^T(M_original a_saved), inheriting the original initial
solve residual; it is not a fresh force query at a changed state.

All compilation inputs were copied into a separate73-file source snapshot
before the final build. This also prevents quoted includes from selecting
mutable ROOT prototypes. The first host compile found a missing namespace
alias; it was fixed before any execution. That initial build failure was
retained separately. The original physics source files were not modified.

## Linear solve: passed at original precision

Authority:
`material_compact_initial_v1_ada_mpi4/run/compact_initial_solve_result.json`.

| Quantity | Actual result |
| --- | ---: |
| GPU / MPI / OpenMP | One RTX5000 Ada / 4 / 8 threads per rank |
| Compact scalar DOF | 1,581,183 |
| Global partial CSR entries | 26,882,512 |
| True original-M equation residual, x | 1.7727218461e-17 |
| True original-M equation residual, y | 1.8656282136e-17 |
| Full outer corrections, x/y | 1 / 1 |
| Polynomial inner iterations over four solves | 194 |
| Computed Schur diagonal | 0.68656726795733848, 0.68656726795733848 |
| Computed Schur off-diagonal | about1.32e-18, 1.25e-18 |
| Rank0 inner CUDA solve time, four solves | 1.12017077 s |
| Rank0 Schur setup | 0.637947824 s |
| Rank0 full solve including original-M residuals | 2.06909983 s |
| Initial velocity coefficient reconstruction difference | 0 |
| Initial position coefficient reconstruction difference | 1.1102230246e-16 |

The Schur determinant is positive and its independently computed
off-diagonal entries agree within the fixed numerical criterion. This is
a computed linear-algebra result, not an interval proof. Border normalization
was applied consistently to matrix columns, coefficients, loads and the
original-M residual; the physical field was not rescaled.

This establishes a working smaller-system CUDA/DD backend. Since the trial
space and matrix size changed, the two-second measurement is not a measured
speedup of the identical old10.8M-DOF operator. It is also not RK4 throughput
or a t1/t10 cost-gate pass.

## Whole-reference initial curvature: failed

After the solve passed, the same native batch lifted its acceleration into
the retained fine reader representation and invoked the original CUDA
curvature kernel on **all2,040,256,512 original positive points**. No mass
integration or accepted evolution interval was replayed. The original
gamma/P0/eps were read from immutable thermal input; no entropy, pressure,
physical heat or trajectory correction was applied.

Authority: `run/full_initial_curvature_result.json` in the same run.

| Quantity | Actual result / unchanged requirement |
| --- | ---: |
| Maximum point chi | **0.00879341071151661 / 0.002** |
| Maximum cell mean chi | 0.000500699726645103 / 0.001 |
| Minimum sampled J | 0.9999999999998741 |
| Maximum velocity error | 1.7763568394e-15 |
| Worst GPU vs CPU-binary80 chi difference | 4.0604103450e-13 |
| Point energy-rate error / U | 6.3548509595e-17 |
| Point angular-rate error / (2K) | 2.1014263360e-20 |
| Rank0 CUDA curvature time | 48.800773 s |
| Native curvature wall | 53.1955941 s |

The point maximum exceeds its unchanged bound by4.3967 times. Passing the
mean, global energy, angular momentum and linear equation checks does not
override this failure. The original full space had max chi0.00016296225454
under the same curvature definition. The CPU witnesses confirm the new
maximum; this is not evidence of a GPU-only arithmetic error.

The worst points in four symmetry-related quadrants have
r=0.50363027078625. For rank0 the saved witness is cell65488, source triangle
1048862, X=-0.30613511148650929, Y=-0.3999059178949142.

## Failure distribution, without another operator evaluation

`analyze_compact_curvature_failure_20261022.py` reads the saved cell summaries,
canonical geometry and compact mass group metadata. Its output is
`material_compact_initial_failure_analysis_v1/analysis.json` and
`failed_cells.csv`.

- All71,080 cells were checked. There are160 point-maximum failures and zero
  cell-mean failures.
- The5,524 protected cells contain24 failures; the65,556 coarse cells contain136.
- Of268 cells adjacent to a core/coarse interface,48 fail. Another112 failed
  cells are not directly adjacent to an interface. Interface localization of
  the largest peak therefore does not establish the sole cause.
- By cell-centroid radius: r<.45 has no failures, max chi.00172185;
  .49–.5 has16 failures; .51–.55 has16; .55–2 has128.
  The2–8 and8+ regions have no failures.

The reference-frame material space, polynomial resolution and residual
freedom changed together. These results reject this candidate; they do not
prove which change alone causes every error or resolve the historical
t6–10 heating mechanism. No physical-time heating history was generated.

## Terminal state, integrity and adopted action

Native exit is2, expected for the actual curvature rejection. Exec13915 and
the four native ranks3180541–3180544 are terminal. Their CUDA use was directly
observed during execution; all four processes were subsequently confirmed
absent. Final source73, input27 and binary3 checks all report OK. The runner,
environment and source snapshot were frozen before launch and not edited
while running.

`material_compact_initial_completion_audit_v1/contract.json` separately scans
all1,138,432 saved curvature records, verifies their point counts, validity,
positive J and finite metrics, and retains output SHA256. Its PASS means
complete evidence collection of a **rejected** candidate, not a physics pass.

Following `ASTRA_COMPACT_SPACE_AFTER_COST_FAILURE_20261022.md`, stop this
fixed space here. Do not change core radius, order, enrichment count, dt or
tolerances to obtain a small passing run. Do not start cubature preprocessing
or long evolution for it. Preserve its native GPU backend, injection, mass,
solved states and failure witnesses for reuse where valid.

One bundled strategic evaluation was requested from the existing Astra agent
after this complete outcome. No new agent or second physical branch was
started. The next decision must address actual t6–10 evolution and material
shear lifetime rather than replacing the objective with more initial-state
certificates. Weight/mass and W0 reduction, shocks, boundary matching and the
longest normalized GIZMO comparison remain unproved.

## Subsequent bundled Astra decision

The requested evaluation completed as
`ASTRA_EXACT_CELL_CONDENSATION_DECISION_20261022.md`, which was read in full.
It supersedes the compact fresh-branch execution plan. The selected next
bundle preserves the original accepted step5 field and all of its interior
x/v and resolved rotation. It eliminates only cell-interior acceleration
unknowns from each fixed-M linear solve, solves the exact cell-boundary
Schur system with the sparse DD backend, and recovers every interior
acceleration before the original RK4 updates. No compact T or rejected
acceleration is applied to the original trajectory.

The same bundle caches immutable selected-reference enriched values,
gradients and generating j0, with unchanged double bytes and actual-stage
state evaluation. Exact skeleton/support dimensions and factor/CSR/cache
memory counts are required before factor construction or another actual
RHS. The skeleton upper limit is2.5M scalar DOF; physical-device memory
requires20% headroom. Setup wall is bounded8h. Actual full residual must
remain below1e-16, and t1<=24h plus t10<=7days remain the continuation cost
gates. The first actual new RHS is step6 from the preserved accepted time,
with one same-RHS old-solver shadow, not another t0 or initial-curvature run.

This is one explicit revision justified by the new sparse-backend evidence,
not permission for more failed-space radius/order/enrichment sweeps. It does
not solve the separate material-shear lifetime or global-map guard; those
remain prerequisites for long evolution. At this report's completion the
new exact-condensation count/factor/cache implementation has not started.
