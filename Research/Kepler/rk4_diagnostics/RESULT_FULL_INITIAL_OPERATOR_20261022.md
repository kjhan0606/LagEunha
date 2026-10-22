# Full native initial state, curvature and cell/rotation ledger — 2026-10-22

The long Kepler objective is ACTIVE. This bundle does not demonstrate suppression
of the historical t6–10 artificial heating, or any accepted full-domain time step.

## Adopted actual evidence

- Original user repository: `/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src`, HEAD `8e25dfe0c2e41a7374c11e2e3102da2a8b028625`, remote `git@github.com:kjhan0606/LagEunha.git`; original source edits/objects and old results preserved.
- Actual initial field: 71,080 source cells, 1,138,432 reference triangles,
  10,816,441 scalar DOFs, owner/ghost MPI4 and OpenMP4 per rank.
- Slurm432854 COMPLETED 0:0: ONE H200 shared by four ranks, 64GiB request,
  actual fixed mass and gravity over 2,040,256,512 positive quadrature points.
  Source/binary/input/object manifests and native initial state are preserved in
  `material_full_initial_v2_h200_mpi4`. Failed432853 was a compile-only missing
  namespace alias, no physical initialization. Pending request was reduced from
  H2004 to H2001 on the same pending job because QOSGrpGRES blocked four GPUs.
- Slurm432855 COMPLETED 0:0: actual initial compression curvature and point work
  over the same full field using saved acceleration, with NO mass solve replay.
- `material_full_initial_audit_v1/contract.json`: saved-output independent all-cell
  mass/owner/topology/raw-block/fixed-matrix audit PASS; no new field queries.
- `material_full_initial_ledger_v1_mpi4/run/full_initial_local_ledger_result.json`:
  actual all-cell hybrid reaction/energy/angular/rotation ledger PASS on lageunha
  in 7.73137443 seconds. Consumes saved M, actual acceleration, raw initial forces
  and curvature metrics, no physical point queries or solve/factorization.

## Numerical results

| Quantity | Actual result | Existing bound |
|---|---:|---:|
| Initial mass equation true residual x/y | 9.79e-17 / 9.79e-17 | <1e-16 |
| Initial global E rate error/U | 1.44e-15 | <1e-8 |
| Initial global angular rate error/(2K) | 6.66e-19 | <1e-11 |
| CPU/CUDA mass shadow relative error | 1.98e-12 | <1e-11 |
| CPU/CUDA initial force shadow relative error | 1.97e-13 | <1e-11 |
| Full maximum initial compression chi | 1.6296225454e-4 | <=2e-3 |
| Full maximum cell mean chi | 4.6798757409e-6 | <=1e-3 |
| Binary80 CPU worst-point curvature difference | 4.30e-15 | <1e-8 |
| Maximum cell E rate closure/U | 2.3034266395e-12 | <1e-8 |
| Maximum cell angular rate closure/(2K) | 1.7623260870e-17 | <1e-11 |
| Zero-trace bubble reaction/load | 2.9914978903e-16 | <1e-11 |
| Owner-summed shared reaction/load | 3.9997963259e-16 | <1e-11 |
| Maximum cell mass error vs retained integrated IC | 5.4860127434e-11 | <1e-8 |

CPU shadow covers the fixed 20 full triangles (five per rank), not every GPU
point. Worst-curvature CPU witness is one actual saved GPU maximum per rank.
The independent initial audit used numeric float equality when labeling ghost
coefficients bitwise; the later native ledger explicitly checks double bits,
including signed zero, and validates strictly increasing global IDs without
unsigned differencing. This supplements that audit without replaying it.

Physical initial U=0.000863999588012891662344; fixed-metric K=8.78884830216854721387.
Resolved intrinsic rotational kinetic energy=7.21604175570218798069e-5. It is a
diagnostic decomposition of actual cell velocities; no extra spin reservoir,
relaxation, damping, steering or heat correction was added. Nonzero individual
cell shared-trace reactions represent exchange and are not artificial heating.

## Measured cost and next actual work

Initialization native108.475352s: basis6.3184s, mass integration22.4229s,
CUDA mass kernel rank0 16.9019s, conditioning0.9900s, and fixed local factorization
PLUS two CG solves72.2680s. The last number is not a measured per-RHS solve cost.
Curvature native19.118655s, CUDA rank0 13.395415s. All four ranks shared one GPU;
per-rank CUDA timings cannot be summed or described as four-card throughput.

`ASTRA_FULL_EVOLUTION_PERFORMANCE_PLAN_20261022.md` is the adopted next bundle;
no further primitive review is requested. New native code now connects immutable
fresh constant-P generating entropy, three-pass changed-state CUDA geometry/EOS/
force, canonical mean pressure and prescribed constant external P0, original
binary80 mass solve, classical RK4, per-cell energy/angular/rotation ledgers,
rollback and full immutable-reference restart. CPU mass fallback is explicit;
GPU DD correction remains a required acceleration task. Local factorization is
parallelized and its timing separated; element preconditioner scratch reused.

Slurm432869 is the submitted first full RK4 and actual restart-comparison bundle
in `material_full_evolution_v1_h200_mpi4` (H2001, MPI4/OMP4,64G,2h). Its submission
is not evidence of a successful evolved stage or accepted step. Check that job's
actual status and saved results before advancing. It consumes the completed full
initial state; it does not rebuild tessellation/basis or repeat fixed mass points.

## Remaining scope

Actual changed-stage pressure and full accepted evolution/restart must pass.
Then one continuous trajectory must reach t1,6,8,10 and properly normalized
longest GIZMO comparison. Do not repeat accepted intervals. Performance must be
measured before a long submission. Moving constant-Pext boundary differs from
the original periodic baseline. Material-map folding/global injectivity,
recipient thermal context on remeshing, dynamic W0 Voronoi reduction and shock
entropy are unresolved; none is proved by these initial-state checks.
