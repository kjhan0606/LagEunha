# Full native material reference connection — 2026-10-22

## Current status

PROGRESS toward the unchanged original long Kepler goal. The original native
MPI/initial-condition/CUDA geometry lifetime now has an isolated opt-in full
material-reference preparation branch. It executed on the preserved71,080-cell
fine_integrated IC, constructed a conforming canonical full reference and compact
streaming quadrature plans, and passed independent saved-output comparison.
A distributed owned/ghost fixed mass operator was also implemented and builds.

No full-disk physical field/mass initialization, actual full RHS, accepted full
step, t6/8/10 or longest-paper evolution is proved yet. The prior actual accepted
RK4/restart result remains the24-cell finite CONTROL, not this full disk.
The user's t>6 artificial-heating issue remains UNSOLVED and goal ACTIVE.

SHARED_CONTEXT was read fully; supplied cwd and original repo were distinguished.
Original repo HEAD8e25dfe0c2e41a7374c11e2e3102da2a8b028625,
origin git@github.com:kjhan0606/LagEunha.git, was confirmed before edits/actions.
Existing dirty Exam/exam.c and Exam/exam_gpu.cu were preserved. All new native
edits/builds are isolated under this diagnostics root. No original reset/edit,
commit/push or user-wide context modification occurred.

## Native entry and actual full-domain execution

- Plan: ASTRA_FULL_DISK_CONNECTION_PLAN_20261022.md (already reviewed; no new
  repeated Astra request this turn).
- Sources: material_full_native.inc, material_full_reference.h/.cpp.
- Isolated Exam/exam.c adds GFS_MATERIAL_PREPARE_REFERENCE before legacy/DG
  branches. It validates fresh time0/W0/no steering/Re-infinity/no heat switch,
  then uses original exam2dUpdateVol and actual raw complete CUDA polygons.
  The explicit preparation branch exits collectively and cannot fall through
  into legacy particle RK scratch. It is not labeled an evolution mode.
- Slurm432845 on H2001/MPI4/OMP4/32GiB, COMPLETED0:0; allocation10s including
  reused native objects and new host build; actual full geometry/reference
  preparation1.82033879s. All four ranks record the same actual H200 UUID;
  this is one GPU card, not four.
- Authoritative result material_full_reference_v2_h200_mpi4/run/material_full_reference.json.
- Original eunha2.o/MPI/main lifecycle and native archives linked by the proven
  original recipe, with isolated Exam archive and fresh original CUDA object.
  Copied native sources, exact header closure, build commands, executable,
  dependency files, actual mapped dynamic libraries and input SHA are retained.
- Slurm432844 compile-only failure is retained: new host builder used the wrong
  runtime UUID API name. Corrected to cudaGetDeviceProperties(...).uuid; original
  C/CUDA objects had built successfully and their exact hashes/source closure
  were reused in432845. No native query or preparation ran in432844.

## Full reference, not a coordinate proximity merge

All raw nonzero face pairs are joined using reciprocal global cell/support-plane
IDs. Exact binary128 dyadic support determinants collapse only proven zero-length
edges and verify shared corner identity. The retained IC has dyadic512 lattice
coordinates bounded by64, so required determinant products fit113 significand
bits. No distance tolerance chooses topology. Periodic image support planes are
retained for geometry but opposite square boundaries are NOT identified as a
shared rotating velocity DOF.

-71,080 real cells;350,144 raw face ends;65,536 exact zero edges collapsed.
-141,792 paired physical faces and1,024 explicit square external faces.
-570,241 refined reference vertices;1,138,432 triangles;1,708,672 edges.
-One centroid fan and one h/2 conforming subdivision, the same reference rule
 used by the retained patch; full reference manifold and Euler characteristic1.
-3,784 rank1 and1,704,888 rank2 canonical edges:10,816,441 prospective scalar DOFs.
-Total area576; minimum original cell area0.0001373291015625.
-Native/canonical vertex maximum shift1.536e-16; native raw/canonical cell area
 relative difference maximum1.683e-15. No original generator is moved.

The first material research path has the explicit moving boundary/Pext=P0
contract selected by Astra. This is a boundary change from the periodic baseline;
it is NOT final matched-boundary/paper stability evidence. No analytic position
or velocity boundary forcing is introduced. Actual boundary work belongs in the
future full evolved-state ledger.

## Streaming and mass implementation

The original64x28 fixed reference rule is preserved in compact per-triangle
GfsEnrichQuadrature plans:6,120,597,504 slots, estimated2,040,256,512 positive
slots. No physical point/field quadrature was evaluated in this preparation.
`material_full_resource_v1/full_reference_memory.json` uses actual compiled ABI
and actual counts. An unstreamed STPoint+GPInput+one GPOutput would require
about1.11TiB; even a full point array alone is about121.6GiB. These are capacity
counts, not measured evolution runtime. No such arrays were allocated.
A131072-point batch is75MiB for those three arrays, before reductions/caches.
Reference mesh and compact plans are durable; complete point/GPOutput dumps
are absent. Global fixed local blocks and factors are each about5.5GiB and
can be partitioned over ranks; no full replicated global CSR is planned.

`material_distributed_mass.h` implements native binary80 owned/ghost layouts,
64-bit peer lengths with bounded MPI chunks, owner-to-ghost pulls, partial
contribution sums to owners, owner-only global inner products, matrix-free
local element matvec, local Cholesky/additive preconditioning, true-residual
PCG at the unchanged1e-16 bound and exact zero-RHS handling. OpenMP rows have
no scatter races. Communication carries shared DOFs, not the full global field.
`material_full_runtime_abi.cpp` builds this interface on lageunha and records
actual sizes. This is code/build evidence ONLY: no full mass matrix or solve
has been constructed or validated, and it is not yet linked into a full RK RHS.

## Independent saved-output audit

`material_full_reference_audit_v1/contract.json`: PASS, on authorized lageunha.
`audit_material_full_reference_v1.py` evaluates no physical fields and no RHS.
It verifies source/input/executable/reused-object hashes, all four actual native
shards and UUIDs, unique global cell ownership, all centers bitwise against
geometry.npz, all positive triangle areas, integer-ID refined manifold/Euler
closure and square external boundaries, compact rule lengths/counts, and all
71,080 cell areas against saved geometry (relative maximum9.8056e-15).
`verified_output.sha256` binds actual outputs; `cell_native_owner.npy` preserves
native MPI cell ownership for the next distributed state construction.
All observed Slurm and SSH/build/audit handles from this turn are terminal/closed.

## Next concrete implementation, without repeating successful work

Consume material_canonical_reference.bin, compact rule plans, native ownership
map and original IC SHA; do not rerun this completed topology constructor or
previous24-cell trajectories/23,339,008 transfer queries. Construct owned local
STTriangle/edge bases and global-to-local DOF IDs, immutable fresh U0/jhat0,
streamed fixed mass blocks and initialized full coefficients. Connect the new
owned/ghost mass operator and actual two-pass streamed GPU field/EOS/stress/
reduction to the original native optional material evolution mode. Then the
first actual full RHS/step must validate original physical/initial-curvature
and true-solve gates; retain all failure inputs. No automatic thresholds/order/
dt scans and no heat repairs.

Only after full state/RHS exists can a real continuous trajectory run, first
bounded by t1 or24h with checkpoints. It must then continue to actual6/8/10 and
KEPLER_GIZMO_COMPLETION_CRITERIA_20261022.md. Preserve the boundary difference,
material-map topology lifetime, recipient-context/remeshing, dynamic W0/Voronoi
reduction and shock entropy limitations as unfinished; no local/preparation
success replaces the user's requested end state.
