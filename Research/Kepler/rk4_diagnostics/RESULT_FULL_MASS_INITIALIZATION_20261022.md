# Full distributed native mass initialization — 2026-10-22

## Actual status

PROGRESS. The full71080-cell native material reference now has distributed
canonical DOFs, real shared-owner communication, actual streamed GPU fixed-mass
integration, normalized local SPD blocks and persisted state/mass operators.
This is not yet a hydrodynamic step. No full-disk RK4, t6/8/10 stability or
GIZMO longevity is proved; the original goal remains ACTIVE/UNACHIEVED.

Authoritative execution: shared_full_initialization_v3_h200_mpi4/contract.json.
Slurm432852 COMPLETED0:0 on H200, allocated42s including hashes/IO; native
23.9359492s. MPI4/OMP8, ONE actual H200 shared across4ranks (all selected GPU0,
identical UUID),128GiB requested host memory. Do not report4GPUs. Selected
NVIDIA H200 NVL UUID hex c5feb54b72313bc24cffb61c4cf3bac3.

## New implementation

All production development remains in isolated
native_full_material_20261022/Exam/material/; old original repo/binaries/results
and user-wide SHARED_CONTEXT remain untouched. Original identity reconfirmed
by SSH syntax: HEAD8e25dfe0c2e41a7374c11e2e3102da2a8b028625,
origin git@github.com:kjhan0606/LagEunha.git. Direct tools run on LagEunha; the
worktree Git metadata is syntax-local. Slurm submission must go through syntax.

- material_parallel_domain.h: full canonical basis/IDs built once, cells sorted
  spatially by generator y/x, whole material cells assigned across4ranks, unique
  DOF owner selected from supporting triangle owners. Scatter compact local
  triangles/coefficients/global IDs; local IDs replace GPU array indices. Real
  physical external edges retain global IDs and are distinguished from MPI cuts.
- Halo: one-time ghost-ID requests, owner pulls and reverse additive reductions.
  Only shared values travel each call. Actual owned-to-ghost ID transfer PASS.
  No replicated full global CSR or per-stage element-array Allreduce.
- Canonical edge Gauss64 nodes/weights cached once with unchanged arithmetic;
  independent edge bases and triangle caches use OpenMP. No degree/rank cutoffs
  or quadrature reduction. Global topology initialization13.168s.
- material_stream_mass.{h,cu}: native fixed quadrature and original st_shape,
  reference density; positive point values generated procedurally within each
  triangle block. Shared128-point tile, paired-binary64 product/sum reductions,
  convert to binary80 local fixed blocks. Static triangle/rule data stay on GPU;
  at most8192 triangle results are transferred per batch. No STPoint/GPOutput
  array for all2.04bpoints, no per-point dumps.
- material_parallel_mass.h: normalized local blocks and packed binary80
  Cholesky, deterministic row-incidence OpenMP matvec/gather, additive local
  block preconditioner and owned-only scalar reductions. Distributed PCG code
  retains true residual<1e-16 and explicit bounded failure, but PCG has NOT yet
  been exercised with a physical full-state RHS. Do not claim that solve gate.
- Native reference-state and fixed-block loaders implemented for continuation
  without recomputing the accepted physical mass quadrature. These are reference
  initialization archives, not yet complete evolved-heat/RK checkpoints.

## Actual full data and checks

Global10816441scalarDOFs,1138432triangles,71080materialcells.
Actual positive GPU mass queries2040256512, exactly the previously constructed
fixed-rule count. Geometry/count-only outputs from the prior turn were not
repeated field queries; this was their first actual full physical mass integral.

| rank | owned DOFs | owned+ghost DOFs | triangles | actual mass points |
|---:|---:|---:|---:|---:|
|0|2703607|2703607|284320|509501440|
|1|2705918|2708489|284896|510626816|
|2|2705880|2708489|284896|510626816|
|3|2701036|2703607|284320|509501440|

Physical integrated mass75.9742060105944259363. Relative difference to the
preserved continuous-IC total5.35054e-17, without normalization or heat changes.
Constant-field mass metric75.9742060105943994852; relative difference from
quadrature mass3.48137e-16. Every local normalized Cholesky SPD. The constant
field matvec exercised reverse shared accumulation on this actual full state.

First full-state CPU shadows were selected at nearest radius0/2/8/17 on each
spatial rank; duplicates yield14distinct actual source triangles. Raw and
normalized witnesses were saved, not only aggregate PASS bits.
The native reported ratio(max absolute difference)/max scale is2.59183e-14 for
raw blocks and9.81735e-15 for normalized blocks. These are NOT max per-source
relative errors: shared_full_mass_audit_v1 independently recomputed each source
comparison from saved matrices, finding actual maximum1.7953736231564463e-13
against unchanged1e-11 native CPU/CUDA shadow bound. Audit PASS. No point,
density, geometry or simulation replay occurred in this saved-output audit.

Per-rank CUDA kernel totals2.2324–2.4241s and mass initialization8.85–9.09s
include concurrent contexts on the same card; do not sum/interpret them as four
independent GPUs. Full execution runtime includes topology, GPU/host IO,
normalization, preconditioner and native state persistence, not evolution cost.

## Persisted artifacts: NEVER REINTEGRATE SUCCESSFUL FULL MASS

shared_full_initialization_v3_h200_mpi4/ contains, for each rank:

- topology_rankN.bin: initial raw-basis reference, canonical global/local IDs,
  ownership, source grouping, physical external edge IDs and raw IC coefficients.
- raw_mass_stream_rankN.bin: completed ordered binary80 raw blocks + triangle
  masses, retained as a durable streamed prefix; all complete in this run.
- normalized_state_rankN.bin: normalized basis metadata and IC coefficients;
  physical initial field preserved by the basis transformation.
- mass_operator_rankN.bin and triangle_mass_rankN.bin: transformed fixed local
  operator and actual integrated material mass. Reuse with reference-state
  loader and rebuild factor/halo metadata when resuming; no quadrature replay.
-14mass_shadow_rankN_sourceS.bin witnesses, device_rankN.json, rankN.json,
  source/input/executable/output SHA manifests and final integrity results.

No dynamic physical heat state or accepted full RK step has been saved yet.
Generating heat context must be bound once in the next real full-state RHS;
never treat a topology or mass checkpoint as a complete evolved state.

## Attempts and resource handling

- v1: direct sbatch on LagEunha failed DNS configuration; NO Slurm job/build/
  physical query. Use ssh syntax for scheduling, not retries on local sbatch.
- v2 Slurm432848: H2004 pending QOSGrpGRES. Reduced pending request to H2001,
  preserving MPI4. It started, but CUDA compile failed because stdint was absent.
  Terminal FAILED1:0,2s, no physical query. Explicit header fixed in v3.
- v3: compiled native/CUDA on authorized LagEunha first, sm90 with fmad off;
  all compilation PASS. SHA-bound prebuilt executable submitted on syntax with
  H2001/MPI4/OMP8/128G. Actual full queries run only once. CPU/CUDA precision
  bounds, polynomial order, mass and heat were not relaxed.

All specific process/session/job handles from this bundle are terminal/closed.
No monitoring loop or live simulation is left. Metadata inspection respected
specific job handles; no repeated cluster/file-tree scans or old result deletion.

## Next actual work (already authorized; no new Astra call needed)

1. Load normalized_state/mass_operator from this completed full initialization.
   Do not rebuild full IC/topology, old23.339m transfer atlas, successful24-cell
   RK4/restart, or these2.040b mass points.
2. Build source-boundary canonical functional/cache from retained triangles,
   separating physical external edges from MPI cuts. Preserve the same finite
   quadrature G/H derivative/work contract. Apply explicit prescribed Pext=P0
   on the moving external boundary; label the change from periodic baseline.
3. Actual streamed CUDA same-stage J/QJdot, evolved heat/mean pressure, volume
   stress/gravity/thermal/rotation diagnostics and shared force assembly;
   generating heat is a function tied to immutable fresh full IC and delta0,
   not a new per-stage P0*J law or a2b-entry entropy buffer. Fixed CPU shadow
   on first full RHS only; semidiscrete and finite RK defects stay separate.
4. Exercise distributed true-residual solve on the first physical RHS and
   connect original Exam/exam.c opt-in native MPI lifetime. Then actual full
   classical RK4, full thermal/rollback/restart state, positive J/nu and map
   topology controls; keep first accepted step as start of one trajectory.
5. Measure actual seconds/step and run bounded t1/24h batch, then same history
   through t6/8/10 and KEPLER_GIZMO_COMPLETION_CRITERIA. Density/heating/orbits,
   annular/global E/L and boundary work are mandatory. Entropy fixed by the
   smooth model is not by itself evidence that false heating is solved.

Dynamic Voronoi W0 reduction, remesh/recipient context representation, shock
entropy support and longest-duration stability remain unproved. No claim of
periodic or paper-equivalent full evolution is supported by this initialization.
