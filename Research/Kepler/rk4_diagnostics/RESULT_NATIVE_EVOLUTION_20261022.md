# Actual native RK4 and restart connection — 2026-10-22

## Status

PROGRESS. The retained finite CONTROL now has actual accepted native evolution
and actual restart continuation. This is the 384-source/24-material-cell patch,
not a full-disk evolution. Historical t6–10 artificial heating is still UNSOLVED;
the user's original longest-duration GIZMO goal remains ACTIVE.

Authoritative original repo was rechecked at HEAD
8e25dfe0c2e41a7374c11e2e3102da2a8b028625, origin
`git@github.com:kjhan0606/LagEunha.git`. The supplied workspace is not that repo.
No original repo edits/reset/commit/push or shared-context modification occurred.

## Actual execution evidence

- `shared_evolution_v2_h200_mpi4/contract.json`: PASS.
- Slurm432843, syn104, `--gres=gpu:H200:1 --mem=32G`, MPI4/OMP4.
  sacct COMPLETED0:0; allocation23s including builds, native execution12.5496192s.
- Same saved finite state SHA74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005;
  all817152 original positive material quadrature nodes/weights/rho/U0 restored
  from the four saved Frozen files, without regenerating the rule or resetting U.
- h2e-4 classical RK4: first accepted step, second uninterrupted step, and the
  same second step after reading the actual accepted restart. Three actual steps,
  twelve actual CUDA RHS per rank. Final CONTROL-relative offset0.0004; this is
  NOT the absolute Kepler age and NOT a fresh IC trajectory.
- Actual CUDA volume forces drive the fixed sparse mass solve. CPU point forces
  are independent shadow/audit values. The same canonical boundary G/H, evolved
  thermal boundary pressure, gravity and hybrid reaction ledger are retained.
- Two persistent GPU contexts per rank (initial and actual restart). Five device
  allocations per context; static topology and point uploads once per context.
  No per-RHS device allocations, topology or point uploads.
- Fixed reference raw mass built once from restored quadrature; CSR/preconditioner
  reused during each path. Restart rebuilds the same CSR from its stored raw mass.
- Three full restart sets (imported, accepted, final) include full shared fields,
  canonical triangle/basis metadata, source grouping, original owner IDs,
  reference quadrature/density/alpha, immutable generating coefficients and
  U0/jhat0, fragment descriptors and input SHA provenance, raw mass, stage base,
  time/stage/epoch and work ledgers. Versioned ABI/checksum/durable fsync.
  Stage failure saves its full input state and restores the last accepted base;
  sparse coefficient witnesses are kept for every actually evaluated stage.

## Physics and numerical results

The twelve-stage independent output audit found these maxima:

| Diagnostic | Observed maximum | Existing bound |
|---|---:|---:|
| global semidiscrete energy-rate error / U | 2.484e-15 | 1e-8 |
| cell energy-rate error / U | 3.059e-13 | 1e-8 |
| global angular-rate error / 2K | 1.642e-17 | 1e-11 |
| cell angular-rate error / 2K | 1.395e-16 | 1e-11 |
| CPU/CUDA force relative difference | 5.236e-13 | 1e-11 |
| true sparse mass equation residual | 6.176e-17 | 1e-16 |
| corrected volume identity relative error | 2.060e-16 | 1e-10 |
| stage EOS consistency relative error | 1.061e-16 | 1e-10 |

Finite-step defects are reported separately from the semidiscrete identities:
maximum |Delta E − external work|/U=4.4451e-11,
|Delta L − external torque|/(2K)=2.5163e-15,
|Delta U − RK integral(Udot)|/U=1.4155e-14.
Mass is identical before/after every accepted step; shared trace and bubble
trace jumps are zero. No energy repair, floors, damping or independent spin
energy was applied. Intrinsic L/I/Krot are diagnostics of the actual full field.

Thermal law: immutable generating U0 and actual jhat0=J0+delta0;
U=U0*(jhat/jhat0)^(1-gamma). The first RHS reproduces all saved U0 bitwise.
This is smooth adiabatic evolution with Udot=-P nudot, not shock entropy support.
Neither current J nor recipient delta is substituted for generating thermal data.

The559 fragment source/domain/target references now have a persisted native
full-field view. A pure relabel preserves coefficients, basis and heat bitwise;
target twelve moments do not reconstruct velocity or independent spin.
After evolution these are materially carried fragments of the remap epoch1,
NOT newly evaluated power-cell overlaps. Target packets are not overwritten
or relabeled as new-time aggregates. Source thermodynamic contexts remain fixed.

## Independent saved-output audit

`shared_evolution_audit_v1/contract.json`: PASS. Driver script
`audit_shared_evolution_v1.py` ran on lageunha. It performs no physical queries
and no simulation. Verified27 source,37 input,1 executable and16 output hashes;
all twelve recorded stages and their original gates; all four rank coefficient
witnesses at every stage; actual later coefficient changes; exact original
finite CONTROL coefficient bytes at first RHS; all four restart stage witnesses
versus uninterrupted stages; equality of the two complete second-step ledgers;
all24 saved starting cell U,K,intrinsic L,I,omega,Krot versus prior saved results.
The native loader also requires checksum-valid full actual restart input and
bitwise state/entropy/provenance plus identical time/work ledgers afterward.

## Cost and failed attempt

Rank0 twelve-stage preparation4.0153s, assembly0.4138s, sparse solve1.0056s,
full cell/trace/physics audits2.6901s, timed accepted restart IO0.1305s.
GPU kernels0.56874s including first-use startup; GPU transfer+kernel0.67926s.
These are patch diagnostics with full CPU shadow, not production throughput.
Do not blindly multiply this high quadrature/CPU audit path into a full disk.

Slurm432842 is retained as BUILD_FAILED (5s): duplicate Audit type in the copied
assembly header, before any query/RK step. v2 removed the duplicate and compares
binary80 mass values numerically (padding bytes are not part of a value).
No physics or gate threshold was relaxed. No successful calculation was replayed.
All observed handles are terminal/closed; no monitoring loop or live simulation
is left from this bundle.

## Next required work

Astra was asked once after this completed major bundle to review the actual
source/results and write ASTRA_FULL_DISK_CONNECTION_PLAN_20261022.md. The bounded review is now complete in
ASTRA_FULL_DISK_CONNECTION_PLAN_20261022.md. The next work is actual full-disk
original CUDA/MPI connection and evolution, not more local frozen overlap proofs. Necessary
issues include restricted-patch boundary conditions, fixed material mesh shear,
full-disk initialization and general cell counts, production GPU reductions and
mass solve cost, and eventual conservative repartition/context re-expression.
Then execute actual t6/8/10 density/heating/orbit/E/L checkpoints and continue to
KEPLER_GIZMO_COMPLETION_CRITERIA_20261022.md. Shock entropy support and a completed
new-time Laguerre remesher are not proved by this smooth material step.

## Adopted next bundle, after Astra review

The plan selects an opt-in mode in the original native main/RK execution
lifetime, using the existing71,080-point fine_integrated fresh continuous IC.
The proven patch is not naively expanded: runtime sizes, streamed GPU
field/EOS/reductions, owned/halo sparse mass matvec and bounded CPU shadows
are required. Preserve an isolated native source/build snapshot and all
legacy binaries/results; do not overwrite the authoritative user's worktree.

The first full material research path has explicit constant prescribed external
pressure Pext=P0 on a moving material boundary. It is NOT an identical periodic
baseline or final matched-paper run: the exact rotating IC is not periodic
across square boundaries, and silently identifying its H1 DOFs is invalid.
This boundary difference remains an open final comparison requirement.
No fictitious fresh target overlaps or final W0/Voronoi success is claimed.

Next execution is one preserved continuous trajectory with physics diagnostics
at.1/.2/1 and actual continuation toward6/8/10, first bounded by t1 or24h.
Full-state initial curvature/positivity or a topology/solver failure is preserved
and stops that run without automatic threshold/dt/order scans. The t1/wall
bound is a checkpointed batch, not a new final objective. No more24-cell
3-step or completed transfer-query replay is required.
