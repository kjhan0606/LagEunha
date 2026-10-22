# Geometry-completed space and original resident GPU/MPI mass connection

2026-10-22. Driver implementation under the same accepted-state compact
transaction, following `ASTRA_COMPACT_CURVED_DOMAIN_ADDENDUM_20261022.md`.
No additional Astra call, physical simulation or source-guard experiment was
made. The previous explanatory turn advanced no calculation; this turn adds
native/CUDA code and executes new finite algebra/API checks on lageunha.

## Geometry completion implemented

`material_compact_geometry.h` defines a fixed accepted exterior trace and
the collar chart G=X+s^4*(Gamma-chord). It evaluates the original frozen trace
and its derivatives, not an analytic orbit or a replacement initial velocity.
The P3 interpolant of G is the polynomial geometry; h=G-I3(G) supplies two
scalar functions per connected component, normalized without discarding
small singular values. Interior patches retain pure P3.

The same local basis and shared coefficient buffer are used for position,
velocity, deformation/velocity gradients, 12x12 mass blocks, symmetric-stress
work and external pressure traction/work/torque. The two geometry coefficient
vectors are dynamic unknowns; no separate spin or heat reservoir is added.
`material_compact_geometry_cuda.h/.cu` expose asynchronous field, stress and
boundary APIs for caller-owned device buffers.

The density tag is explicitly parameter-reference-volume density. Initial
det(DG) belongs in the source-to-reference projection. Mass integration does
not multiply it a second time. Curved mass assembly requires an explicitly
supplied positive functional; straight degree9 exactness is not a certificate
for curved nonpolynomial collars.

These are components. Rounded endpoint equality is not a proof of canonical
event equality/coverage. Actual chart injectivity, exact component gluing,
geometry-mode function rank, quadrature/error bounds, MPI geometry ownership
and checkpoint serialization remain to be connected. The current adapter
retains two independent candidate modes; it has no tolerance-based rank
reduction. A genuine lower-rank case needs exact algebraic handling or rejects
the candidate. The stress helper does not supply the complete entropy/EOS
or nonlocal geometric-pressure volume derivatives.

## Mass solver compatibility problem resolved in the adapter

The original native and resident GPU preconditioners require positive local
Cholesky factors. Added component modes can have zero/dependent support on
an individual patch while the assembled physical mass is positive. Factoring
every physical local block would therefore reject valid assembled systems.

`material_compact_mass.h` reuses the original Ownership and Mass interfaces.
It constructs positive local PRECONDITIONER blocks from each physical P3
principal block and the assembled 2x2 component-mode metric divided by its
support count. Bulk/mode coupling is omitted only in that preconditioner.
Before any matrix apply or device upload, the original Mass.block is replaced
by the untouched authoritative physical blocks, including every coupling.
No mass, kinetic energy or physical equation is regularized.

This allows the unchanged original DD GPU/resident inner solver and binary80
outer true-residual check to operate on the completed space. It avoids a new
global dense solve. Component IDs are compared exactly across MPI ranks.
The finite-functional positive 2x2 check is necessary, not a substitute for
continuous geometry-function rank or positive-quadrature proof. Prepared
objects are heap-owned/nonmovable because original Mass stores references;
Ownership and Prepared must be destroyed before MPI_Finalize.

## Actual new runtime evidence

These are bounded synthetic algebra/API fixtures, not physical Kepler runs,
accepted-source guards or conservative-remap feasibility results. No saved
checkpoint, completed quadrature producer, old RHS or accepted interval was
queried again.

`material_compact_geometry_algebra_v1_ada`:

- 34 CPU/Ada points, including apex and curved exterior endpoints.
- Two independent polynomial completion functions, arbitrary dynamic mode
  velocities, positive reference density and a finite mass functional.
- Rigid-rotation gradient error 7.13e-16.
- Mass/point K and L differences 2.13e-19 and 3.18e-18.
- Symmetric-stress energy/torque residuals 1.11e-17 / 1.57e-17.
- Exterior work/torque residuals 3.45e-24 / 3.90e-24.
- CPU/GPU maximum normalized difference 2.89e-15.
- PASS, native0/runner0, wall0.288s. Session59744 closed0.

`material_compact_mass_resident_v1_ada_mpi2`:

- Two MPI ranks, 22 global owned DOFs, two deliberately rank-deficient local
  physical blocks. The original physical-block Cholesky rejects both.
- Adapter preserves every physical block entry; original GPU matvec differs
  from the native authoritative product by exactly0 on this fixture.
- Original resident solver: three inner iterations, one outer correction,
  true relative residual1.8565069529462421e-21, solution difference3.39e-21.
- Exact binary80-to-DD conversion, 1152 actual sparse ghost-pair bytes.
- PASS, observed SSH execution85176 closed0, wall2.304s. See its explicit
  runner_observation.json for the absent shell-tail receipt (consistent with
  mpirun consuming SSH-shell stdin). No rerun occurred.
  Future SSH-shell mpirun commands must redirect stdin from /dev/null so MPI
  does not consume the remaining shell script.

The first fixture uses a synthetic polynomial curve, not an actual enriched
step5 exterior. Neither fixture verifies source embedding, target coverage,
long-time accuracy, real transfer feasibility or sustained production cost.

## Build and integrity receipts

`material_compact_rezone_build_v5`: original source/context plus new geometry
compiled as native MPI/OpenMP; new geometry CUDA objects compiled for AdaSM89
and H200SM90. Relocatable component links include the preserved original GPU
objects. Its native TU is subsequently superseded by v6 for the mass adapter;
the frozen v5 geometry objects are reused unchanged.

`material_compact_rezone_build_v6`: new native source/mass adapter compiled,
then linked with those unchanged CUDA objects and the original Ada/H200 GPU
objects. Original mf_initial_/mf_evolved_ APIs resolve in both component links.
Frozen/root sources match, source/object/retained-object/dependency integrity
checks pass, 32 compiler project dependencies recorded. Session1620 closed0;
final receipt69926 closed0. This is not a full transaction executable.

Both fixture source/binary/input/dependency manifests were checked. All new
tool sessions are closed; no new Slurm job, simulation or monitor was launched.
Original user repository HEAD8e25dfe, completed source/rules/outputs and shared
global context remain unchanged.

## Actual physical state and next work

Accepted state remains step5 at t=0.00011263377801805954, prefix
`material_full_positive_transition_v2_ada_mpi4/run/continuous_step5/accepted_r`.
No compact representation was committed, no new accepted time interval ran.

Next integrate the actual source-certified curved power-cell partition and
finite W-linked mass law, positive overlap quadrature, rho/e and shared p/L/K
projection, field/entropy/error ledgers, geometric-pressure area derivatives,
MPI geometry/restart and complete native controller. Only a passing transfer
permits the next physical RK4 step. Include collar/remap/guard costs in the
existing amortized t1<=24h gate. Do not use these fixture PASS results to
promote a physical state or restart old10.8M-space long production.

t6/8/10 and the longest properly normalized GIZMO test remain unverified;
artificial shear heating is unresolved and the full objective stays active.
