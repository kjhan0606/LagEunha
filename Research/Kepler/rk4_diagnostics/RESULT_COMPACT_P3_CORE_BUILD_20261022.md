# Compact P3 core implementation and native object build

2026-10-22. Adopted next major plan:
`ASTRA_CONSERVATIVE_COMPACT_REZONING_PLAN_20261022.md`. Driver implemented the
initial reusable components. This is code/build progress; no compact physical
transaction, state transfer, checkpoint import execution, or new step has run.
Accepted trajectory remains step5, t=0.00011263377801805954.

## Implemented components

- `material_compact_p3.h`: a straight-fan shared P3 scalar registry with canonical
  vertex/edge/interior IDs, orientation-aware edge DOFs, positive reference
  geometry, manifold incidence and disk Euler prerequisite checks. The registry
  takes actual current canonical power-cell polygons and their area centroids;
  it does not generate analytic-orbit sites, merge coordinates by tolerance,
  or inherit the old h/2 source subdivision.
- Pure P3 shape/value/gradient evaluation uses the original native coefficient
  and physical-point ABI. There are no calls to IC residual evaluation in this
  compact evaluation path. Current positions, velocities, J and actual velocity
  gradients are resolved. No spin reservoir is introduced.
- Independent reference mass-density and physical internal-energy density use
  ten Bernstein P3 coefficients per triangle. Negative controls are rejected,
  never clamped. Nonnegative controls alone do not prove strictly positive
  continuous physical fields; that still belongs to the transaction gates.
- A positive 6x5 Duffy/Gauss mass rule integrates total-degree9 polynomials in
  exact arithmetic. Thus rho(P3)*N(P3)*N(P3) can build a consistent symmetric
  mass matrix without the old nonpolynomial IC basis. Floating arithmetic
  errors and nonpolynomial evolved EOS/gravity still require actual checks.
- `material_compact_checkpoint.h`: exact v3 accepted-payload reader, pinning
  owner/count/time/step, immutable physical and positive-rule binding, checksum,
  original binary80 ledger ABI and coefficient ABI. Physical and numerical
  budgets are retained independently, with no k1 or initial-force replay.
  Caller must supply the independently verified binding/generation; the
  actual step5 import has NOT yet been executed.
- `material_compact_p3_cuda.h` / `material_compact_p3.cu`: asynchronous CUDA
  evaluation over caller-owned epoch-specific device buffers, returning the
  physical map/velocity/gradient and reference density/energy at supplied
  samples. This is a field sampler, not a completed EOS/force operator.

## Actual build evidence

All builds ran on authorized lageunha, not the login node. Source snapshots
and failure logs are preserved. `material_compact_rezone_build_v1` stopped at
host compilation because the old header undefines its private ST_HD macro and
the new header needed `<algorithm>`. Driver fixed the new header with its own
host/device macro and explicit include; old files were not changed.

`material_compact_rezone_build_v2` finished exit0:

- GNU13 C++17 host header syntax check PASS, including v3 ABI static asserts.
- CUDA13 sm89 Ada object build PASS (`compact_ada.o`).
- CUDA13 sm90 H200 object build PASS (`compact_h200.o`).
- Frozen/root new sources match. Source/object integrity checks PASS.
- Original frozen native Exam header dependency hashes recorded.
- Sessions39728 CLOSED1,24139 CLOSED0,57431 CLOSED0.

`build_contract.json` explicitly records full native application link=false,
runtime=false, actual step5 import=false, projection=false, geometry proof=false,
compact commit=false, new step=false and long Kepler validation=false. No GPU
was reserved or physical calculation launched for these compilation checks.

## Remaining work in this one major transaction

1. Bind full original step5 basis/ownership, fixed M, immutable generating
   entropy and current corrected thermal context to the exact imported state.
2. Establish actual full-source injectivity and coverage as a transaction
   prerequisite. The current exterior is a material curved boundary. A straight
   power-cell target registry must establish the same physical domain; silently
   clipping to the initial square is not a proven substitute. Existing finite
   CONTROL packets and sampled positive J cannot establish this prerequisite.
3. Generate actual transported sites, canonical target registry and approved
   finite weight-linked mass law; assemble actual source/target overlaps and
   conservative m/q/Q/p/L/K/U/Phi packets. No obsolete packet is reused as the
   current field.
4. Solve the fixed-space positive density/energy fit and coupled shared-velocity
   p/L/K projection with finite sparse-solver budgets. Preserve physical K and U
   separately; no heat/spin reservoir, arbitrary high-frequency mode or repair.
5. Verify field/entropy approximation and remaining T10 error budgets, then
   connect the same variational pressure/work law and native mass/RK4/restart.
6. Only a passing full transaction permits the next genuine native step and
   amortized t1<=24h cost decision, including necessary remap/guard cost.

There is no separate new small guard experiment or long old-space continuation.
The complete planned batch is bounded to24h and preserves old accepted outputs
on rejection. Original CUDA/MPI/OpenMP infrastructure is reused; the original
user repository and all completed frozen binaries/results remain unchanged.
No global shared context was edited. Full t>6 and longest normalized GIZMO
stability are still the goal, not established by these builds.
