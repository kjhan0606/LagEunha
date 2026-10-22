# Existing material timestep provenance — native build, 2026-10-22

Goal remains ACTIVE: historical t6–10 artificial heating and the longest
normalized GIZMO comparison are NOT validated. Current accepted trajectory
remains step4,t=8.2355880701376864e-5. The previous goal turn made progress
through the independent8192-record prefix audit, then verified the original
producer and single dependency live. Neither completed work is repeated here.

## Changes made while the original preprocessing continues

ASTRA_MATERIAL_SHEAR_LIFETIME_AND_GUARD_20261022.md section3 recommends recording
the actual minimum source/branch in the existing reduction, rather than inferring
it from the origin analytic match. Implemented optional host-only provenance:

- `material_full_dt_witness.h`, connected in both streamed and cached evolved
  host paths. Enabled only by `GFS_MATERIAL_SAVE_DT_WITNESS=1`; default is off.
- The original eta/gradient, eta/Omega, and altitude/sound candidates and the
  existing e.dt_bound update expressions remain unchanged. Only their argmin is
  recorded. Canonical source ID and unique owner rank break exact global ties.
- Output includes all three candidate values, limiting branch/source/cell/rank,
  actual block extrema, reference centroid radius, coefficient epoch, stage,
  rule generation and whether the reduced rule is active.
- It records the triangle giving the existing block-extrema bound, NOT an
  unsaved extremum-point location. Min length and max sound need not occur at
  the same point. It does not claim an enriched spectral CFL or authorize a
  larger dt from a selected subset. The paired oracle's full-node constraints
  remain unchanged. No physical coefficient, heat, entropy, M or ledger is reset.

The root changes are ONLY future host code. The active producer, its frozen
68source closure, both already built pending transition binaries, dispatcher,
immutable rule generation and prepared highk1 are untouched. This adds no
force evaluation, guard experiment, new job/monitor, or Astra review.

## Actual compile/link and narrow verification

Frozen build `material_dt_witness_build_v1`:

1. Actual full native C++ driver including both force paths compiled on lageunha
   with original MPI/toolchain and precise FP options; source hashes pass.
2. A single synthetic4-rank metadata/serialization fixture checks deterministic
   global ties, owner broadcast, branch output, exact unchanged original dt and
   physical ledger, and default-off no-op. It passes. The result explicitly
   labels synthetic metadata and physical_runtime_validation_complete=false.
3. Original H200 native archive was copied into this isolated directory and its
   C++ driver member replaced. Its CUDA ABI/source are unchanged. The original
   native link succeeded, dynamic libraries resolve, output/archive hashes are
   saved in native_binary.sha256. The inherited ifx option warning is recorded;
   actual link exit is0. The binary was NOT physically executed.

Actual synthetic result is `metadata_fixture/verification.json`; its dt=.125,
source10,owner3,branchfull-gradient are FIXTURE VALUES, never a Kepler result.
The real dt limiting source remains UNCONFIRMED until a future necessary new
RHS uses the opt-in diagnostic. Do not replay a saved RHS merely to fill it.
H200 native full build now exists, not a validated physical continuation.
SSH build handle19696 and link handle5012 are CLOSED0.

## Preserved live pipeline and required next evidence

The unchanged producer launcher2166274/four GPU ranks and dispatcher2351718
were verified live by explicit ps during this turn. Current small progress
snapshot: 468792/1138432triangles (41.18%), prefix8.169x
point reduction. Projected total wallfromstart~5.2–5.4h, NOT remaining time.
A progress file alone is not liveness proof. No complete rule or next accepted
physical time is claimed.

Keep the existing H200-first single atomic-claim dispatch. On actual full ready
completion, run the already built independent whole-record audit (limit0),
inspect the actual paired stage1–4/endpoint/v3 restart and selected-cost result,
and adopt step5 only if its real gates pass. Saved highk1 must be reused.
If current paired accuracy/cost fails, follow the adopted strategic cost stop;
do not keep optimizing this fixed10.8M-DOF space or sweep dt/order/tolerances.
If it passes, the conditional continuous-map/thermodynamic/output contracts
still have to be closed before long production. This diagnostic changes none
of those gates. Shared global context and original HEAD8e25dfe/remoteLagEunha
were read/revalidated and were not modified.
