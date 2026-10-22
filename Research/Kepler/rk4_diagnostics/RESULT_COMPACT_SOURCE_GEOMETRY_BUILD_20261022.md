# Full-source import and geometry components of compact transaction

2026-10-22. Progress under the one adopted
`ASTRA_CONSERVATIVE_COMPACT_REZONING_PLAN_20261022.md` bundle. Native code was
implemented and compiled, but no new physical calculation or separate small
guard experiment was launched. Current accepted state remains step5 at
t=0.00011263377801805954; all completed source data and checkpoints persist.

## Implemented source connection

`material_compact_source.h` connects the real full 71,080-cell / 1,138,432-
triangle / 10,816,441-scalar-DOF source to the exact v3 checkpoint reader.
It loads the original saved basis, preserves generating coefficients and
immutable thermal corrections, and restores actual evolved coefficients.
It retains the original fixed-M path; no mass or initial RHS is rebuilt.

The positive-rule generation comes from the four completed independent WHOLE
audit streams. The code parses uint64 stream IDs exactly, without conversion
through double, then reproduces the original native generation hash. Exact
checkpoint binding includes the original generating manifest and the same
positive-rule directory/generation. The complete transaction must verify the
prelaunch SHA256 manifest covering these actual input files.

The original owned/ghost MPI operator checks all four position/velocity
components against their owners without replacing mismatched checkpoint bits.
The original reference receipt is bound to full cell/triangle/vertex/edge
counts, square topology/area and actual gamma/P0/eps. Actual native gamma is
1.6666669845581055; do not replace it with an assumed exact5/3.

`Model::layout` must be destroyed before MPI_Finalize, as in the original
Ownership API. The eventual full transaction controller must use a scoped
model lifetime and collective abort/rejection for rank-local import failures.

## Implemented current-source geometry checks

`material_compact_source_geometry.h` routes canonical edge descriptors to
their original owner, verifies exact function parameters and global polynomial/
enrichment DOF identity across ranks, and checks oriented interior/exterior
incidence. It does not infer exact shared traces from sampled agreement or
merge slightly different coordinates.

For every actual source triangle it binds the entire enriched current map,
including edge and bubble terms, to the accepted coefficient hash/time/epoch.
It reuses the existing binary128 real-function enclosures. The tight Taylor
bound is tried first; if its analytic domain is not available, the existing
general real enclosure is the one finite fallback. There is no physical mesh
change, precision/order search or unbounded proof refinement.

The sufficient CURRENT-state condition is global
sup||Dchi-I||_infinity<1 on the already constructed convex reference square,
with exact shared canonical traces. Segment integration then gives a global
lower Lipschitz bound and excludes multiple preimages. This is conditional on
the original reference construction and the existing <=8ulp engineering libm
premise, not a new unconditional proof of CPU/GPU elementary functions.
An unresolved triangle makes the global bounds unbounded and the certificate
false; it is not reported as actual folding. The resulting source bounds are
retained in memory for this exact epoch rather than treated as reusable for a
different state. No complete target-domain or corrected-thermal proof follows
from this source-only condition.

## Missing endpoint context identified

Step5 endpoint stores physical cell aggregates and low/high differences, but
not the absolute full per-source endpoint correction context. The saved stage4
context belongs to its RK stage state and is not interchangeable with the
accepted endpoint. Code explicitly leaves `current_context_bound=false` and
blocks physical transfer until the actual accepted-geometry/rule context is
derived. It must not be fitted from mean U, copied from stage4, or reset to P0.
Only missing current context may be reconstructed; accepted intervals and the
already consumed step4 high k1 are not replayed.

The current exterior is a curved material boundary. Actual target power-cell
domain coverage must be established; clipping to the initial square is not a
proved replacement. Neither this source certificate nor a receipt alone is
an implemented compact remap.

## Actual build evidence

`material_compact_rezone.cpp` exposes full-source import and geometry component
entry points. It deliberately has no standalone partial-experiment main; the
complete transfer/EOS/conditional-step controller remains to be assembled.

On authorized lageunha, `material_compact_rezone_build_v3` completed exit0:

- Native MPI/GNU17/OpenMP source object PASS (`compact_native.o`).
- Ada sm89 and H200 sm90 CUDA field-sampler objects PASS.
- Frozen/root sources match; source/object integrity reports PASS.
- 26 compiler project dependencies recorded in `dependencies.sha256`.
- Sessions59650 and80230 CLOSED0. No tool session or calculation remains live.

`build_contract.json` records actual import execution=false, source proof
execution=false, thermal-context reconstruction=false, target coverage=false,
projection=false, full application link=false, compact commit=false and new
step=false. This is build evidence, not physical or runtime validation.

Remaining implementation is actual current thermal context, canonical current
target geometry and weight-linked mass, conservative overlap packets, positive
rho/e fit, coupled shared p/L/K projection, field/entropy error ledgers and
compatible EOS/force/mass/RK4/restart. Run the complete bounded transaction
once assembled, then apply actual accuracy and amortized-cost gates.
Original user source/results and shared global context were not modified.
t6/8/10 and the longest normalized GIZMO comparison remain incomplete.
