# Conservative overlap packet contract — 2026-10-22

## v5 update: bounded-memory packet stream

Refactored the common validation and moment ledger into `Validator` and
`StreamAccumulator`. An overlap producer can append one point/physical sample
at a time, across arbitrary chunks, while enforcing one global ordered-ID,
point-count and area-error budget. Finalization checks the independent area
marginals and source/target m,q,Q,p,L,K,U,Phi totals. The old vector functions
now delegate to the same stream implementation. The native degree-one/fixed-M
entrypoint can construct this accumulator from the actual target registry and
chart, without constructing a complete packet vector.

Memory now scales with source-triangle and target-cell/patch ledgers and ID
lists, rather than overlap point count. Source-field evaluation, chunk
producer, deterministic MPI reduction, and projection remain to be connected.
No claim of full-state memory feasibility follows from the four-point fixture.

The synthetic two-chunk fixture compared buffered/streamed moments and
rejected duplicate IDs across a boundary, a global point-budget breach,
premature finalization, and double finalization. A10 job 434721 completed
exit0 (3516 KiB max RSS); the native rezone object compiled. Frozen files and
receipt: `material_compact_overlap_stream_v5/build_contract_v5.json`. Accepted
step5 and its physical time remain unchanged; no actual transfer or step ran.

## v4 update: original identities and shared moment ledger

The original v2 validator incorrectly treated `Point.target_cell` and
`Point.source_triangle` as dense array indices. The target registry instead
uses original, potentially noncontiguous cell IDs. The validator now accepts
independent sorted source-triangle and target-cell ID lists plus a target
patch-to-cell map. It checks membership, patch ownership, and ordered-list
hashes against the transaction seal. The hashes guard accidental mixups; they
are not cryptographic source-file verification. It also limits each point's
coordinate uncertainty and the packet plus independent marginal area errors
under the supplied finite budgets.

`material_compact_moments.h` now consumes the *same validated packet and source
field stream* once to assemble outgoing source-triangle and incoming
target-cell m, qx/qy, Qxx/Qxy/Qyy, px/py, L, K, U, and Phi. Both directions use
the source physical coordinate and the identical per-packet contribution. It
requires matching packet IDs and transaction seals, finite strictly positive
physical density/internal-energy density, and algebraic source/target totals
within long-double roundoff. It does not certify continuous physical
quadrature, source-field evaluation, or the future target P3 projection.
The native rezone component now exposes a degree-one/fixed-M-law-gated entry
point that derives target original cell IDs and patch ownership from the
actual registry/chart before invoking that ledger. It compiled as an object;
it has no transaction caller yet.

The synthetic fixture uses original source IDs 101/103 and target cell IDs
10/20. It rejects unregistered/changed IDs, a chart assigned to another cell,
stale field provenance, incorrect sample IDs, negative rho/e, and inflated
coordinate/area errors. A10 job 434703 passed the initial ID/moment version;
434704 passed the strengthened sealed-ID/error-budget version (exit 0, 3516
KiB max RSS). The final ledger reports area 1, mass 2.5, heat 5, kinetic
energy 66.25, maximum source/target reduction gap 8.67362e-19. Frozen sources,
binary, job script and receipt are under
`material_compact_overlap_packets_v4/build_contract_v4.json`.

The original accepted state remains step5 at t=0.00011263377801805954. No
actual source state was imported, no actual overlap or physical fields were
constructed, and no projection, restart commit, or new RK4 step ran. Next
priority is the single bounded actual source/target overlap and positive
field-transfer transaction, including native chart rounding and physical
quadrature-error bounds. The synthetic ledger is only a prerequisite. The
old vector API holds a complete packet stream in memory; v5 supplies a bounded
chunked alternative. Production still requires deterministic cross-rank
reduction before any full-state run.

For that evaluator, the accepted CUDA RHS uses reference density from
`gfs_enrich_reference`, corrected current volume factor
`j=J_current+delta_cell`, and corrected generating factor
`j0=J_generating+generating_delta_cell`. Thus packet physical density must be
`rho_reference/j`; its physical heat density must be
`[P0*j0/(gamma-1)*(j/j0)^(1-gamma)]/j`. The source velocity comes from the
current accepted map, and softened potential is
`-1/sqrt(x*x+y*y+eps*eps)`. The actual evaluator must bind these to the
accepted source state and certify the transformation/error bounds. Treating
`rho_reference` as physical rho or using bare J would break the transfer.

## Historical v2 scope

Added `material_compact_overlap.h` as the common data and validation contract
for a future source/target conservative transfer. One packet represents one
positive physical-area quadrature point and stores both source-triangle and
target-chart coordinates, their physical-position enclosures, owner IDs, and
the immutable source/target/weight/fixed-mass-law seals.

The validator checks finite positive area weights, finite nonnegative error
bounds, closed-triangle barycentrics, ordered unique packet IDs, source/target
owner ranges, overlap of the two physical-point enclosures, a finite global
quadrature-error budget, exact agreement with independently supplied source,
target, weight and fixed-mass-law seals, and independent source-triangle and
target-patch area marginals. At v2 the packet type only specified a shared
input stream. The v4 algebraic moment ledger above now consumes one stream;
actual field projection operators remain unimplemented.

## Evidence and limits

`check_compact_overlap_packets.cpp` exercises a four-packet synthetic partition
and rejection cases for negative weights, duplicate IDs, inconsistent inverse
locations, missing marginal area, stale source epoch, and out-of-domain target
coordinates/patches. The first version passed on A10 in Slurm job `434695`;
after that run, provenance equality and independent area-error inputs were
added. The strengthened version passed on A10 in job `434696` (exit 0, 3516 KiB
max RSS). Its source, binary and job receipts are recorded in
`material_compact_overlap_packets_v2/build_contract_v2.json`.

The packet validator does **not** construct source/target intersections, invert
the actual charts, evaluate physical fields, perform density/energy or shared
momentum projection, or commit a restart. Area and marginal agreement is only a
discrete measure check: `continuous_physical_coverage` is deliberately always
false. The chart disk/Jacobian/degree gates remain separate prerequisites.

Accepted state remains step5 / `t=0.00011263377801805954`; no accepted-source
import, projection, commit, or new physical step was run. The active fixed-mass
law remains `mref / saved Vsource`, coefficient 1, immutable through the W
solve. Selected point-functional mass must not be substituted for `mref`.

## Next implementation gate

Implement the bounded actual overlap constructor against the exact source
triangle maps and the sealed target charts. It must emit the same positive
packet stream for both directions, bind each source/target marginal to its
independent positive geometry integral, retain interval and quadrature error
bounds, and reject unresolved intersections. Then implement positive
reference-measure rho/e transfer and the coupled shared p/L/K projection.
Do not invoke the actual accepted source until these pieces, compact
commit/restart, and the full transaction cost estimate are connected.
