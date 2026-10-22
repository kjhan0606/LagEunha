# Compact power cells and curved-domain boundary topology

2026-10-22. Continued under the existing authorization. Driver implemented
the geometry primitives; no new Astra call or actual-source experiment.
Shared context read fully; original repository identity rechecked:
`Kepler_stencil5x5_Nx256_A100_call33720_20261020/src`, HEAD `8e25dfe`.
Original sources, CUDA kernels, completed outputs and accepted state retained.

## Implemented

- `material_compact_target_registry.h`: fieldwise MPI1/2/4 gathering of small
  transported-site/frozen-exterior geometry. Binary80 mass is transported
  losslessly as two doubles; native struct padding is not used. Exact epoch,
  time bits and four source hashes bind the packets. Unique cell/edge owners,
  exact shared accepted vertex x/v bits and one oriented exterior cycle are
  checked. This replicates small geometry, not field DOFs or a dense mass solve.
- `material_compact_power_curve.h`: canonical sorted site-pair power planes,
  with actual generator enclosures and finite weights. Evaluates the SAME
  canonical source trace using the original functions, not fitted boundary
  curves. Positive/negative range exclusion, strict derivative/bracket and
  interval Newton isolate transverse roots. Root identities retain plane,
  source edge and ordinal. Coincident/endpoint/uncertain cases remain
  unresolved; no endpoint snap or approximate zero merge is performed.
- `material_compact_power_cell.h`: bounded plane power cells using a rectangle
  enclosing the current exterior. This is NOT old-square clipping. A spatial
  tree omits nodes only after a strict bound on every omitted power inequality.
  Canonical interior vertex identity is the sorted site triple: all three
  incident cells evaluate the SAME two defining planes. Bounding vertices and
  two-line/root definitions remain implicit, with Q128 enclosures. Four-way
  concurrence is unresolved, not perturbed. Work/corner budgets are explicit.
- `material_compact_domain.h`: builds the oriented boundary graph of
  `Omega* intersect PowerCell_i`. Retains straight-edge portions inside Omega
  and original-curve portions inside the cell. Point classification uses a
  fixed ray and complete certified curve events; degeneracies reject rather
  than trigger ray/parameter scans. Shared root identities join exterior
  pieces without dropping strips. Graphs may have multiple components and
  valid two-arc curved caps. The interior-generator shortcut is disabled
  until bound to the actual global source geometry certificate.

These are geometry components, not a complete transaction executable. New
native wrappers are linked into `material_compact_rezone.cpp`. Root and corner
enclosures must not be converted to new exact rounded chart geometry. Positive
chart/Jacobian and coordinate-mode rank guards remain necessary.

## Finite checks on authorized lageunha

| Fixture | Evidence | Native/runner | Wall |
| --- | --- | --- | --- |
| target_registry_v1 MPI1/2/4 | mass low bits and all curve fields retained; empty owner rank; 9 malformed cases rejected each | 0/0 | 1.840/1.880/1.936 s |
| power_curve_v1 MPI1 | 9 cases; three crossings; weight sign/gauge and pair reversal; 44 boxes; largest root width 1.81527e-17 | 0/0 | 1.781 s |
| power_cell_v2 OpenMP2/MPI1 | 25 cells x2 weight choices; 882 plane tests vs1200 exhaustive; 39 excluded nodes; exhaustive inequalities and96 shared vertices agree; area enclosure width1.61778e-31 | 0/0 | 1.905 s |
| domain_v1 MPI1 | 5 curved-domain cells;24 arcs,8 curved,8 shared straight; full exterior joins; two-arc cap; two disconnected components; broken closure/unbound proof rejected | 0/0 | 1.983 s |

The power-cell construction portion took0.0234913s on2 OpenMP threads for the
small fixture. None of these timings extrapolates the actual full transaction.
The curved-domain fixture used68 root calls/242 boxes/32 winding calls. These
are synthetic geometries, not accepted source import, physical validation,
restart, new time evolution or a proof of actual step5 target coverage.

`material_compact_rezone_build_v9` through `_v12` all have passing native
MPI/OpenMP objects and original Ada/H200 GPU component links. Existing v8
geometry/thermo and v3 sampler CUDA objects are unchanged/reused, along with
the unchanged original native CUDA object. Latest v12 records39 project
dependencies; source/object/dependency/retained-input SHA checks pass and all
frozen compact sources equal the current root. Earlier frozen versions retain
their own scope, before subsequent components were added.

Preserved engineering failures: v9's first strict C++17 compile rejected an
ORIGINAL Q numeric suffix. GNU++17 attempt2 passed without editing that source.
The serial power_cell_v1 result JSON reports PASS, but its shell runner had an
unmatched quote AFTER execution and exited2; its exact native-return receipt
was not captured. Files and `runner_receipt.json` are preserved. The separate
new OpenMP2 fixture v2 has complete measured native/runner0 receipts. No actual
source/physical interval was replayed for either issue.

## Follow-up: retain canonical generator precision (v13)

The first adapter retained a native double position with its cache error.
Using only that ball in power geometry imposes a binary64 rounding floor on
new root/corner enclosures. `Site` now ALSO stores native-high/canonical-low
and a tail enclosure for both position and velocity, with an explicit
`canonical_bound` flag. These come from the SAME original sf source evaluation;
the native high values and original cache balls remain. Vertex/edge/interior
paths do not reset an orbit or velocity. Sites without this provenance retain
their original uncertainty; a missing flag never implies a tighter ball.

The fieldwise wire now has5 integer tags/20 doubles per site. Power planes,
spatial exclusions and the certified-generator interior shortcut use the
refined canonical position. This is a metadata/API refinement, not an accepted
checkpoint ABI change. Curve/native GPU kernels and evolution remain unchanged.

New-header regressions on lageunha all pass native/runner0:

- target_registry_v2 MPI4: explicit precision flags, low words and tail bounds
  survive, alongside all prior mass/curve/closure/9-malformed-case checks.
- source_location_v3 MPI1: all6 source chart cases return canonical precision;
  position/velocity tail errors <1e-25; native fields still agree8.882e-16.
- power_curve_v2 MPI1:10 cases. Original uncertain q remains unresolved;
  an explicitly supplied independently refined SAME q can certify the root.
- power_cell_v3 OpenMP2/MPI1: refined q through tree/planes,882vs1200tests,
 96shared vertices and exhaustive checks PASS; area width2.60890e-28.
- domain_v2 MPI1: refined q through all5 cells and full curved exterior joins,
  two-arc cap and two-component case PASS, unchanged24arc topology.

Latest frozen build is `material_compact_rezone_build_v13`: native MPI/OpenMP,
original Ada/H200 component links and retained-input/source/object/dependency
integrity pass;39 project dependencies and current-root equality are recorded.
All original CUDA objects are reused unchanged. This follow-up still performs
NO actual source import/context/extraction, actual target/transfer or new step.

After the connected topology bundle, the existing Astra planning/evaluation
agent was asked ONCE to review false acceptance conditions and fix the next
finite chart strategy plus immutable mass-law parameters. Driver remains the
coder. No new agent or primitive-by-primitive review was requested. The review
is now complete; its concrete decision is described below.

## Astra decision and chart connectivity (v14)

Read `ASTRA_COMPACT_TOPOLOGY_AND_CHART_DECISION_20261022.md`. Astra found no
obvious direction reversal or unresolved-to-success branch in the reviewed
geometry primitives, without elevating this to an actual-state proof. It
adopted mref=previous accepted fixed-M cell mass, mu=mref/Vsource with coefficient
1, sealed before W Newton and never refitted to W0 overlap mass. Vsource is the
saved canonical physical source area. A manifest and actual law solve are NOT
implemented/executed merely because this parameter choice is adopted.

Driver added `material_compact_chart_topology.h` following the ONE selected
candidate: each curved arc gets ONE canonical parameter midpoint; certify a
simple positive chord polygon; select ears in fixed canonical-node order;
only ears with curved edges get a three-way centroid fan. Bulk uses the compact
centroid fan. Each collar patch has at most one original curved edge. Midpoints
and centroids keep parent function identities/enclosures, not new exact
rounded coordinates. These reference nodes do not steer source generators.
Holes, crossings, collinear/uncertain predicates and exhausted finite budgets
reject. Native curved-J/coverage/rank certification remains unimplemented.
Patch.component is a LOCAL graph-cycle index, not an assigned global FEM
component. Assemble shared-edge connectivity before mode IDs: at most2 scalar
geometry modes per GLOBAL connected reference component, not2 per cell or
per local cut cycle. The current physical Omega is one connected component.

The actual-domain and actual-chart wrappers now require actual source geometry
binding even if the interior-anchor shortcut is unused. Each spatial trial
owns an immutable W copy and lossless scalar weight seal. Bounded cells,
domain graphs and chart candidates carry the source/weight/epoch provenance;
stale source or W data cannot pass the checked next entrypoint.

New `material_compact_chart_topology_v1_mpi1` passes native/runner0: the5 curved
cells plus the two-arc cap produce51 candidate patches,18 curved patches and9
midpoints. Every patch has <=1 curved edge and internal fan edges appear twice
with opposite orientation. Cap needs only its one prescribed midpoint, giving
3 positive chord patches. Open graphs, unbound source, stale source/W and a
self-crossing chord polygon reject. The nonconvex-domain test still retains2
boundary components; this does NOT claim that those components have passed
the later chart-J/coverage/rank proof.

`material_compact_power_cell_v4_omp2_mpi1` passes with immutable weight-copy
checks, original882vs1200 exhaustive inequalities and all prior geometry tests.
Latest frozen buildv14 passes native MPI/OpenMP and unchanged original
Ada/H200 component links,40 dependencies and source/object/input/root equality.
No CUDA kernel changed or was reexecuted for these CPU geometry checks.

IMPORTANT mref provenance discovery: the old accepted endpoint writer
`material_full_evolution_host.h` physical_cells CSV mass column is p[3], the
selected quadrature point mass; it is NOT a direct fixed-M contraction. Do not
silently label EndpointBaseline.mass or Site.source_mass as the adopted mref.
The full transaction must read the unchanged `Model.fixed_mass_path` payload
(magic0x4d464d4153533031/version1, nt*ST_MAX*ST_MAX binary80 entries). Per source
triangle, constant1 has first10 P3 coefficients1 and residual coefficients0,
so fixed-M mass is sum of the first10x10 entries; sum these by source cell.
Keep the endpoint/point-mass discrepancy in its existing numerical budget.
The saved canonical area remains the correct Vsource. This reader/manifest
binding is a next implementation item, not an actual source calculation here.

## Physical state and remaining work

Authoritative accepted state STILL step5 / t=0.00011263377801805954:
`material_full_positive_transition_v2_ada_mpi4/run/continuous_step5/accepted_r`.
No actual source guard, endpoint reconstruction, generator extraction,
target-domain construction, screened-law solve, overlap/projection, compact
commit or new RK4 step has run. No new simulation/Slurm job/monitor was launched.

Next: construct actual component charts from the implicit boundary events,
certify same-domain coverage/J>0 and geometry-mode rank, build positive actual
source-target overlaps, seal finite immutable W-linked mass-law parameters,
fit positive rho/e and shared p/L/K with independent field/entropy budgets,
then restart/full controller and a conditional actual native step. All source
triangle interfaces in the overlap are their current curved functions too;
exterior topology alone does not supply conservative overlap packets.

The ONE bounded full transaction remains <=24h, with actual amortized t1<=24h
including guard/remap/collar costs. No tiny actual-source guard experiment,
old k1 reuse, completed interval/preprocessing replay or old10.8M-space long
production. Kepler Re=infinity, no steering/heat repair/spin reservoir. t6/8/10
and the longest normalized GIZMO test remain open; shear heating is unresolved.
