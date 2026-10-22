# Compact-space topology preflight, following measured production cost failure

2026-10-22. This is implementation progress under
`ASTRA_COMPACT_SPACE_AFTER_COST_FAILURE_20261022.md`, not physical evolution.
The long-term t6–10 heating/stability and longest-GIZMO objectives remain
active and unachieved.

## Adopted decision and preserved results

The actual positive-rule paired step5 passed accuracy and restart but failed
the adopted production cost gate. Its last accepted physical time remains
`0.00011263377801805954`; all four checkpoint shards are preserved. Read
`RESULT_ACTUAL_POSITIVE_TRANSITION_COST_GATE_20261022.md` and the retrospective
cell/output audit for that completed calculation. No unchanged-space long
continuation, new physical trajectory, or optional diagnostic was launched.

The existing Astra agent supplied one new bundled strategic decision:
retain the same source cells and protect the existing h/2 subdivision for
cells intersecting r<0.5; use parent fan P3 outside, with constrained shared
traces and two initial-velocity residual functions per connected component.
This changes the trial space, requiring a separately identified discretization
branch rather than silently projecting the evolved step5 field or repairing
its heat. Initial weak-operator restriction, curvature, physical accuracy,
memory and actual timing gates still have to be implemented and passed.

A separately appearing `ASTRA_CONSERVATIVE_COMPACT_REZONING_PLAN_20261022.md`
was explicitly disclaimed by the active Astra reviewer. Its author is not
established by this run. It was neither overwritten nor executed. The selected
plan is `ASTRA_COMPACT_SPACE_AFTER_COST_FAILURE_20261022.md`; do not launch
both alternatives. Existing reports/checkpoints in the shared ROOT were
preserved; the new count preflight writes only its own files.

## One completed native structural check

`inspect_compact_reference_topology_20261022.cpp` reads the original saved
canonical reference directly. It checks every stored group of four children,
source-cell identity, shared edge/midpoint identity, manifold connectivity,
and bitwise agreement with the original midpoint construction. It counts
the proposed space without evaluating physical point fields or forces.

| Structural quantity | Actual count |
| --- | ---: |
| Original source cells | 71,080 |
| Original fine triangles | 1,138,432 |
| Recovered parent triangles | 284,608 |
| Protected core source cells | 5,524 |
| Protected parent triangles | 22,168 |
| Active candidate triangles | 351,112 |
| Active vertices before trace constraints | 176,189 |
| Active edges before trace constraints | 527,540 |
| Coarse/fine interfaces | 240 |
| Fine trace slave variables | 1,200 |
| Polynomial scalar DOF after constraints | 1,581,181 |
| Connected components | 1 |
| Polynomial plus two residual scalar DOF | **1,581,183** |
| Reduction from original 10,816,441 scalar DOF | **6.8407 times** |
| Fixed candidate upper limit | 2,500,000 |

The five trace slaves at each interface are its midpoint vertex and four
fine-edge interior variables; the two coarse-edge interior variables remain
masters. This is a topology/constraint count, not an assembled injection,
bordered M factorization or owned/ghost partition measurement. The latter
must still be constructed. DOF reduction is not measured RK throughput.

Triangles ambiguous under inexpensive radius/diameter bounds were resolved
using exact dyadic polynomial predicates (804 such parent triangles), with
checked error-free binary64 expansions and FMA products. The code rejects
overflow or subnormal product terms instead of silently discarding them.
Core classification is fixed from the original reference, not fitted to
pressure results. No core radius or polynomial-order sweep was performed.

The native execution returned exit0 after 0.596283277 seconds. Its only first
compile attempt failed because Boost headers were absent; before any native
execution the dependency was replaced with local expansion predicates and
the rebuilt program completed once. See the execution record, build log,
source/binary/input SHA256 and topology JSON under
`material_compact_topology_preflight_v1`.

## New finding: topological nesting does not prove exact geometric nesting

The original construction rounds a long-double midpoint to binary64. Among
the recovered parent edges, 1,684 have a midpoint different from the exact
coordinate mean; 1,576 have a nonzero exact collinearity determinant.
The largest midpoint coordinate rounding observed in long-double arithmetic
is 5.5511151231257827e-17.

This is not an observed physical fold or a measured heating mechanism.
It does mean that successful connectivity recovery alone cannot justify
an exact polynomial nesting/identity-map claim. The actual coefficient
injection and shared traces must account for these stored coordinates and
measure their effect on initial x/v, volume, M, work and torque. The old
fixed mass and generating thermal context must not be quietly rebuilt to
hide a discrepancy. No additional enrichment, relaxed tolerance or geometry
repair is authorized by this count result.

## Next required implementation, without a new approval step

Construct the actual constrained injection and rank2 residual border from
the retained native basis; produce real owned/ghost and matrix storage
counts, restrict the authoritative stored M/initial force, and apply the
unchanged full initial curvature/physical gates. Reuse the positive-rule
and resident-DD work where the same-function restriction permits it.
Only a candidate passing these gates can enter the native full-disk branch.

The first actual steady costs must satisfy t1<=24h and t10<=7days under the
adopted plan, without assuming a larger timestep or free future guards.
Full t6/8/10 density/pressure/orbit/energy/rotation checks remain the target;
this structural PASS supplies no evidence of their success. Material-map
lifetime, conservative future remapping, W0/Voronoi reduction, shocks,
matched boundaries and the longest normalized GIZMO comparison remain open.
