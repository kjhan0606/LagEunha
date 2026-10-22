# Full native reference connection — 2026-10-22

PROGRESS after the actual native RK4/restart bundle. This is full-domain input
and reference construction, not full-domain hydrodynamic evolution. The user's
original t6/8/10 and GIZMO longest-duration objective remains ACTIVE/UNPROVED.

## Authoritative evidence

- ASTRA_FULL_DISK_CONNECTION_PLAN_20261022.md: next major plan delivered and read;
  do not call Astra again for routine implementation of this bundle.
- shared_full_inputs_v1/contract.json: inspected all71080 stored native GPU
  polygons and their original generator metadata, exact match to the preserved
  fine_integrated IC SHA46253652acfdd282904f77b1ac230708aa29a315aca7dcc96c89cb07342d2437.
  No new tessellation, physical field evaluation or IC generation.
- Exact coordinate incidence:284608 nonzero directed polygon edges,
  141792 paired material interfaces,1024 external box edges, ZERO unpaired
  interior or nonmanifold edges.65536 exactly zero-length native edges omitted;
  no tolerance/near-coordinate merging. Domain area576. Interior support IDs
  were not present in this old file; exact matching plus actual manifold/domain
  checks sufficed for this particular retained structured/refined geometry.
- Native C++ build_full_reference_v1.cpp/material_full_reference_geometry.h
  produced shared_full_reference_v1/full_reference.bin,47259664bytes,
  SHA518e06c1cbef7106fedce680c4b5a4a7be265a48e8eee60a5df244f9d12f69f7.
  Reuses retained native geometry and the same one-level4-way fan refinement.
  1138432triangles,570241vertices,1708672edges,2048external refined edges.
  Positive triangle area and exact fixed-rule geometry checks PASS; total
  area576.000000000000004385,minimum triangle area4.283642519641983e-6.
  Native construction/counting time0.4191s on authorized lageunha,OMP4.
- shared_full_reference_v1/contract.json counts the actual fixed64x28 rule:
  logical queries6120597504 and positive points2040256512. These are geometry
  COUNTS; no6.12b query loop or2.04b field integration has been executed.
- shared_full_reference_adoption_v1/contract.json: actual binary artifact import
  and compact full-reference packet PASS. The same canonical Factory topology
  C++ template instantiates for the compact packet; the global full Factory/mass
  has NOT been constructed. Root material_shared_trace_factory.h's topology
  signature is now templated, leaving the previously executed frozen source
  and successful old outputs untouched. New compact record56bytes; all full
  reference records63752192bytes. The packet carries geometry and cell IDs,
  not legacy per-triangle velocity maps or a fake copied559fragment state.

## Concrete resource implications

P3 scalar DOFs5126017; enriched scalar DOFs between9111553 and10820225 pending
actual radial/nonradial edge rank assignment. One global raw binary80 ST18
mass block array would be5901631488bytes before transformed blocks/factors/CSR.
Native sizeof STPoint64,GPOutput504. Fully materializing points and two outputs
at all positive nodes would require2187154980864bytes (~2.0TiB), and the literal
old CPU audit path is inappropriate. A65536-point window uses70254592bytes for
these three arrays. Production must stream reference quadrature and use GPU
reductions plus distributed owned/ghost mass/force operations, as planned.
These are sizing facts, not measured full hydrodynamic seconds/step.

## Original source identity and isolated implementation directory

Direct tools are currently running on LagEunha. The shared original src/.git
references /tmp/lageunha-user-master/.git/worktrees/src3, which exists on syntax
but not lageunha. SSH syntax successfully reconfirmed original toplevel,
HEAD8e25dfe0c2e41a7374c11e2e3102da2a8b028625 and origin
 git@github.com:kjhan0606/LagEunha.git. No repair/reset is needed or performed.

native_full_material_20261022/ is an isolated39-file source closure from the
actual original native program and its Exam GPU/source includes, with exact
original_source.sha256 and PROJECT_IDENTITY.json. Old original sources,
objects, executables and successful results remain intact. This directory has
not yet been linked and does NOT yet implement a material evolution mode.
NEXT_IMPLEMENTATION.md specifies the actual remaining integration work.

## Immediate remaining implementation

1. Canonical full DOF/global IDs, spatial triangle/DOF owners and shared halos;
   distributed fixed mass/local preconditioner and accurate true residual solve.
2. Persistent streamed CUDA geometry/thermal/mean-pressure/force reductions;
   immutable generating entropy function, same-stage geometric pressure work,
   sparse CPU shadow on first full RHS only. No full per-point result transfers.
3. Completed C++ material mode connected to original Exam/exam.c opt-in branch
   and original MPI lifecycle; runtime/restart and GPU ownership manifests.
4. Actual continuous full-state initialization and first step, preserved as the
   first part of the trajectory. Explicit Pext=P0 moving boundary is this first
   research configuration; it is a boundary change from the periodic historical
   baseline and cannot be reported as an exact baseline/paper match.
5. The same accepted trajectory toward t1/24h first batch then t6/8/10, with
   Eulerian density/heating/orbit/E/L diagnostics and map injectivity controls.

No successful local3step, native polygons or23.339m transfer queries should be
rerun. Current tool/process handles for both full input programs are CLOSED,
exit0. No simulation is live. Shock entropy, conservative remesh/context
re-expression, dynamic Voronoi W0 reduction and GIZMO longevity remain unproved.
