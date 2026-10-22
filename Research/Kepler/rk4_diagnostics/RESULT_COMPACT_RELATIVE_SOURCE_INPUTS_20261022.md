# Relative geometry precision and actual-source target input adapters

2026-10-22. Progress under the same compact transaction. Driver codes; no new
Astra call, actual source guard experiment or physical simulation was launched.

`material_compact_arithmetic.h` reuses the original CUDA two-double arithmetic
on host/device. `relative_trace` evaluates the frozen source curve relative
to its endpoint line using P3 curvature coefficients and the original scalar
residual functions. It avoids absolute Gamma/chord cancellation before small
mode normalization. The same small interpolation values now feed patch and
canonical boundary bases. Endpoint differences remain recorded: no snapping,
mode removal, exterior flattening or domain adjustment is introduced.

`material_compact_relative_v1_ada`:387 samples, coordinates about12, mode
scale1.45519e-11, independent binary128 reference, PASS(native/runner0),0.285s.
Old normalized error4.04611e-4 -> new3.39809e-15; mode-value error8.77989e-18;
CPU/Ada difference7.10543e-15. Endpoint residual1.23862e-15 was retained.
This verifies the tested finite evaluations, not actual subinterval event
identity/domain coverage. Root/event enclosures and exact shared charts remain
required before physical transfer.

`material_compact_predicate.h` supplies an error-free expansion/FMA orientation
sign, including exact zero;1280 independent binary128 checks pass. Its current
engineering input range covers the ORIGINAL reference square (abs<=16,
nonzero abs>=2^-128). It is not a physical bound on long-time positions and
must not be blindly reused outside that numerical range for later epochs.

`material_compact_source_target.h` implements actual-source target INPUTS:

- Reads original native generator coordinates from the exact reference shard;
  centered coordinates are checked algebraically, not advanced analytically.
- Locates each generator in its own original reference cell using exact
  orientation and deterministic minimum source-triangle ID ties.
- Vertex/edge/interior paths evaluate the actual accepted coefficient field,
  with canonical source-function enclosures. No IC orbit/velocity reset.
- Preserves actual accepted source-functional mass separately from legacy
  native mass and any future screened-law reference mass/mu.
- Extracts each owned exterior curve once, with canonical vertex/function
  identities, current frozen coefficients and whole-source bounding boxes.
- Requires the actual source geometry/current context and coefficient epoch;
  checks global71080 generators and the original exterior-edge count.

This input extraction is implemented/compiled, NOT run on the accepted source.
Original raw geometry shards must be included in the eventual input manifest.
The caller must provide collective failure handling around rank-local I/O.

New source-location fixturev1 failed because its coarse triangles crossed
the origin and violated the UNCHANGED original radial-Taylor basis domain
abs(L+Q)<0.75. Failure/source/binary remain preserved. Only fixture coordinates
were moved to a valid domain; production basis/physical parameters unchanged.
Fixturev2 passes6 vertex/edge/interior/near-edge queries using an evolved
affine/shear field, including shared-face ties. Field difference8.882e-16;
native/canonical enclosure error4.996e-16; native/runner0,wall1.955s.
No actual source import or target_source_inputs invocation occurred.

The NEW-header RHS regression `material_compact_operator_v2_ada_mpi2` passes
with enriched curves and original resident CUDA/MPI mass solve. Eidentity/U
4.66848e-13 after mass solve, true mass residual5.70068e-20; native/runner0,
wall2.849s. This does not extrapolate production cost or long-time stability.

Buildv8 passes nativeMPI/OpenMP, rebuilt geometry+thermal AdaSM89/H200SM90
objects and original GPU component links.38 project dependencies and source,
object, input, fixture integrity/root equality checked. Sessions71074,58198,
15320,29444,32211,12318,51964,4440 closed0;51244 closed1(failedfixturev1).
All handles terminal; no new simulation/job/dispatcher/monitor.

Accepted state STILL step5/t0.00011263377801805954 at
material_full_positive_transition_v2_ada_mpi4/run/continuous_step5/accepted_r.
Original repository HEAD8e25dfe, old kernels/checkpoints/results and shared
global context unchanged. No source proof/context/actual target extraction,
weight-law solution, conservative projection, compact commit or new step.

Next assemble the small global transported-site/exterior registry, actual
curved power partition, finite immutable W-linked mass law and positive actual
overlaps, rho/e/shared p/L/K projections with independent field/entropy/error
budgets, chart/rank/positivity guards, restart and bounded full controller.
Only a passing actual transfer permits the next RK4 step and unchanged
amortized t1<=24h gate including collar/remap/guard. t6/8/10 and the longest
normalized GIZMO test remain open; artificial shear heating is unresolved.
