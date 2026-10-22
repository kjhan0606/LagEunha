# Shared canonical function, native cache and frozen heat: connection results

2026-10-21. Continuation of the active long Kepler objective. Previous goal
turn was PROGRESS (compatible geometric pressure/heat and actual nu
curvature). This turn is also PROGRESS: native function/cache/thermal and
real inverse contracts are implemented and verified. Entire target68526
integration and actual long evolution remain UNACHIEVED.

## State and implemented native APIs

Same SAVED finite shared_trace_rhs_v1_mpi4/shared_state.bin,24 physical
cells/384 source triangles/3797 scalar shared DOFs. Original coefficients,
reference mesh, mass inputs, quadrature, degree, precision and physical
criteria are retained. No trajectory, new acceleration, old-packet
reintegration, state fit, averaging, heat reset or steering.

material_shared_function.h compiles the canonical mathematical reference:
exact affine barycentrics and rational P3 interpolation nodes, analytic
v0=(-Y,X)*(X^2+Y^2+eps^2)^(-3/4), canonical edge remainder r=v0-I3(v0),
E=s^4 r(lambda_B/s), zero-trace interior residual and the stored invertible
basis normalization. The SAME stored coefficients are condensed
algebraically into P3 + B*v0 + sum(C_e*E_e); no sampled fit or replacement.
Native rounded qnode/gradient caches are evaluation approximations, not
independent exact local functions. Their errors are explicitly enclosed.

material_shared_native_ball.h mirrors the actual native radial series,
series stop, shared STJet basis operation graph, cache constants, gradients,
coordinate/lambda inputs, normalization, coefficient sum, J and inverse.
All13 native field values must agree numerically with the actual native
helper at every original saved quadrature point. Independent Q128
canonical function values must be inside native forward bounds. IEEE
nearest arithmetic and <=8ulp elementary-function models are explicit
engineering assumptions, not formal CPU/GPU libm proofs.

FrozenThermal owns a CONST deep copy of the generating geometry/state,
P0/gamma/epoch. For THIS saved adiabatic control only, immutable off-node
physical energy per reference volume is

    u_ref(X)=P0/(gamma-1)*J_GENERATING(X)^(1-gamma).

The generating state is not the later current geometry. Original physical
U is READ from the protected frozen_physical_rank*.bin point tables; the
recipe is checked against those original values. It is not used to
regenerate current simulation heat. thermal_generating_state.bin and
thermal_function.json preserve its exact identity. General evolved
thermal reconstruction/restart is still incomplete. Geometric push-forward
rho0/J versus EOS rho0/(J+delta) remains a required transfer error term;
no off-node density or thermal-rate integral proof is claimed here.

## Full original native field/thermal verification

shared_function_v2_mpi4, CPU/OpenMP4×MPI4 on authorized lageunha:

| Quantity | Result |
|---|---:|
| Original positive points checked |817152 /817152|
| Sources checked |384 /384|
| Native13-field operation-graph value mismatches |0|
| Canonical function outside native+oracle enclosure |0|
| Frozen physical U operation-graph mismatches |0|
| Largest observed field error /forward bound |.99093760779744|
| Max source physical U same-rule relative difference |1.00205973814e-16|
| Max source forward U bound /source U |2.03895827926e-14|
| Positive off-original-node frozen-heat queries |98304|

No observed-error multiplier enlarges a bound. The <=8ulp engineering
model is inherited, not calibrated to these results. U comparisons are
same-rule comparisons, not analytic quadrature-remainder certificates.
No new CUDA transfer kernel or GPU error certificate is claimed; previous
native CUDA physics remains unchanged, and this reference-function stage
runs in C++/MPI/OpenMP with Q128 software arithmetic.

The first coarse real source bound split small v0-I3(v0) into two large
terms. It gave max L71.7563782 and127 inverse failures, despite all actual
field/thermal checks passing. Native run correctly exited2; its full
witness and output are retained. v1 was a compile-only failure from a C++
multiple-auto declaration of different scalar types; v2 fixed it without
changing physics or rerunning a trajectory.

## Same-function global derivative bound and real inverse

material_shared_real_bound.h retains the small Taylor residual
q=v0-T3(v0) rather than losing cancellation:

    v0-I3(v0)=q-I3(q).

Its absolute q/gradient bounds follow the original binomial series with
|s|<.75, coefficient tail bound c4*|s|^4/(1-|s|) and analytic derivative
bound. The P3 polynomial and interpolation gradients use an exact nodal
Bernstein control-net conversion and its derivative convex hull. The edge
extension derivative uses

    grad E=s^3*[4 grad(s) r + (grad(lambda_B)-t grad(s))*r'],

with s,t in[0,1]. It remains bounded at the opposite vertex. This is a
REAL triangle derivative proof, not holomorphic continuation across that
vertex. Every upper-bound operation is rounded outward in Q128 under the
stated arithmetic model. No source state, guard, mesh, degree, quadrature,
precision, or contraction criterion was changed.

shared_real_bound_v1_mpi{1,4} both exit0:

| Quantity | Result |
|---|---:|
| Sources passing whole real source L<1 |384 /384|
| Maximum whole real source L |.0210777156723553|
| Newly computed inside inverse proofs |32512|
| Unchanged retained inside inverse proofs |65792|
| Total off-original-node inside inverse proofs |98304|
| MPI1/4 source bound CSV |bitwise identical|
| Caller checkpoint versus original saved input |bitwise identical|

The real inverse uses actual residual and L to define an invariant box:
R=2*residual/(1-L), strictly inside the exact source triangle. Newton finds
a candidate; the real contraction and invariance prove its enclosed
preimage. There is no radius/guard scan. These inside inverse certificates
do not classify arbitrary outside points or certify complex Cauchy domains.

The full817152 point proof is reused through protected source/output/input
SHA and byte comparison of canonical/native model headers. Only the
previously missing32512 inverse queries are recalculated. No healthy
source's field integration is repeated to repair an unrelated global bound.

## Same target68526 classification

344 sources excluded by a whole-image power-plane witness,4 completely
inside,36 requiring a true curved source/target intersection. The old
problematic sources remain intersection candidates:

| Source | Cell | New whole real L |
|---|---:|---:|
|118|67942|.0184377680310|
|246|68526|.0186656755916|
|254|68526|.0194802288335|

They pass the new REAL inside inverse test. Their previous COMPLEX inverse
and density-branch failures are not marked solved or reused as successes.
No target packet has been published. All12 target moments, curved-domain
coverage, quadrature remainder and total physical error budget still need
the new canonical function API, not the old polynomial+v0 compiler.

## Provenance, timing and next work

Source/build/binary/input/output SHA preserved; fresh frozen runners and
new result directories refuse overwrite. Function full-point native run
158.67s at MPI4×OMP4; real-bound/missing-inverse stage38.48s at MPI4×OMP4,
141.22s at MPI1×OMP4, compilation excluded. These are Q128 reference tests
on a finite24-cell control, not full-grid Kepler evolution speed.

No authoritative repository source changes/commit/push. All handles60821,
42314,37875,39347 CLOSED; no new Slurm jobs or monitor loops remain.

Next: canonical edge value/first/second derivatives; clipped physical
source/target domains; actual whole target68526 all384 sources/all12 moments
with unchanged physical budget and transactional publish. The E=s^4r(t)
reference is C3 at its opposite vertex, not holomorphic there. Cartesian
Cauchy rectangles must not cross its s=0 singular chart. Source-corner
Duffy charts/geometrically determined partitions must preserve the SAME
stored function and complete domain coverage. No refit, basis swap, IC
reset or guard/precision search is authorized by these results.

Target integration must carry the immutable off-node thermal generator and
rho0/J versus EOS rho0/(J+delta), keep old U lower/budget, reject unproved
terms, and leave caller unchanged on any failure. Then full target set,
weights/cold reexpression/actual native RK4 and historical t6/8/10 plus
GIZMO longest-duration Kepler. Full objective remains ACTIVE, no blocker.
