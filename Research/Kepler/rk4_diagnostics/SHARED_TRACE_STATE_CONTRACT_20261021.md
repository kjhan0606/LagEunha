# Prospective shared-trace native state: implemented contract

2026-10-21. Driver implementation of
`ASTRA_SHARED_TRACE_FACTORY_PLAN_20261021.md`. Same24 physical cells and384
reference triangles, independent fresh representation control from the
original continuous IC. The rejected old packet was NOT refit or replaced.

## Mathematical reference and trace ownership

Reference vertices are the original binary64 coordinates, interpreted as
exact real inputs. Barycentric coordinates refer to the exact affine
reference geometry. Floating gradients/caches approximate evaluation of
that reference; they are not independent geometry/state parameters.

Initial v0(X,Y)=(-Y,X)*(X^2+Y^2+eps^2)^(-3/4), with unchanged eps=.047.
I3 uses the ten canonical rational barycentric P3 nodes, including the
four nodes j/3 of each oriented edge. It does not use acceleration, gravity
or a future orbit as enrichment. Original rho/P are retained.

For canonical edge endpoints A/B, r_e(t)=v0(A+t(B-A))-I3(v0)|e. Let
s=lambda_A+lambda_B and t=lambda_B/s. Fixed extension:

    E_e(r_e)=s^4 r_e(t), s>0; E_e=0 at the opposite vertex.

On its own edge s=1; on either other edge t=0 or1 and r_e=0. The homogeneous
degree4 extension is C3 at the opposite vertex for smooth r_e; H1 and
trace conformity do not require analytic continuation there. Define

    b_T=v0-I3(v0)-sum_e E_e(r_e).

It vanishes on every edge. Subtracting a multiple of local P3 interior
bubble N9 and applying an invertible two-component normalization preserves
that property. These basis transformations depend on original geometry/IC
mass metrics, not force residuals or evolved states. No old state is
projected, averaged or made to match new acceleration.

Scalar DOFs:

-237 shared vertices.
-1232 shared P3 edge-node values (two for each616 edges).
-384 local P3 interior values. Total shared P3 space1853.
-1176 canonical residual edge values:56 radial edges have rank1,
  560 nonradial edges rank2. Rank follows exact A cross B zero/nonzero,
  not a sampled singular-value cutoff.
-768 interior residual values, two per triangle.

Total3797 scalar DOFs. Previous broken space2088 and its nominal conforming
space1865 are retained as prior results. The new abstract space contains
the ideal component-compatible old P3+v0 residual space through the stated
edge/interior decomposition; independent edge/interior coefficients enlarge
it. This is not an exact embedding/refit of the old separately rounded
broken packet. Every dynamic coefficient remains free to evolve.

Each DOF has one global ID, parent/type and owner `global_id % MPI_ranks`.
The owner initializes/updates once; coefficient scatter broadcasts one
shared table. Each triangle stores only indices into that table. No
per-cell dynamic trace coefficient copies exist.

## Numerical evaluation and an important proof obligation

`material_shared_trace_batch.h` evaluates this fixed reference. Canonical
qnode caches approximate analytic IC/Taylor-remainder nodal evaluations.
They are **evaluation caches, not independently exact mathematical basis
coefficients**. The mathematical basis is the analytic specification above.
Any future target proof/compiler must enclose those analytic nodal constants
and cache/evaluator errors; it must NOT interpret each cached double as an
unrelated exact local trace function. That interval/compiler connection is
not implemented or certified by the present engineering checks.

Boundary x/v are evaluated in one canonical node/edge order with a common
evaluation origin. Bubble boundary value0 is its algebraically known trace.
The mathematical C0 argument is shared functions + shared DOFs, not observed
small jumps. Interior CPU/GPU evaluation differs by ordinary arithmetic;
no claim of exact floating continuity in every independent volume query.

Evaluation may subtract an origin for precision but never stores independently
rounded centered dynamic coefficients. Checkpoints store global coefficients
once, fixed reference/basis metadata, IDs/owners and epoch. On load, an
incomplete, duplicate or wrong-epoch table cannot overwrite caller state.
Basis caches are loaded/validated, and IC state is not regenerated on restart.

## Native operator

Position and velocity use the SAME space. Mass is integral N_a N_b dm.
Raw element Gram matrices are normalized algebraically and assembled by
sparse index scatter into102181 global CSR entries. Element18-or-smaller
Cholesky blocks supply an additive preconditioner. No global dense Schur.

The user's native `gfs_material12_stress_at_point` is applied in12+6 blocks;
thermal work and gravity power count once. Native CUDA evaluates the SAME
shared coefficient table, basis and original symmetric Piola stress. No
viscosity, damping, heat top-up, independent spin/kinetic reservoir, or
generator steering is added. Re=infinity.

Current finite state is ONLY h2e-4 control from newly solved acceleration,
not an accepted RK4 timestep. Physical thermal control remains adiabatic
P0*J^-gamma, with positive U checked independently of bulk kinetic energy.
Full local hybrid traction/resolved rotation-rate verification is still
pending because the maximum initial curvature gate failed first.

## Present gate and diagnostic

Initial IC/trace/checkpoint/native CPU/CUDA gates pass. Sparse equations and
global work/angular checks pass. Max chi=.003078319 exceeds unchanged.002;
every-cell mean max=.000109007 is below.001. No field/degree/mesh/order/
precision/tolerance scan follows this failure.

Saved-state witnesses isolate nearly all worst-point error to the nominal
uniform-pressure weak load. Its norm is7.74e-12 of full load; pressure-only
acceleration contributes chi=.003078319030, whereas gravity remainder max
chi=5.50001e-5. Direct normalized mass vs transformed raw Gram differs only
2.68e-17 in relative max entry and1.21e-17 in represented acceleration
action. Thus a raw-Gram transform discrepancy does not explain this failure.

Next operator decision must preserve exact constant-pressure nullspace AND
consistent physical thermal work/volume geometry. Subtracting a pressure
constant from force alone is not authorized by these results.
