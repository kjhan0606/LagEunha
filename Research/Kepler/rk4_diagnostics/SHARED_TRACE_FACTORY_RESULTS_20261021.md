# Shared-trace factory and native sparse RHS: initial pass, curvature fail

2026-10-21. Same original continuous IC /24 physical cells /384 source
triangles. New scalar space3797, no old state averaging/refit. See
`SHARED_TRACE_STATE_CONTRACT_20261021.md` for its exact functional reference,
ownership, cache/error proof obligation and relation to old spaces.

## Implemented and verified

`shared_trace_factory_v1_mpi{1,4}` both exit0. Actual CPU/OpenMP/CUDA
evaluation of817152 original positive nodes:

| Quantity | Result |
|---|---:|
| Max initial velocity error |8.8818e-16|
| Max cell K error / physical U |6.3378e-14|
| Max physical U relative error |4.6199e-17|
| Max L error /2K |4.3993e-16|
| Sampled minimum J |.9999999999999963|
| CPU/CUDA force relative difference |7.4067e-13|
| Shared edges / canonical trace queries |536 /9112|
| Canonical x/v jump, bubble trace, reversal difference |0 /0 /0|

Single-owner global state checkpoint roundtrip passes. Wrong epoch,
duplicate global ID and missing coefficient reject without changing caller.
MPI1/4 coefficient VALUES are identical, and canonical trace CSV is bitwise
identical. The coefficient binary files differ in49 signs of zero; no
nonzero numerical difference. Do not label those files bitwise identical.

Times1.21s (MPI1) /1.30s (MPI4), OMP4 each. These are finite24-cell controls,
not timings of a full128^2 evolution.

`shared_trace_rhs_v1_mpi4` implements sparse mass/RHS (102181 CSR entries),
element Cholesky additive preconditioner and original native Piola
stress-work. Both acceleration components converge in44 iterations;
actual equation relative residuals6.174e-17. No global dense Schur.

Owner-once h2e-4 finite CONTROL, with original physical thermal law:

-Sampled min J=.9999999959746803; no invalid nodes.
-Global energy rate residual /U=5.8588e-16.
-Global angular residual /2K=2.7360e-18.
-CPU/CUDA force difference7.4020e-13; point stress-work identity2.4592e-14.
-Canonical trace/bubble0; finite shared-state checkpoint/transaction passes.

This demonstrates initial representation and global native algebra. It does
not establish all cell-wise hybrid traction/rotation rates, a complete cold
transfer, an accepted timestep or long Kepler stability.

## Decisive unchanged curvature failure

Native RHS exit2. Max initial chi=.003078319055743 >.002.
Max of every physical cell's mean chi=.000109007057684 <.001.
The maximum criterion remains failed. Full bundle pass is false and this
state is not promoted to transfer/evolution. Physical conservation alone
does not waive the curvature gate.

`shared_trace_witness_v1` reproduces this from SAVED initial checkpoint and
saved acceleration, without a mass/RHS solve or coefficient update:

-Source5, cell17894, X=9.609458364476886, Y=-5.437500014486225.
-Radius~11.04, outside Rout8 in the low-density atmosphere.
-Smallest barycentric coordinate6.18078919e-7, near a reference vertex/edge.
-Acceleration error8.9739e-10; gradient error1.3754e-6.
-Divergence-of-acceleration error−1.372158149e-6.
-Velocity-gradient-square error−1.24745e-16.

The previous ideal-conforming fine P3 test on this same cell reported
chi9.8844e-8. This new worst failure is distinct evidence; it is not itself
a diagnosis of the historical inner-disk t6–10 collapse.

## Read-only / linear-load diagnosis, no parameter scan

`shared_trace_witness_v2` compares the old raw-Gram normalization against
direct integration of the SAME normalized basis at the SAME817152 nodes:

-Relative max entry difference2.6782e-17.
-Represented acceleration action difference1.2105e-17.
-Maximum cell kinetic difference /U1.2549e-16.

This rules out a significant normalized mass-transform mismatch as the
explanation. It is not a whole floating-error/condition-number certificate.

`shared_trace_witness_v3` adds a linear decomposition of the SAME initial
weak load on the saved shared state (not an alternative physical trajectory):

| Quantity | Result |
|---|---:|
| Pressure-only load norm /full load norm |7.7430e-12|
| Pressure-only acceleration max chi |.003078319030321|
| Gravity remainder max chi |5.50001e-5|
| Pressure fraction of divergence error at saved worst point |.9999999916507|

Thus the nearly uniform pressure weak-load residue accounts for essentially
all this maximum curvature failure. Gravity remainder meets the old bound.
The evidence points to constant-pressure cancellation and its amplification
in acceleration derivatives. It does not authorize arbitrary zeroing of
loads, pressure/thermal floors or altered physical stress.

No mesh, basis degree, quadrature order, precision, guard, U lower or
tolerance scan was performed. All frozen code/input/binary/output SHA
records and original result directories remain. No running script was edited.

## Next bundled decision — final corrected plan received

Existing Astra evaluated the bundle. Its first fixed-Pref proposal was:

    Fsigma = Q[(P-Pref)cof(F)grad(N)] + Pref*dV_boundary/dx
    Udot = -Q[(P-Pref)J div(v)] - Pref*Vdot_boundary.

Driver identified that fixed Pref alone only preserves its own constant
mode; arbitrary constant P leaves the old geometric quadrature defect.
Astra withdrew that claim and selected the corrected compatible discrete
volume/pressure operator in `ASTRA_GEOMETRIC_PRESSURE_PLAN_20261021.md`:

    nu_q=w_q J_q+(w_q/sum w)*(V_boundary-Q[J]),
    P_q=(gamma-1)U_q/nu_q,
    F=sum P_q*dnu_q/dx,
    Udot_q=-P_q*nudot_q.

This explicitly changes thermodynamic discrete volume while preserving
saved material mass/physical U and kinematic state. All nu must be positive,
and nu versus physical wJ / new EOS pressure versus original material
pressure must meet existing accuracy gates. It is not a heat reset or an
old-U overwrite. Actual derivatives of the SAME boundary-volume functional
must define reference forces/work/torque. The corrected plan is authorized
for driver implementation next; no constant-P operator code has been added
yet. The completed current fixture remains a curvature failure.

Original CUDA/MPI/OpenMP body and all older states remain preserved.
Real t6/8/10 and longest GIZMO-duration stability remain UNACHIEVED.
